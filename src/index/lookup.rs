// Save & Load the vector of file names
// Working with controller.path_vec: Vec<String>
// Lookup file format
// id\tpath\tinteger\tfloat
// id\tpath\tn_res\tplddt

use std::io::Write;
use std::fs::File;
use std::io::BufWriter;

use memmap2::Mmap;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rayon::str::ParallelString;

use crate::utils::log::{log_msg, FAIL};
use crate::utils::pod_cache::{store_and_map, CacheRecord, PodCache};

/// One lookup entry, in the fixed layout the cache is mapped as.
///
/// `#[repr(C)]`, 40 bytes, no padding: the field order puts the four 8-byte
/// fields first and pairs the three 4-byte ones, so the bytes on disk are
/// exactly the value and the table can be cast straight out of the mapping.
/// `name_offset` is relative to the start of the name blob.
#[repr(C)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LookupRecord {
    pub id: u64,
    pub nres: u64,
    pub db_key: u64,
    pub name_offset: u64,
    pub plddt: f32,
    pub name_len: u32,
}

// SAFETY: `#[repr(C)]`, four u64s followed by an f32 and a u32, so 40 bytes with
// no padding and 8-byte alignment; every field is an integer or a float, so
// every bit pattern is a valid value.
unsafe impl CacheRecord for LookupRecord {
    const MAGIC: &'static [u8; 8] = b"FDLOOKUP";
    // v1/v2 were the parse-on-load text and record formats; v3 is mapped and
    // cast, with a padded header so the record table is 8-aligned.
    const VERSION: u32 = 3;
}

/// One lookup entry as its consumers want it.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LookupEntry<'a> {
    pub name: &'a str,
    pub id: usize,
    pub nres: usize,
    pub plddt: f32,
    pub db_key: usize,
}

/// The lookup table of an index, held as a mapping rather than as parsed data.
///
/// Loading one is an `mmap`, a header check and a scan of the record table --
/// no per-entry allocation, no text parsing, no sorting. Names are `&str`
/// slices of the mapping, resolved when a caller asks for them, so the entries
/// a query never reports are never touched.
///
/// The previous representation was a `Vec<(String, usize, usize, f32, usize)>`,
/// which on a 53.7 M-entry lookup meant 53.7 M simultaneously live `String`s:
/// 8.1 GB of live heap for 1.3 GB of names, and roughly three quarters of the
/// load time spent in the allocator rather than on the data.
pub struct LookupTable {
    cache: PodCache<LookupRecord>,
}

impl LookupTable {
    /// Map a cache and check it describes the current lookup file.
    ///
    /// Returns `None` for anything the mapping layer rejects, and additionally
    /// for a record table that is not structurally sound. That second check
    /// costs one parallel pass over the records -- tens of milliseconds even at
    /// 50 M entries -- and is not optional: a torn write leaves all-zero
    /// records, whose every field is individually in range, and an `id` of 0
    /// would then quietly attribute those hits to the first structure in the
    /// index. Anything rejected here falls back to parsing the text.
    fn open(cache_path: &str, lookup_path: &str) -> Option<Self> {
        let cache = PodCache::<LookupRecord>::map(cache_path, lookup_path)?;
        let count = cache.len();
        let names_len = cache.names_blob().len() as u64;
        let sound = cache.records().par_iter().all(|record| {
            // An empty name is not data: a lookup entry is a file path. And an
            // id at or past the entry count is used to index a vector of
            // exactly that length.
            record.name_len > 0
                && record.id < count as u64
                && record.name_offset.saturating_add(record.name_len as u64) <= names_len
        });
        if !sound {
            return None;
        }
        Some(LookupTable { cache })
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.cache.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.cache.is_empty()
    }

    #[inline]
    pub fn records(&self) -> &[LookupRecord] {
        self.cache.records()
    }

    /// The name of entry `index`, borrowed from the mapping.
    #[inline]
    pub fn name(&self, index: usize) -> &str {
        let record = &self.records()[index];
        self.cache.name_at(record.name_offset, record.name_len)
    }

    /// Entry `index`, with its name resolved.
    #[inline]
    pub fn entry(&self, index: usize) -> LookupEntry<'_> {
        let record = &self.records()[index];
        LookupEntry {
            name: self.cache.name_at(record.name_offset, record.name_len),
            id: record.id as usize,
            nres: record.nres as usize,
            plddt: record.plddt,
            db_key: record.db_key as usize,
        }
    }

    /// Every name, in entry order.
    pub fn names(&self) -> impl Iterator<Item = &str> + '_ {
        (0..self.len()).map(move |i| self.name(i))
    }

    /// The owned tuples the previous representation handed out. For callers
    /// that genuinely need to own the data, and for tests; a query path should
    /// use `entry` and leave the names in the mapping.
    pub fn to_owned_vec(&self) -> Vec<(String, usize, usize, f32, usize)> {
        (0..self.len()).map(|i| {
            let e = self.entry(i);
            (e.name.to_string(), e.id, e.nres, e.plddt, e.db_key)
        }).collect()
    }
}

fn lookup_cache_path(path: &str) -> String {
    format!("{}.cache", path)
}

/// Parse the text lookup into cache records, with the names copied into one
/// contiguous blob.
///
/// The parse itself is parallel; the blob is assembled in a second serial pass
/// because each name's offset depends on the lengths of all the names before
/// it. Neither pass allocates per entry.
fn parse_lookup_text(content: &str) -> (Vec<LookupRecord>, Vec<u8>) {
    // (id, nres, plddt, db_key, name offset within `content`, name length)
    let parsed: Vec<(u64, u64, f32, u64, u64, u32)> = content.par_lines().map(|line| {
        let mut split = line.split('\t');
        let id = split.next().unwrap().parse::<u64>().unwrap();
        let name = split.next().unwrap();
        let nres = split.next().unwrap().parse::<u64>().unwrap();
        let plddt = split.next().unwrap().parse::<f32>().unwrap();
        // The fifth column was added later; an older lookup reuses the id.
        let db_key = match split.next() {
            Some(text) => text.parse::<u64>().unwrap(),
            None => id,
        };
        let name_offset = name.as_ptr() as usize - content.as_ptr() as usize;
        (id, nres, plddt, db_key, name_offset as u64, name.len() as u32)
    }).collect();

    let names_len: usize = parsed.iter().map(|p| p.5 as usize).sum();
    let mut names = Vec::with_capacity(names_len);
    let mut records = Vec::with_capacity(parsed.len());
    let source = content.as_bytes();
    for (id, nres, plddt, db_key, text_offset, name_len) in parsed {
        let start = text_offset as usize;
        let name_offset = names.len() as u64;
        names.extend_from_slice(&source[start..start + name_len as usize]);
        records.push(LookupRecord {
            id, nres, db_key, name_offset, plddt, name_len,
        });
    }
    (records, names)
}

pub fn save_lookup_to_file(
    path: &str, path_vec: &Vec<String>, numeric_id_vec: &Vec<usize>, 
    optional_int_vec: Option<&Vec<usize>>, optional_float_vec: Option<&Vec<f32>>,
    numeric_db_key_vec: Option<&Vec<usize>>,
) {
    assert_eq!(path_vec.len(), numeric_id_vec.len());
    if optional_int_vec.is_some() {
        assert_eq!(path_vec.len(), optional_int_vec.unwrap().len());
    }
    if optional_float_vec.is_some() {
        assert_eq!(path_vec.len(), optional_float_vec.unwrap().len());
    }
    if numeric_db_key_vec.is_some() {
        assert_eq!(path_vec.len(), numeric_db_key_vec.unwrap().len());
    }
    
    // Save the vector of file names to a file
    let mut file = BufWriter::new(File::create(path).expect(&log_msg(FAIL, "Unable to create the lookup file")));
    for i in 0..path_vec.len() {
        let mut numeric_db_key = numeric_id_vec[i];
        if numeric_db_key_vec.is_some() {
            // Set numeric_db_key as db_key_vec
            numeric_db_key = numeric_db_key_vec.unwrap()[i];
        }

        let line = match (optional_int_vec, optional_float_vec) {
            (Some(int_vec), Some(float_vec)) => {
                format!("{}\t{}\t{}\t{}\t{}\n", numeric_id_vec[i], path_vec[i], int_vec[i], float_vec[i], numeric_db_key)
            },
            (Some(int_vec), None) => {
                format!("{}\t{}\t{}\t{}\t{}\n", numeric_id_vec[i], path_vec[i], int_vec[i], 0.0, numeric_db_key)
            },
            (None, Some(float_vec)) => {
                format!("{}\t{}\t{}\t{}\t{}\n", numeric_id_vec[i], path_vec[i], 0, float_vec[i], numeric_db_key)
            },
            (None, None) => {
                format!("{}\t{}\t{}\t{}\t{}\n", numeric_id_vec[i], path_vec[i], 0, 0.0, numeric_db_key)
            }
        };
        file.write_all(line.as_bytes()).expect(&log_msg(FAIL, "Unable to write the lookup file"));
    }
}

/// Load the lookup table of an index.
///
/// A valid cache next to the lookup file is mapped and used as-is. Otherwise
/// the text is parsed once, the cache is written, and the result is served from
/// that same mapping -- so every reader takes exactly one code path and there
/// is no second representation to keep in step.
pub fn load_lookup_from_file(path: &str) -> LookupTable {
    let cache_path = lookup_cache_path(path);
    if let Some(table) = LookupTable::open(&cache_path, path) {
        return table;
    }
    let file = std::fs::File::open(path).expect(&log_msg(FAIL, "Unable to open the lookup file"));
    let mmap = unsafe {
        Mmap::map(&file).expect(&log_msg(FAIL, "Unable to mmap the lookup file"))
    };
    let content = unsafe { std::str::from_utf8_unchecked(&mmap) };
    let (records, names) = parse_lookup_text(content);
    let cache = store_and_map(&cache_path, path, &records, &names).expect(
        &log_msg(FAIL, "Unable to build the lookup cache")
    );
    let table = LookupTable { cache };
    // The freshly built cache has to pass the same checks a mapped one does; if
    // it does not, the parse produced something this loader would refuse, and
    // saying so beats serving it.
    debug_assert!(
        table.records().iter().all(|r| r.name_len > 0 && r.id < table.len() as u64),
        "the lookup text produced records the cache loader would reject"
    );
    table
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::pod_cache::HEADER_SIZE;
    use std::mem::size_of;

    /// Offsets of the header fields this module's tests patch.
    const COUNT_OFFSET: u64 = 16;

    // Unique lookup path in a temporary directory, so that tests running in
    // parallel never share a cache.
    fn temp_lookup_path(tag: &str) -> (String, String) {
        let unique = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)
            .unwrap().as_nanos();
        let path = format!(
            "{}/folddisco_{}_{}.lookup", std::env::temp_dir().to_string_lossy(), tag, unique
        );
        let cache_path = lookup_cache_path(&path);
        (path, cache_path)
    }

    fn write_test_lookup(path: &str, names: &Vec<String>, ids: &Vec<usize>) {
        let nres = ids.iter().map(|i| 100 + i * 7).collect::<Vec<_>>();
        let plddt = ids.iter().map(|i| 50.0 + *i as f32).collect::<Vec<_>>();
        let db_key = ids.iter().map(|i| 1000 + i).collect::<Vec<_>>();
        save_lookup_to_file(path, names, ids, Some(&nres), Some(&plddt), Some(&db_key));
    }

    /// Overwrite bytes of the cache in place.
    /// Only working within tests and should not be used in production code.
    fn patch_cache(cache_path: &str, offset: u64, bytes: &[u8]) {
        use std::io::{Seek, SeekFrom};
        let mut cache = std::fs::OpenOptions::new().write(true).open(cache_path).unwrap();
        cache.seek(SeekFrom::Start(offset)).unwrap();
        cache.write_all(bytes).unwrap();
    }

    /// The layout promises `LookupRecord` makes to `CacheRecord`.
    #[test]
    fn test_record_layout() {
        assert_eq!(size_of::<LookupRecord>(), 40);
        assert_eq!(size_of::<LookupRecord>() % 8, 0);
        assert_eq!(HEADER_SIZE % 8, 0);
    }

    #[test]
    fn test_save_and_load_lookup() {
        let (path, cache_path) = temp_lookup_path("save_and_load");
        let path_vec = vec!["path1.pdb".to_string(), "path2.pdb".to_string(), "path3.pdb".to_string()];
        let numeric_id_vec = vec![0, 1, 2];
        let nres_vec = Some(vec![100, 200, 5000]);
        let plddt_vec = Some(vec![50.0, 60.0, 70.0]);
        let numeric_db_key_vec = Some(vec![100, 110, 200]);
        let expected_lookup = vec![
            ("path1.pdb".to_string(), 0, 100, 50.0, 100),
            ("path2.pdb".to_string(), 1, 200, 60.0, 110),
            ("path3.pdb".to_string(), 2, 5000, 70.0, 200)
        ];
        save_lookup_to_file(&path, &path_vec, &numeric_id_vec, nres_vec.as_ref(), plddt_vec.as_ref(), numeric_db_key_vec.as_ref());

        let loaded = load_lookup_from_file(&path);
        assert_eq!(loaded.to_owned_vec(), expected_lookup);
        // Field accessors agree with the tuples.
        assert_eq!(loaded.name(2), "path3.pdb");
        assert_eq!(loaded.entry(1).db_key, 110);
        assert_eq!(loaded.names().collect::<Vec<_>>(), vec!["path1.pdb", "path2.pdb", "path3.pdb"]);

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&cache_path);
    }

    #[test]
    fn test_lookup_cache_roundtrip() {
        let (path, cache_path) = temp_lookup_path("cache_roundtrip");
        let names = (0..1000).map(|i| format!("dir{}/entry_{}.pdb", i % 7, i)).collect::<Vec<_>>();
        let ids = (0..1000).collect::<Vec<_>>();
        write_test_lookup(&path, &names, &ids);

        // First load parses the text and writes the cache
        let from_text = load_lookup_from_file(&path).to_owned_vec();
        assert_eq!(from_text.len(), 1000);
        assert!(std::path::Path::new(&cache_path).is_file());

        // The mapped cache decodes into exactly what the text path returned
        let from_cache = LookupTable::open(&cache_path, &path).expect("cache should map");
        assert_eq!(from_cache.to_owned_vec(), from_text);
        assert_eq!(load_lookup_from_file(&path).to_owned_vec(), from_text);

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&cache_path);
    }

    /// A lookup without the fifth column, as older indices have.
    #[test]
    fn test_lookup_without_db_key_column_reuses_the_id() {
        let (path, cache_path) = temp_lookup_path("four_column");
        std::fs::write(&path, "0\ta.pdb\t100\t50.0\n7\tb.pdb\t200\t60.0\n").unwrap();
        let loaded = load_lookup_from_file(&path);
        assert_eq!(loaded.entry(0).db_key, 0);
        assert_eq!(loaded.entry(1).db_key, 7);
        assert_eq!(loaded.entry(1).id, 7);
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&cache_path);
    }

    #[test]
    fn test_lookup_cache_corrupted_falls_back_to_text() {
        let (path, cache_path) = temp_lookup_path("cache_corrupt");
        let names = vec!["a.pdb".to_string(), "b.pdb".to_string()];
        let ids = vec![0, 1];
        write_test_lookup(&path, &names, &ids);
        let expected = load_lookup_from_file(&path).to_owned_vec();

        // Truncated in the middle of the record table
        let cache = std::fs::OpenOptions::new().write(true).open(&cache_path).unwrap();
        cache.set_len(HEADER_SIZE as u64 + 10).unwrap();
        drop(cache);
        assert!(LookupTable::open(&cache_path, &path).is_none());
        assert_eq!(load_lookup_from_file(&path).to_owned_vec(), expected);

        // Bad magic
        let mut cache = std::fs::OpenOptions::new().write(true).open(&cache_path).unwrap();
        cache.write_all(b"NOTACACH").unwrap();
        drop(cache);
        assert!(LookupTable::open(&cache_path, &path).is_none());
        assert_eq!(load_lookup_from_file(&path).to_owned_vec(), expected);

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&cache_path);
    }

    #[test]
    fn test_lookup_rewritten_within_one_mtime_tick_is_not_served_from_cache() {
        // A lookup rewritten inside one filesystem timestamp tick (1 s on ext3, HFS+,
        // most NFS servers) leaves the cache with an mtime that is not older. The
        // cache then has to be rejected on the file's identity, not on its clock:
        // serving it would report every hit under the wrong name and pLDDT.
        let (path, cache_path) = temp_lookup_path("cache_sametick");
        write_test_lookup(&path, &vec!["a.pdb".to_string()], &vec![0]);
        load_lookup_from_file(&path);
        let cache_modified = std::fs::metadata(&cache_path).unwrap().modified().unwrap();

        write_test_lookup(&path, &vec!["totally_different.pdb".to_string()], &vec![1]);
        let lookup = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
        lookup.set_times(std::fs::FileTimes::new().set_modified(cache_modified)).unwrap();
        drop(lookup);
        assert_eq!(
            std::fs::metadata(&path).unwrap().modified().unwrap(), cache_modified,
            "the test needs both mtimes equal to mean anything"
        );

        assert!(LookupTable::open(&cache_path, &path).is_none());
        assert_eq!(
            load_lookup_from_file(&path).to_owned_vec(),
            vec![("totally_different.pdb".to_string(), 1, 107, 51.0, 1001)]
        );

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&cache_path);
    }

    #[test]
    fn test_lookup_cache_with_zeroed_count_is_rejected() {
        // A zeroed count passes magic, version and the minimum length check, and
        // would decode as an empty index: the search then silently finds nothing.
        let (path, cache_path) = temp_lookup_path("cache_count0");
        let names = (0..5).map(|i| format!("prot_{}.pdb", i)).collect::<Vec<_>>();
        write_test_lookup(&path, &names, &(0..5).collect());
        let expected = load_lookup_from_file(&path).to_owned_vec();
        assert_eq!(expected.len(), 5);

        patch_cache(&cache_path, COUNT_OFFSET, &0u64.to_le_bytes());
        assert!(LookupTable::open(&cache_path, &path).is_none());
        assert_eq!(load_lookup_from_file(&path).to_owned_vec(), expected);

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&cache_path);
    }

    /// This is why the load does a pass over the records rather than only
    /// casting them: every field of an all-zero record is individually in
    /// range, and an `id` of 0 would attribute those hits to entry 0.
    #[test]
    fn test_lookup_cache_with_a_zeroed_record_is_rejected() {
        let (path, cache_path) = temp_lookup_path("cache_zerohole");
        let names = (0..5).map(|i| format!("prot_{}.pdb", i)).collect::<Vec<_>>();
        write_test_lookup(&path, &names, &(0..5).collect());
        let expected = load_lookup_from_file(&path).to_owned_vec();

        let record_2 = HEADER_SIZE + 2 * size_of::<LookupRecord>();
        patch_cache(&cache_path, record_2 as u64, &[0u8; size_of::<LookupRecord>()]);
        assert!(LookupTable::open(&cache_path, &path).is_none());
        assert_eq!(load_lookup_from_file(&path).to_owned_vec(), expected);

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&cache_path);
    }

    /// `id` is used to index a vector of exactly `len()` entries, so one past
    /// the end has to be caught at load rather than panicking mid-search.
    #[test]
    fn test_lookup_cache_with_an_out_of_range_id_is_rejected() {
        let (path, cache_path) = temp_lookup_path("cache_badid");
        let names = (0..5).map(|i| format!("prot_{}.pdb", i)).collect::<Vec<_>>();
        write_test_lookup(&path, &names, &(0..5).collect());
        let expected = load_lookup_from_file(&path).to_owned_vec();

        // Record 3's id field is the first 8 bytes of the record.
        let record_3_id = HEADER_SIZE + 3 * size_of::<LookupRecord>();
        patch_cache(&cache_path, record_3_id as u64, &99u64.to_le_bytes());
        assert!(LookupTable::open(&cache_path, &path).is_none());
        assert_eq!(load_lookup_from_file(&path).to_owned_vec(), expected);

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&cache_path);
    }

    #[test]
    fn test_lookup_cache_older_than_lookup_is_ignored() {
        let (path, cache_path) = temp_lookup_path("cache_stale");
        write_test_lookup(&path, &vec!["a.pdb".to_string()], &vec![0]);
        load_lookup_from_file(&path);
        assert!(std::path::Path::new(&cache_path).is_file());

        // Rewrite the lookup and backdate the cache, so the cache holds the
        // previous content and must be discarded
        write_test_lookup(&path, &vec!["b.pdb".to_string()], &vec![1]);
        let lookup_modified = std::fs::metadata(&path).unwrap().modified().unwrap();
        let cache = std::fs::OpenOptions::new().write(true).open(&cache_path).unwrap();
        cache.set_times(std::fs::FileTimes::new().set_modified(
            lookup_modified - std::time::Duration::from_secs(10)
        )).unwrap();
        drop(cache);
        assert!(LookupTable::open(&cache_path, &path).is_none());
        assert_eq!(
            load_lookup_from_file(&path).to_owned_vec(),
            vec![("b.pdb".to_string(), 1, 107, 51.0, 1001)]
        );

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&cache_path);
    }
}
