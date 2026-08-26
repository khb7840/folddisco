// Save & Load the vector of file names
// Working with controller.path_vec: Vec<String>
// Lookup file format
// id\tpath\tinteger\tfloat
// id\tpath\tn_res\tplddt

use std::io::Write;
use std::fs::File;
use std::io::BufWriter;

use memmap2::Mmap;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use rayon::str::ParallelString;

use crate::utils::log::{log_msg, print_log_msg, FAIL, WARN};

// Binary cache of a parsed lookup file, saved next to it as <lookup>.cache.
// The record table has a fixed stride so the cache is decoded with rayon just
// like the mmap + par_lines text path; a serial reader would be slower than
// text parsing whenever more than a couple of threads are available.
// Layout (little-endian):
//   magic [u8; 8] | version u32 | count u64 | count * 40 byte records | names blob
// Record: id u64 | nres u64 | db_key u64 | plddt f32 | name_len u32 | name_offset u64
// name_offset is relative to the start of the names blob.
const LOOKUP_CACHE_MAGIC: &[u8; 8] = b"FDLOOKUP";
const LOOKUP_CACHE_VERSION: u32 = 1;
const LOOKUP_CACHE_HEADER_SIZE: usize = 20;
const LOOKUP_CACHE_RECORD_SIZE: usize = 40;

fn lookup_cache_path(path: &str) -> String {
    format!("{}.cache", path)
}

fn is_cache_fresh(lookup_path: &str, cache_path: &str) -> bool {
    let lookup_modified = std::fs::metadata(lookup_path).and_then(|meta| meta.modified());
    let cache_modified = std::fs::metadata(cache_path).and_then(|meta| meta.modified());
    match (lookup_modified, cache_modified) {
        (Ok(lookup_time), Ok(cache_time)) => cache_time >= lookup_time,
        _ => false,
    }
}

// Decode the cache in parallel. Returns None on anything unexpected (missing,
// truncated, bad magic or version, out-of-range name) so that the caller falls
// back to parsing the text lookup instead of returning partial data.
fn load_lookup_from_cache(path: &str) -> Option<Vec<(String, usize, usize, f32, usize)>> {
    let file = File::open(path).ok()?;
    let mmap = unsafe { Mmap::map(&file).ok()? };
    if mmap.len() < LOOKUP_CACHE_HEADER_SIZE || &mmap[..8] != LOOKUP_CACHE_MAGIC {
        return None;
    }
    if u32::from_le_bytes(mmap[8..12].try_into().unwrap()) != LOOKUP_CACHE_VERSION {
        return None;
    }
    let count = u64::from_le_bytes(mmap[12..20].try_into().unwrap()) as usize;
    let names_offset = count.checked_mul(LOOKUP_CACHE_RECORD_SIZE)?
        .checked_add(LOOKUP_CACHE_HEADER_SIZE)?;
    if mmap.len() < names_offset {
        return None;
    }
    let records = &mmap[LOOKUP_CACHE_HEADER_SIZE..names_offset];
    let names = &mmap[names_offset..];
    records.par_chunks_exact(LOOKUP_CACHE_RECORD_SIZE).map(|record| {
        let id = u64::from_le_bytes(record[0..8].try_into().unwrap()) as usize;
        let nres = u64::from_le_bytes(record[8..16].try_into().unwrap()) as usize;
        let db_key = u64::from_le_bytes(record[16..24].try_into().unwrap()) as usize;
        let plddt = f32::from_le_bytes(record[24..28].try_into().unwrap());
        let name_len = u32::from_le_bytes(record[28..32].try_into().unwrap()) as usize;
        let name_offset = u64::from_le_bytes(record[32..40].try_into().unwrap()) as usize;
        let name = names.get(name_offset..name_offset.checked_add(name_len)?)
            .and_then(|bytes| std::str::from_utf8(bytes).ok())?;
        Some((name.to_string(), id, nres, plddt, db_key))
    }).collect::<Option<Vec<_>>>()
}

fn save_lookup_cache(path: &str, lookup: &[(String, usize, usize, f32, usize)]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writer.write_all(LOOKUP_CACHE_MAGIC)?;
    writer.write_all(&LOOKUP_CACHE_VERSION.to_le_bytes())?;
    writer.write_all(&(lookup.len() as u64).to_le_bytes())?;
    let mut name_offset = 0u64;
    for (name, id, nres, plddt, db_key) in lookup {
        writer.write_all(&(*id as u64).to_le_bytes())?;
        writer.write_all(&(*nres as u64).to_le_bytes())?;
        writer.write_all(&(*db_key as u64).to_le_bytes())?;
        writer.write_all(&plddt.to_le_bytes())?;
        writer.write_all(&(name.len() as u32).to_le_bytes())?;
        writer.write_all(&name_offset.to_le_bytes())?;
        name_offset += name.len() as u64;
    }
    for (name, _, _, _, _) in lookup {
        writer.write_all(name.as_bytes())?;
    }
    writer.flush()
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

// pub fn load_lookup_from_file(path: &str) -> (Vec<String>, Vec<usize>, Vec<usize>, Vec<f32>) {
//     let mut path_vec: Vec<String> = Vec::new();
//     let mut numeric_id_vec: Vec<usize> = Vec::new();
//     let mut integer_vec: Vec<usize> = Vec::new();
//     let mut float_vec: Vec<f32> = Vec::new();
//     let file = std::fs::File::open(path).expect(&log_msg(FAIL, "Unable to open the lookup file"));
//     let reader = std::io::BufReader::new(file);
//     for line in reader.lines() {
//         let line = line.expect(&log_msg(FAIL, "Unable to read the lookup file"));
//         let line_vec: Vec<&str> = line.split("\t").collect();
//         numeric_id_vec.push(line_vec[0].parse::<usize>().unwrap());
//         path_vec.push(line_vec[1].to_string());
//         integer_vec.push(line_vec[2].parse::<usize>().unwrap());
//         float_vec.push(line_vec[3].parse::<f32>().unwrap());
//     }
//     (path_vec, numeric_id_vec, integer_vec, float_vec)
// }
pub fn load_lookup_from_file(path: &str) -> Vec<(String, usize, usize, f32, usize)> {
    let cache_path = lookup_cache_path(path);
    if is_cache_fresh(path, &cache_path) {
        if let Some(cached_lookup) = load_lookup_from_cache(&cache_path) {
            return cached_lookup;
        }
    }
    let file = std::fs::File::open(path).expect(&log_msg(FAIL, "Unable to open the lookup file"));
    let mmap = unsafe { Mmap::map(&file).expect(&log_msg(FAIL, "Unable to mmap the lookup file")) };
    let content = unsafe { std::str::from_utf8_unchecked(&mmap) };
    let loaded_lookup = content.par_lines().map(|line| {
        let mut split = line.split("\t");
        let id = split.next().unwrap().parse::<usize>().unwrap();
        let name = split.next().unwrap().to_string();
        let nres = split.next().unwrap().parse::<usize>().unwrap();
        let plddt = split.next().unwrap().parse::<f32>().unwrap();
        let db_key = split.next().unwrap_or(&id.to_string()).parse::<usize>().unwrap();
        (name, id, nres, plddt, db_key)
    }).collect::<Vec<_>>();
    // A read-only index directory is a normal deployment, so a failed cache
    // write must not fail the load.
    if let Err(err) = save_lookup_cache(&cache_path, &loaded_lookup) {
        print_log_msg(WARN, &format!("Unable to write the lookup cache {}: {}", &cache_path, err));
    }
    loaded_lookup
}


#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn test_save_and_load_lookup() {
        let path = "data/lookup_test.lookup";
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
        // Save the data to a file
        save_lookup_to_file(path, &path_vec, &numeric_id_vec, nres_vec.as_ref(), plddt_vec.as_ref(), numeric_db_key_vec.as_ref());

        // Load the data from the file
        let loaded_lookup = load_lookup_from_file(path);
        // Check that the loaded data is the same as the original data
        assert_eq!(loaded_lookup, expected_lookup);

        // Clean up the test file
        // std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_lookup_cache_roundtrip() {
        let (path, cache_path) = temp_lookup_path("cache_roundtrip");
        let names = (0..1000).map(|i| format!("dir{}/entry_{}.pdb", i % 7, i)).collect::<Vec<_>>();
        let ids = (0..1000).collect::<Vec<_>>();
        write_test_lookup(&path, &names, &ids);

        // First load parses the text file and writes the cache
        let from_text = load_lookup_from_file(&path);
        assert_eq!(from_text.len(), 1000);
        assert!(std::path::Path::new(&cache_path).is_file());

        // The cache decodes into exactly what the text path returned
        let from_cache = load_lookup_from_cache(&cache_path).expect("cache should decode");
        assert_eq!(from_cache, from_text);
        assert_eq!(load_lookup_from_file(&path), from_text);

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&cache_path);
    }

    #[test]
    fn test_lookup_cache_corrupted_falls_back_to_text() {
        let (path, cache_path) = temp_lookup_path("cache_corrupt");
        let names = vec!["a.pdb".to_string(), "b.pdb".to_string()];
        let ids = vec![0, 1];
        write_test_lookup(&path, &names, &ids);
        let expected = load_lookup_from_file(&path);

        // Truncated in the middle of the record table
        let cache = std::fs::OpenOptions::new().write(true).open(&cache_path).unwrap();
        cache.set_len(LOOKUP_CACHE_HEADER_SIZE as u64 + 10).unwrap();
        drop(cache);
        assert!(load_lookup_from_cache(&cache_path).is_none());
        assert_eq!(load_lookup_from_file(&path), expected);

        // Bad magic
        let mut cache = std::fs::OpenOptions::new().write(true).open(&cache_path).unwrap();
        cache.write_all(b"NOTACACH").unwrap();
        drop(cache);
        assert!(load_lookup_from_cache(&cache_path).is_none());
        assert_eq!(load_lookup_from_file(&path), expected);

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
        assert!(!is_cache_fresh(&path, &cache_path));
        assert_eq!(load_lookup_from_file(&path), vec![("b.pdb".to_string(), 1, 107, 51.0, 1001)]);

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&cache_path);
    }
}