//! Caches whose on-disk layout *is* the in-memory layout, so loading one is an
//! `mmap` plus a header check rather than a parse.
//!
//! The record table is a plain array of a `#[repr(C)]` POD type, laid out at an
//! 8-aligned offset, so it is handed to callers as a `&[R]` cast straight out of
//! the mapping. Nothing per entry is allocated, decoded or sorted at load time:
//! whatever ordering or auxiliary index the reader needs is computed once when
//! the cache is built. Variable-length names live in a single blob at the end
//! and are resolved only for the entries a caller actually asks about.
//!
//! ```text
//! magic [u8; 8] | version u32 | pad u32 | count u64 |
//! source length u64 | source mtime nanos u64 | names length u64 |   (48 bytes)
//! count * size_of::<R>() bytes of records                          (8-aligned)
//! names length bytes of UTF-8                                      (name blob)
//! ```
//!
//! ## What is validated, and why that is enough
//!
//! Every load checks the magic, the version, the length and mtime of the file
//! the cache was built from, and that the file size is exactly what the header's
//! three counts imply. The source length and mtime together identify the exact
//! input: an mtime comparison alone accepts a stale cache whenever the source is
//! rewritten inside one filesystem timestamp tick, which would serve wrong names
//! and keys with no signal at all. The size check pins the counts, so a
//! corrupted `count` cannot decode as a short-but-well-formed cache -- a zeroed
//! count would otherwise look like an empty database and make a search silently
//! find nothing. Anything that fails sends the caller back to parsing the
//! source, which is always still there.
//!
//! There is deliberately no checksum over the record table: computing one would
//! cost a full pass over the file, which is the entire thing this cache exists
//! to avoid, and after the size and identity checks no reachable write path
//! leaves a table that is corrupt yet passes. There is also deliberately no
//! temp-file-and-rename: cache content is a pure function of the source, so two
//! processes writing the same cache emit identical bytes from offset 0, and a
//! reader that catches a write in progress sees a short file and rejects it.
//! Both arguments hold only while the bytes stay deterministic -- putting a
//! timestamp or a thread id in the body would break them and bring atomicity
//! back into scope.
//!
//! ## What is not validated
//!
//! Names are checked for UTF-8 when they are resolved, not up front, because
//! validating a multi-gigabyte blob on every load would defeat the purpose. A
//! name that is not UTF-8 reads back as empty rather than panicking.
//!
//! Like every other mapping in this crate, the cache is mapped read-only and
//! assumed not to be rewritten underneath a live reader; the identity check
//! catches a source that moved on, not a cache actively being overwritten by
//! something that is not this code.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::mem::{align_of, size_of};

use memmap2::Mmap;

/// Header size, a multiple of 8 so that the record table that follows is
/// aligned for any of the record types below.
pub const HEADER_SIZE: usize = 48;

/// Marks a type as safe to reinterpret from arbitrary cache bytes.
///
/// # Safety
/// Implementors must be `#[repr(C)]`, contain no padding bytes, have an
/// alignment of at most 8, a size that is a multiple of 8 (so consecutive
/// records stay aligned), and be valid for every bit pattern -- i.e. be built
/// only out of integer and floating-point fields, with no enums, references,
/// `bool`s or `NonZero` types.
pub unsafe trait CacheRecord: Copy {
    /// Distinguishes cache kinds, so one kind's file is never read as another's.
    const MAGIC: &'static [u8; 8];
    const VERSION: u32;
}

/// Where a cache's bytes live: a mapping of the cache file, or an in-memory
/// image for when the cache could not be written to disk.
///
/// The in-memory image is held as `u64` words rather than bytes because the
/// record table is cast out of it: a `Vec<u8>` is only guaranteed 1-aligned, and
/// most allocators happen to return more, which is not something an alignment
/// invariant should rest on. `len` is the logical byte length, which is what the
/// header's size check is compared against; the word buffer rounds up past it.
enum Backing {
    Mapped(Mmap),
    Owned { words: Vec<u64>, len: usize },
}

impl Backing {
    #[inline]
    fn bytes(&self) -> &[u8] {
        match self {
            Backing::Mapped(mmap) => mmap,
            // SAFETY: `words` holds at least `len` bytes, and any bit pattern is
            // a valid `u8`.
            Backing::Owned { words, len } => unsafe {
                std::slice::from_raw_parts(words.as_ptr() as *const u8, *len)
            },
        }
    }

    /// Copy a serialised image into an 8-aligned word buffer.
    fn own(image: Vec<u8>) -> Self {
        let len = image.len();
        let mut words = vec![0u64; len.div_ceil(8)];
        // SAFETY: the destination holds `words.len() * 8 >= len` bytes.
        unsafe {
            std::ptr::copy_nonoverlapping(
                image.as_ptr(), words.as_mut_ptr() as *mut u8, len,
            );
        }
        Backing::Owned { words, len }
    }
}

/// Length and modification time of `path`, as stored in a cache header.
///
/// `None` when the file is gone or carries a timestamp outside the range this
/// encoding covers, in which case no cache is written or trusted.
pub fn source_identity(path: &str) -> Option<(u64, u64)> {
    let meta = std::fs::metadata(path).ok()?;
    let mtime = meta.modified().ok()?
        .duration_since(std::time::UNIX_EPOCH).ok()?
        .as_nanos();
    Some((meta.len(), u64::try_from(mtime).ok()?))
}

/// A loaded cache: a record table and a name blob, both borrowed from one
/// mapping.
pub struct PodCache<R: CacheRecord> {
    backing: Backing,
    count: usize,
    names_offset: usize,
    names_len: usize,
    _record: std::marker::PhantomData<R>,
}

impl<R: CacheRecord> PodCache<R> {
    /// Map `cache_path` and check it was built from the current `source_path`.
    ///
    /// Returns `None` for anything unexpected -- missing, truncated, wrong magic
    /// or version, a source that no longer matches, a size that disagrees with
    /// the header -- so the caller rebuilds from the source instead of returning
    /// stale or partial data.
    pub fn map(cache_path: &str, source_path: &str) -> Option<Self> {
        let file = File::open(cache_path).ok()?;
        let mmap = unsafe { Mmap::map(&file).ok()? };
        Self::from_backing(Backing::Mapped(mmap), source_path)
    }

    fn from_backing(backing: Backing, source_path: &str) -> Option<Self> {
        let bytes = backing.bytes();
        if bytes.len() < HEADER_SIZE || &bytes[..8] != R::MAGIC {
            return None;
        }
        if u32::from_le_bytes(bytes[8..12].try_into().unwrap()) != R::VERSION {
            return None;
        }
        let count = u64::from_le_bytes(bytes[16..24].try_into().unwrap()) as usize;
        let src_len = u64::from_le_bytes(bytes[24..32].try_into().unwrap());
        let src_mtime = u64::from_le_bytes(bytes[32..40].try_into().unwrap());
        let names_len = u64::from_le_bytes(bytes[40..48].try_into().unwrap()) as usize;
        if source_identity(source_path)? != (src_len, src_mtime) {
            return None;
        }

        let names_offset = count.checked_mul(size_of::<R>())?.checked_add(HEADER_SIZE)?;
        if bytes.len() != names_offset.checked_add(names_len)? {
            return None;
        }
        // Guaranteed by HEADER_SIZE and the size_of::<R>() contract, but a
        // misaligned cast would be undefined behaviour, so it is checked rather
        // than assumed.
        if bytes.as_ptr() as usize % align_of::<R>() != 0 || HEADER_SIZE % align_of::<R>() != 0 {
            return None;
        }

        Some(PodCache {
            backing, count, names_offset, names_len,
            _record: std::marker::PhantomData,
        })
    }

    /// Build a cache in memory rather than on disk, for when the cache
    /// directory is not writable.
    fn from_image(image: Vec<u8>, source_path: &str) -> Option<Self> {
        Self::from_backing(Backing::own(image), source_path)
    }

    /// A cache with no entries, for a reader with no database attached.
    pub fn empty() -> Self {
        PodCache {
            backing: Backing::own(vec![0u8; HEADER_SIZE]),
            count: 0,
            names_offset: HEADER_SIZE,
            names_len: 0,
            _record: std::marker::PhantomData,
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.count
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// The whole record table, cast out of the mapping.
    #[inline]
    pub fn records(&self) -> &[R] {
        if self.count == 0 {
            // `from_raw_parts` demands an aligned, non-null pointer even for a
            // zero length, and an empty cache has no record table to point at.
            return &[];
        }
        let bytes = self.backing.bytes();
        let start = HEADER_SIZE;
        // SAFETY: `from_backing` established that `bytes` is aligned for R, that
        // HEADER_SIZE is a multiple of R's alignment, and that the mapping holds
        // at least `count * size_of::<R>()` bytes from HEADER_SIZE onwards. R is
        // `CacheRecord`, whose contract is that every bit pattern is a valid
        // value. The mapping outlives the returned slice, which borrows `self`.
        unsafe {
            std::slice::from_raw_parts(
                bytes[start..].as_ptr() as *const R,
                self.count,
            )
        }
    }

    /// The raw name blob. Slices of it are handed out by the typed wrappers.
    #[inline]
    pub fn names_blob(&self) -> &[u8] {
        &self.backing.bytes()[self.names_offset..self.names_offset + self.names_len]
    }

    /// Resolve one name. Empty when the range is out of bounds or the bytes are
    /// not UTF-8; both mean a cache this loader should not have accepted, and
    /// neither is worth panicking a search over.
    #[inline]
    pub fn name_at(&self, offset: u64, len: u32) -> &str {
        let start = offset as usize;
        let end = start + len as usize;
        let blob = self.names_blob();
        if end > blob.len() {
            return "";
        }
        std::str::from_utf8(&blob[start..end]).unwrap_or("")
    }
}

/// Serialise a cache: `records` in the order readers will see them, and `names`
/// as the blob the records' offsets point into.
pub fn write_cache<R: CacheRecord, W: Write>(
    sink: &mut W, source_len: u64, source_mtime: u64,
    records: &[R], names: &[u8],
) -> std::io::Result<()> {
    sink.write_all(R::MAGIC)?;
    sink.write_all(&R::VERSION.to_le_bytes())?;
    sink.write_all(&0u32.to_le_bytes())?; // pad, keeps the counts 8-aligned
    sink.write_all(&(records.len() as u64).to_le_bytes())?;
    sink.write_all(&source_len.to_le_bytes())?;
    sink.write_all(&source_mtime.to_le_bytes())?;
    sink.write_all(&(names.len() as u64).to_le_bytes())?;
    // SAFETY: R is `CacheRecord`, so it is `#[repr(C)]` and padding-free; its
    // bytes are exactly its value and can be written as-is.
    let record_bytes = unsafe {
        std::slice::from_raw_parts(
            records.as_ptr() as *const u8,
            std::mem::size_of_val(records),
        )
    };
    sink.write_all(record_bytes)?;
    sink.write_all(names)?;
    sink.flush()
}

/// Write a cache next to its source and map it back, falling back to an
/// in-memory image when the directory is not writable.
///
/// A read-only index directory is a normal deployment, so a failed cache write
/// must not fail the load; it only costs the next process the parse again.
pub fn store_and_map<R: CacheRecord>(
    cache_path: &str, source_path: &str,
    records: &[R], names: &[u8],
) -> Option<PodCache<R>> {
    let (source_len, source_mtime) = source_identity(source_path)?;
    // A 4 MB buffer instead of the default 8 KB: these files reach several GB,
    // and the records go out in one `write_all` of many hundreds of MB.
    let written = File::create(cache_path).and_then(|file| {
        let mut writer = BufWriter::with_capacity(4 << 20, file);
        write_cache(&mut writer, source_len, source_mtime, records, names)
    });
    match written {
        Ok(()) => {
            if let Some(cache) = PodCache::map(cache_path, source_path) {
                return Some(cache);
            }
            // Written but not mappable: fall through to the in-memory image
            // rather than leaving the caller with nothing.
            crate::utils::log::print_log_msg(
                crate::utils::log::WARN,
                &format!("Wrote but could not map the cache {}; using memory", cache_path),
            );
        }
        Err(err) => {
            crate::utils::log::print_log_msg(
                crate::utils::log::WARN,
                &format!("Unable to write the cache {}: {}", cache_path, err),
            );
        }
    }
    let mut image = Vec::new();
    write_cache(&mut image, source_len, source_mtime, records, names).ok()?;
    PodCache::from_image(image, source_path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[repr(C)]
    #[derive(Copy, Clone, Debug, PartialEq)]
    struct TestRecord {
        key: u64,
        name_offset: u64,
        name_len: u32,
        _pad: u32,
    }
    unsafe impl CacheRecord for TestRecord {
        const MAGIC: &'static [u8; 8] = b"FDTESTC1";
        const VERSION: u32 = 1;
    }

    /// The alignment and padding promises `CacheRecord` makes.
    #[test]
    fn record_layout_is_castable() {
        assert_eq!(size_of::<TestRecord>() % 8, 0);
        assert!(align_of::<TestRecord>() <= 8);
        assert_eq!(HEADER_SIZE % 8, 0);
    }

    fn temp_paths(tag: &str) -> (String, String) {
        let unique = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)
            .unwrap().as_nanos();
        let dir = std::env::temp_dir().to_string_lossy().to_string();
        (format!("{}/fd_pod_{}_{}.src", dir, tag, unique),
         format!("{}/fd_pod_{}_{}.cache", dir, tag, unique))
    }

    fn sample() -> (Vec<TestRecord>, Vec<u8>) {
        let names: Vec<u8> = b"alpha beta gamma".to_vec();
        let records = vec![
            TestRecord { key: 7, name_offset: 0, name_len: 5, _pad: 0 },
            TestRecord { key: 9, name_offset: 6, name_len: 4, _pad: 0 },
            TestRecord { key: 11, name_offset: 11, name_len: 5, _pad: 0 },
        ];
        (records, names)
    }

    #[test]
    fn roundtrip_through_a_file() {
        let (src, cache) = temp_paths("roundtrip");
        std::fs::write(&src, b"source bytes").unwrap();
        let (records, names) = sample();

        let mapped = store_and_map(&cache, &src, &records, &names)
            .expect("cache should store and map");
        assert_eq!(mapped.len(), 3);
        assert_eq!(mapped.records(), records.as_slice());
        assert_eq!(mapped.name_at(0, 5), "alpha");
        assert_eq!(mapped.name_at(6, 4), "beta");
        assert_eq!(mapped.name_at(11, 5), "gamma");

        // And again from a cold map, i.e. the path a second process takes.
        let reopened = PodCache::<TestRecord>::map(&cache, &src).expect("should map");
        assert_eq!(reopened.records(), records.as_slice());

        let _ = std::fs::remove_file(&src);
        let _ = std::fs::remove_file(&cache);
    }

    #[test]
    fn a_source_that_changed_is_rejected() {
        let (src, cache) = temp_paths("stale");
        std::fs::write(&src, b"source bytes").unwrap();
        let (records, names) = sample();
        store_and_map(&cache, &src, &records, &names).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_some());

        // Same length, different content: only the mtime moves.
        std::fs::write(&src, b"SOURCE BYTES").unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        let _ = std::fs::remove_file(&src);
        let _ = std::fs::remove_file(&cache);
    }

    #[test]
    fn truncated_wrong_magic_and_wrong_version_are_rejected() {
        let (src, cache) = temp_paths("corrupt");
        std::fs::write(&src, b"source bytes").unwrap();
        let (records, names) = sample();
        store_and_map(&cache, &src, &records, &names).unwrap();
        let good = std::fs::read(&cache).unwrap();

        // Truncated: the size check fails.
        std::fs::write(&cache, &good[..good.len() - 3]).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        // Shorter than a header.
        std::fs::write(&cache, &good[..12]).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        // Wrong magic, i.e. another cache kind's file under this name.
        let mut wrong_magic = good.clone();
        wrong_magic[..8].copy_from_slice(b"FDOTHER1");
        std::fs::write(&cache, &wrong_magic).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        // Wrong version.
        let mut wrong_version = good.clone();
        wrong_version[8..12].copy_from_slice(&99u32.to_le_bytes());
        std::fs::write(&cache, &wrong_version).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        // A zeroed count must not decode as an empty database.
        let mut zero_count = good.clone();
        zero_count[16..24].copy_from_slice(&0u64.to_le_bytes());
        std::fs::write(&cache, &zero_count).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        // Nor may a wildly large one.
        let mut huge_count = good.clone();
        huge_count[16..24].copy_from_slice(&u64::MAX.to_le_bytes());
        std::fs::write(&cache, &huge_count).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        let _ = std::fs::remove_file(&src);
        let _ = std::fs::remove_file(&cache);
    }

    #[test]
    fn a_missing_cache_is_not_an_error() {
        let (src, cache) = temp_paths("missing");
        std::fs::write(&src, b"source bytes").unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());
        let _ = std::fs::remove_file(&src);
    }

    #[test]
    fn names_are_bounds_and_utf8_checked_at_resolution() {
        let (src, cache) = temp_paths("names");
        std::fs::write(&src, b"source bytes").unwrap();
        let records = vec![TestRecord { key: 1, name_offset: 0, name_len: 2, _pad: 0 }];
        let mapped = store_and_map(&cache, &src, &records, &[0xff, 0xfe]).unwrap();
        assert_eq!(mapped.name_at(0, 2), "");   // not UTF-8
        assert_eq!(mapped.name_at(0, 99), "");  // past the blob
        let _ = std::fs::remove_file(&src);
        let _ = std::fs::remove_file(&cache);
    }

    #[test]
    fn an_empty_cache_has_no_records() {
        let empty = PodCache::<TestRecord>::empty();
        assert!(empty.is_empty());
        assert_eq!(empty.len(), 0);
        assert_eq!(empty.records(), &[]);
        assert_eq!(empty.names_blob(), &[] as &[u8]);
        assert_eq!(empty.name_at(0, 3), "");
    }

    /// An unwritable directory must degrade to an in-memory image, not fail.
    /// The image has to be 8-aligned for the record cast to be sound.
    #[test]
    fn an_unwritable_cache_path_falls_back_to_memory() {
        let (src, _) = temp_paths("readonly");
        std::fs::write(&src, b"source bytes").unwrap();
        let (records, names) = sample();
        let cache = "/proc/definitely-not-writable/fd.cache";
        let mapped = store_and_map(cache, &src, &records, &names)
            .expect("should fall back to an in-memory image");
        assert_eq!(mapped.records(), records.as_slice());
        assert_eq!(mapped.name_at(0, 5), "alpha");
        assert_eq!(mapped.records().as_ptr() as usize % align_of::<TestRecord>(), 0);
        let _ = std::fs::remove_file(&src);
    }
}
