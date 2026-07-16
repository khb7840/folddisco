// Save & Load the vector of file names
// Working with controller.path_vec: Vec<String>
// Lookup file format
// id\tpath\tinteger\tfloat
// id\tpath\tn_res\tplddt

use std::io::{BufReader, Error, ErrorKind, Read, Write};
use std::fs::File;
use std::io::BufWriter;

use memmap2::Mmap;
use rayon::iter::ParallelIterator;
use rayon::str::ParallelString;

use crate::utils::log::{log_msg, FAIL};

const LOOKUP_CACHE_MAGIC: &[u8; 8] = b"FDLKP001";

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

fn read_u32<R: Read>(reader: &mut R) -> Result<u32, Error> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_u64<R: Read>(reader: &mut R) -> Result<u64, Error> {
    let mut bytes = [0u8; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

fn read_f32<R: Read>(reader: &mut R) -> Result<f32, Error> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(f32::from_le_bytes(bytes))
}

fn save_lookup_cache(path: &str, lookup: &[(String, usize, usize, f32, usize)]) -> Result<(), Error> {
    let mut writer = BufWriter::new(File::create(path)?);
    writer.write_all(LOOKUP_CACHE_MAGIC)?;
    writer.write_all(&(lookup.len() as u64).to_le_bytes())?;

    for (name, id, nres, plddt, db_key) in lookup {
        let name_bytes = name.as_bytes();
        writer.write_all(&(name_bytes.len() as u32).to_le_bytes())?;
        writer.write_all(name_bytes)?;
        writer.write_all(&(*id as u64).to_le_bytes())?;
        writer.write_all(&(*nres as u64).to_le_bytes())?;
        writer.write_all(&plddt.to_le_bytes())?;
        writer.write_all(&(*db_key as u64).to_le_bytes())?;
    }
    writer.flush()
}

fn load_lookup_from_cache(path: &str) -> Result<Vec<(String, usize, usize, f32, usize)>, Error> {
    let mut reader = BufReader::new(File::open(path)?);

    let mut magic = [0u8; 8];
    reader.read_exact(&mut magic)?;
    if &magic != LOOKUP_CACHE_MAGIC {
        return Err(Error::new(ErrorKind::InvalidData, "Invalid lookup cache header"));
    }
    let entry_count = read_u64(&mut reader)? as usize;
    let mut lookup = Vec::with_capacity(entry_count);

    for _ in 0..entry_count {
        let name_len = read_u32(&mut reader)? as usize;
        let mut name_bytes = vec![0u8; name_len];
        reader.read_exact(&mut name_bytes)?;
        let name = String::from_utf8(name_bytes)
            .map_err(|_| Error::new(ErrorKind::InvalidData, "Invalid UTF-8 in lookup cache"))?;
        let id = read_u64(&mut reader)? as usize;
        let nres = read_u64(&mut reader)? as usize;
        let plddt = read_f32(&mut reader)?;
        let db_key = read_u64(&mut reader)? as usize;
        lookup.push((name, id, nres, plddt, db_key));
    }

    Ok(lookup)
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
        match load_lookup_from_cache(&cache_path) {
            Ok(cached_lookup) => return cached_lookup,
            Err(err) => eprintln!(
                "{}",
                log_msg(
                    FAIL,
                    &format!("Lookup cache load failed (fallback to text parsing): {}", err)
                )
            ),
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

    if let Err(err) = save_lookup_cache(&cache_path, &loaded_lookup) {
        eprintln!(
            "{}",
            log_msg(
                FAIL,
                &format!("Failed to write lookup cache {}: {}", cache_path, err)
            )
        );
    }
    loaded_lookup
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};
    use std::thread::sleep;
    use std::time::Duration;

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
    fn test_lookup_cache_is_created_and_refreshed() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = format!(
            "{}/folddisco_lookup_cache_{}.lookup",
            std::env::temp_dir().to_string_lossy(),
            unique
        );
        let cache_path = lookup_cache_path(&path);

        let path_vec_1 = vec!["a.pdb".to_string()];
        let numeric_id_vec_1 = vec![0];
        let nres_vec_1 = Some(vec![10]);
        let plddt_vec_1 = Some(vec![50.0]);
        let db_key_vec_1 = Some(vec![100]);
        save_lookup_to_file(
            &path,
            &path_vec_1,
            &numeric_id_vec_1,
            nres_vec_1.as_ref(),
            plddt_vec_1.as_ref(),
            db_key_vec_1.as_ref(),
        );
        let loaded_1 = load_lookup_from_file(&path);
        assert_eq!(loaded_1, vec![("a.pdb".to_string(), 0, 10, 50.0, 100)]);
        assert!(std::path::Path::new(&cache_path).is_file());

        sleep(Duration::from_millis(5));
        let path_vec_2 = vec!["b.pdb".to_string()];
        let numeric_id_vec_2 = vec![1];
        let nres_vec_2 = Some(vec![20]);
        let plddt_vec_2 = Some(vec![80.0]);
        let db_key_vec_2 = Some(vec![200]);
        save_lookup_to_file(
            &path,
            &path_vec_2,
            &numeric_id_vec_2,
            nres_vec_2.as_ref(),
            plddt_vec_2.as_ref(),
            db_key_vec_2.as_ref(),
        );
        let loaded_2 = load_lookup_from_file(&path);
        assert_eq!(loaded_2, vec![("b.pdb".to_string(), 1, 20, 80.0, 200)]);

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&cache_path);
    }
}