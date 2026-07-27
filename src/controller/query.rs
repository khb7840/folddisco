// File: query.rs
// Created: 2023-12-22 17:00:50
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use rustc_hash::FxHashMap as HashMap;
use crate::geometry::core::{GeometricHash, HashType};
use crate::index::indextable::FolddiscoIndex;
use crate::utils::convert::{is_aa_group_char, map_one_letter_to_u8_vec};
use crate::utils::combination::CombinationIterator;
use crate::utils::log::{log_msg, print_log_msg, FAIL, WARN};
use super::expand::{FeatureExpander, ToleranceConfig};
use super::feature::get_single_feature;
use super::io::read_compact_structure;
use crate::structure::core::CompactStructure;

/// Ceiling on the hashes generated for a single residue pair.
///
/// Amino acid substitutions multiply with the geometric neighbourhood, so a
/// wildcard on both residues of a pair (`400` alternatives) together with a wide
/// expansion radius could generate tens of thousands of hashes for one pair. The
/// neighbourhood is enumerated nearest-first, so cutting it off here drops the
/// least useful variants.
const MAX_HASHES_PER_PAIR: usize = 4096;

/// Calculate IDF (Inverse Document Frequency) for a given GeometricHash
/// Uses index to determine hash count
pub fn calculate_idf_for_hash(
    hash: &GeometricHash,
    index: &Option<&FolddiscoIndex>,
    total_structures: f32,
) -> f32 {
    // Use index to get hash counts
    if let Some(ref idx) = index {
        let entries = idx.get_entries(hash.as_u32());
        let hash_count = entries.len();
        if hash_count > 0 {
            return (total_structures / (hash_count as f32)).log2();
        }
    }
    // Default IDF if hash not found in either index
    0.0
}

pub fn parse_threshold_string(threshold_string: Option<String>) -> Vec<f32> {
    if threshold_string.is_none() {
        return Vec::new();
    }
    let threshold_string = threshold_string.unwrap();
    // Remove whitespace
    let threshold_string = threshold_string.replace(" ", "");
    let mut thresholds: Vec<f32> = Vec::new();
    for threshold in threshold_string.split(',') {
        let threshold = threshold.parse::<f32>().expect(
            &log_msg(FAIL, "Failed to parse threshold")
        );
        thresholds.push(threshold);
    }
    thresholds
}

// Doesn't support duplicate hash
// If hash is already in the hash_collection, skip.
fn insert_binned_hash(
    hash_collection: &mut HashMap<GeometricHash, ((usize, usize), bool, f32)>,
    feature: &Vec<f32>, indices: (usize, usize), hash_type: HashType,
    nbin_dist: usize, nbin_angle: usize, multiple_bin: &Option<Vec<(usize, usize)>>,
    is_primary: bool, idf: f32,
) {
    if let Some(multiple_bin) = multiple_bin {
        for (nbin_dist, nbin_angle) in multiple_bin.iter() {
            let hash_value = if *nbin_dist == 0 || *nbin_angle == 0 {
                GeometricHash::perfect_hash_default(feature, hash_type)
            } else {
                GeometricHash::perfect_hash(feature, hash_type, *nbin_dist, *nbin_angle)
            };
            if hash_collection.contains_key(&hash_value) {
                continue;
            } else {
                hash_collection.insert(hash_value, (indices, is_primary, idf));
            }
        }
    } else {
        let hash_value = if nbin_dist == 0 || nbin_angle == 0 {
            GeometricHash::perfect_hash_default(feature, hash_type)
        } else {
            GeometricHash::perfect_hash(feature, hash_type, nbin_dist, nbin_angle)
        };
        if hash_collection.contains_key(&hash_value) {
            return;
        } else {
            hash_collection.insert(hash_value, (indices, is_primary, idf));
        }
    }
}

/// Amino acid pairs to query for one residue pair, besides the observed one.
///
/// The alternatives of both sides are crossed with each other and with the observed
/// residue, so `A164:H,A200:ND` reaches (His, Asp), (His, Asn), (His, observed),
/// (observed, Asp) and (observed, Asn). Unknown one-letter codes map to 255, which
/// is not encodable, and are dropped.
fn substitution_variants(
    i_idx: usize, j_idx: usize, observed: (f32, f32),
    substitution_map: &HashMap<usize, Vec<u8>>,
) -> Vec<(f32, f32)> {
    let sub_i = substitution_map.get(&i_idx);
    let sub_j = substitution_map.get(&j_idx);
    if sub_i.is_none() && sub_j.is_none() {
        return Vec::new();
    }
    let alternatives = |observed_aa: f32, subs: Option<&Vec<u8>>| -> Vec<f32> {
        let mut out = vec![observed_aa];
        if let Some(subs) = subs {
            for aa in subs {
                // 255 marks an unknown one-letter code and is not encodable
                if *aa as usize >= 20 {
                    continue;
                }
                let aa = *aa as f32;
                if !out.contains(&aa) {
                    out.push(aa);
                }
            }
        }
        out
    };
    let alt_i = alternatives(observed.0, sub_i);
    let alt_j = alternatives(observed.1, sub_j);

    let mut variants = Vec::with_capacity(alt_i.len() * alt_j.len());
    for aa_i in &alt_i {
        for aa_j in &alt_j {
            if *aa_i == observed.0 && *aa_j == observed.1 {
                continue; // the observed pair is inserted separately
            }
            variants.push((*aa_i, *aa_j));
        }
    }
    variants
}

pub fn make_query_map(
    path: &String, query_residues: &Vec<(u8, u64)>, hash_type: HashType,
    nbin_dist: usize, nbin_angle: usize, multiple_bin: &Option<Vec<(usize, usize)>>,
    tolerance: &ToleranceConfig,
    amino_acid_substitutions: &Vec<Option<Vec<u8>>>, distance_cutoff: f32, serial_query: bool,
    index: &Option<&FolddiscoIndex>,
    total_structures: f32,
) -> (HashMap<GeometricHash, ((usize, usize), bool, f32)>, Vec<usize>, HashMap<(u8, u8), Vec<(f32, usize)>>) {
    let (compact, _) = read_compact_structure(path).expect("Failed to read compact structure");
    make_query_map_from_structure(
        &compact, query_residues, hash_type, nbin_dist, nbin_angle, multiple_bin,
        tolerance, amino_acid_substitutions, distance_cutoff, serial_query,
        index, total_structures,
    )
}

/// Query map for a torsion-angle ENM conformer ensemble.
///
/// The query structure is wiggled along its low-frequency torsional normal modes and
/// every conformer's hashes are merged into one query, so a single search covers the
/// whole ensemble. This is the "union hashes" strategy: it costs one search rather
/// than one per conformer, and it was the only one of the three sampling strategies
/// that measured better than a plain search per unit of runtime.
///
/// The original structure is hashed first, so its hashes keep `is_primary` and their
/// own IDF. Hashes contributed only by a conformer are marked non-primary and inherit
/// the primary IDF of the same residue pair, which keeps the rare-hash filter in
/// `count_query` comparing like with like.
///
/// Residue indices and the observed distance map come from the original structure —
/// residue matching must be done against the real query, not a wiggled copy.
pub fn make_query_map_with_ensemble(
    path: &String, query_residues: &Vec<(u8, u64)>, hash_type: HashType,
    nbin_dist: usize, nbin_angle: usize, multiple_bin: &Option<Vec<(usize, usize)>>,
    tolerance: &ToleranceConfig,
    amino_acid_substitutions: &Vec<Option<Vec<u8>>>, distance_cutoff: f32, serial_query: bool,
    index: &Option<&FolddiscoIndex>, total_structures: f32,
    num_confs: usize, target_rmsd: f32, nma_modes: usize,
) -> (HashMap<GeometricHash, ((usize, usize), bool, f32)>, Vec<usize>, HashMap<(u8, u8), Vec<(f32, usize)>>) {
    let (compact, _) = read_compact_structure(path).expect("Failed to read compact structure");

    let (mut hash_collection, indices, observed_distance_map) = make_query_map_from_structure(
        &compact, query_residues, hash_type, nbin_dist, nbin_angle, multiple_bin,
        tolerance, amino_acid_substitutions, distance_cutoff, serial_query,
        index, total_structures,
    );
    if num_confs == 0 {
        return (hash_collection, indices, observed_distance_map);
    }

    // Primary IDF per residue pair, for the hashes the conformers add
    let mut primary_idf: HashMap<(usize, usize), f32> = HashMap::default();
    for (edge, is_primary, idf) in hash_collection.values() {
        if *is_primary {
            primary_idf.insert(*edge, *idf);
        }
    }

    let ensemble = match crate::structure::nma::generate_ensemble(
        &compact, num_confs, target_rmsd, nma_modes
    ) {
        Ok(ensemble) => ensemble,
        Err(err) => {
            print_log_msg(WARN, &format!(
                "Torsion-ENM sampling failed ({}); searching the original query only", err
            ));
            return (hash_collection, indices, observed_distance_map);
        }
    };

    // Skip the first entry: generate_ensemble returns the original structure there and
    // it is already hashed above.
    for conformer in ensemble.iter().skip(1) {
        let (conformer_map, _, _) = make_query_map_from_structure(
            conformer, query_residues, hash_type, nbin_dist, nbin_angle, multiple_bin,
            tolerance, amino_acid_substitutions, distance_cutoff, serial_query,
            index, total_structures,
        );
        for (hash, (edge, _, conformer_idf)) in conformer_map {
            if hash_collection.contains_key(&hash) {
                continue;
            }
            let idf = primary_idf.get(&edge).copied().unwrap_or(conformer_idf);
            hash_collection.insert(hash, (edge, false, idf));
        }
    }
    (hash_collection, indices, observed_distance_map)
}

/// Same as `make_query_map` for a structure already in memory, which is what the
/// torsion-ENM conformer ensemble produces.
pub fn make_query_map_from_structure(
    compact: &CompactStructure, query_residues: &Vec<(u8, u64)>, hash_type: HashType,
    nbin_dist: usize, nbin_angle: usize, multiple_bin: &Option<Vec<(usize, usize)>>,
    tolerance: &ToleranceConfig,
    amino_acid_substitutions: &Vec<Option<Vec<u8>>>, distance_cutoff: f32, serial_query: bool,
    index: &Option<&FolddiscoIndex>,
    total_structures: f32,
) -> (HashMap<GeometricHash, ((usize, usize), bool, f32)>, Vec<usize>, HashMap<(u8, u8), Vec<(f32, usize)>>) {
    let mut hash_collection = HashMap::default();
    let mut observed_distance_map: HashMap<(u8, u8), Vec<(f32, usize)>> = HashMap::default();
    
    // Convert residue indices to vector indices
    let mut indices = Vec::new();
    let mut query_residues = query_residues.clone();
    let mut amino_acid_substitutions = amino_acid_substitutions.clone();

    if query_residues.is_empty() {
        // Iterate over all residues and set to query_residues
        for i in 0..compact.num_residues {
            let chain = compact.chain_per_residue[i];
            let residue_index = compact.residue_serial[i];
            query_residues.push((chain, residue_index));
            amino_acid_substitutions.push(None);
        }
    }

    let mut substitution_map: HashMap<usize, Vec<u8>> = HashMap::default();
    
    for (i, (chain, ri)) in query_residues.iter().enumerate() {
        let index = if serial_query { Some(*ri as usize) } else { compact.get_index(&chain, &ri) };
        if let Some(index) = index {
            // convert u8 array to string
            let _residue: String = compact.get_res_name(index).iter().map(|&c| c as char).collect();
            indices.push(index);
            if let Some(substitution) = amino_acid_substitutions[i].clone() {
                substitution_map.insert(index, substitution);
            }
        }
    }
    // Amino acid positions inside the feature vector, if this hash type encodes them
    let aa_indices = hash_type.amino_acid_index().map(|idx| (idx[0], idx[1]));
    let mut expander = FeatureExpander::new(hash_type, nbin_dist, nbin_angle, tolerance);
    // Make combinations
    let comb_iter = CombinationIterator::new(indices.len());
    let mut feature = vec![0.0; 9];
    let mut variant = vec![0.0; 9];
    comb_iter.for_each(|(i, j)| {
        if i == j {
            return;
        }
        let is_feature = get_single_feature(
            indices[i], indices[j], &compact, hash_type, distance_cutoff, &mut feature
        );

        if is_feature {
            // Gather observed distance & aa pairs.
            let aa_dist_info = compact.get_list_amino_acids_and_distances(indices[i], indices[j]);
            if let Some(aa_dist_info) = aa_dist_info {
                // Check if the pair is already in the map
                let aa_pair = (aa_dist_info.0, aa_dist_info.1);
                if observed_distance_map.contains_key(&aa_pair) {
                    observed_distance_map.get_mut(&aa_pair).unwrap().push((aa_dist_info.2, indices[i]));
                } else {
                    observed_distance_map.insert(aa_pair, vec![(aa_dist_info.2, indices[i])]);
                }
            }

            // Insert observed hash - calculate IDF first
            let observed_hash = if nbin_dist == 0 || nbin_angle == 0 {
                GeometricHash::perfect_hash_default(&feature, hash_type)
            } else {
                GeometricHash::perfect_hash(&feature, hash_type, nbin_dist, nbin_angle)
            };
            let idf = calculate_idf_for_hash(&observed_hash, index, total_structures);

            insert_binned_hash(
                &mut hash_collection, &feature, (indices[i], indices[j]),
                hash_type, nbin_dist, nbin_angle, multiple_bin, true, idf
            );

            // Amino acid alternatives requested with `<residue>:<ALT>`. They are
            // applied to the observed geometry *and* to every geometric neighbour
            // below: a substitution and a distance tolerance have to compose,
            // otherwise `164:H` only matches a His sitting in exactly the bins the
            // original residue happened to fall into.
            let aa_variants = match aa_indices {
                Some((a0, a1)) => substitution_variants(
                    indices[i], indices[j], (feature[a0], feature[a1]), &substitution_map
                ),
                None => Vec::new(),
            };
            if let Some((a0, a1)) = aa_indices {
                variant.copy_from_slice(&feature);
                for (aa_i, aa_j) in &aa_variants {
                    variant[a0] = *aa_i;
                    variant[a1] = *aa_j;
                    insert_binned_hash(
                        &mut hash_collection, &variant, (indices[i], indices[j]),
                        hash_type, nbin_dist, nbin_angle, multiple_bin, false, idf
                    );
                }
            }

            // Tolerance neighbourhood, nearest first
            let per_neighbor = 1 + aa_variants.len();
            let mut generated = per_neighbor;
            expander.for_each_neighbor(&feature, |neighbor| {
                variant[..neighbor.len()].copy_from_slice(neighbor);
                insert_binned_hash(
                    &mut hash_collection, &variant, (indices[i], indices[j]),
                    hash_type, nbin_dist, nbin_angle, multiple_bin, false, idf
                );
                if let Some((a0, a1)) = aa_indices {
                    for (aa_i, aa_j) in &aa_variants {
                        variant[a0] = *aa_i;
                        variant[a1] = *aa_j;
                        insert_binned_hash(
                            &mut hash_collection, &variant, (indices[i], indices[j]),
                            hash_type, nbin_dist, nbin_angle, multiple_bin, false, idf
                        );
                    }
                }
                generated += per_neighbor;
                generated + per_neighbor <= MAX_HASHES_PER_PAIR
            });
        }
    });
    (hash_collection, indices, observed_distance_map)
}

pub fn parse_query_string(query_string: &str, mut default_chain: u8) -> (Vec<(u8, u64)>, Vec<Option<Vec<u8>>>) {
    let mut query_residues = Vec::new();
    let mut amino_acid_substitutions = Vec::new();

    if query_string.is_empty() {
        return (query_residues, amino_acid_substitutions);
    }
    if !default_chain.is_ascii_alphabetic() {
        default_chain = b'A';
    }
    // Remove whitespace
    let query_string = query_string.replace(" ", "");
    for segment in query_string.split(',') {
        let (chain, rest) = if let Some(first) = segment.chars().next() {
            // NOTE: 2025-01-15 15:55:19
            // Current querying doesn't support chain ID with more than 1 character
            if first.is_ascii_alphabetic() {
                (first as u8, &segment[1..])
            } else {
                (default_chain, segment)
            }
        } else {
            (default_chain, segment)
        };

        let (range_part, subst_part) = match rest.split_once(':') {
            Some((r, s)) => {
                let sub_vec = s
                    .chars()
                    .filter(|c| is_aa_group_char(*c))
                    .flat_map(|c| map_one_letter_to_u8_vec(c))
                    .collect::<Vec<_>>();
                (r, Some(sub_vec))
            }
            None => (rest, None),
        };

        if range_part.contains('-') {
            let (start_str, end_str) = range_part.split_once('-').expect("Invalid range");
            let start = start_str.parse::<u64>().expect("Invalid start residue");
            let end = end_str.parse::<u64>().expect("Invalid end residue");
            for r in start..=end {
                query_residues.push((chain, r));
                amino_acid_substitutions.push(subst_part.clone());
            }
        } else {
            let residue_num = range_part.parse::<u64>().expect("Invalid residue");
            query_residues.push((chain, residue_num));
            amino_acid_substitutions.push(subst_part);
        }
    }

    (query_residues, amino_acid_substitutions)
}



// ADD TEST
#[cfg(test)]
mod tests {
    use super::*;
    
    fn zinc_finger_query_map(
        tolerance: &ToleranceConfig, substitutions: Vec<Option<Vec<u8>>>,
    ) -> HashMap<GeometricHash, ((usize, usize), bool, f32)> {
        let path = String::from("query/1G2F.pdb");
        let query_residues = vec![(b'F', 207), (b'F', 212), (b'F', 225)];
        let (hash_collection, _indices, _observed_dist_map) = make_query_map(
            &path, &query_residues, HashType::PDBTrRosetta, 16, 4, &None,
            tolerance, &substitutions, 20.0, false, &None, 1000.0
        );
        hash_collection
    }

    #[test]
    fn test_make_query_map() {
        let no_substitution = vec![None; 3];
        let hash_collection = zinc_finger_query_map(
            &ToleranceConfig::default_query(), no_substitution.clone()
        );
        let exact = hash_collection.values().filter(|(_, is_primary, _)| *is_primary).count();
        // At most one observed hash per ordered residue pair; two pairs sharing a
        // hash collapse into a single entry.
        assert!(exact > 0 && exact <= 6, "primary hashes: {}", exact);
        assert!(hash_collection.len() > exact);
    }

    #[test]
    fn wider_tolerance_only_adds_hashes() {
        let no_substitution = vec![None; 3];
        let tight = zinc_finger_query_map(
            &ToleranceConfig::new(vec![0.5], vec![5.0], 0.0, 1), no_substitution.clone()
        );
        let loose = zinc_finger_query_map(
            &ToleranceConfig::new(vec![0.5], vec![5.0], 0.0, 2), no_substitution.clone()
        );
        let elastic = zinc_finger_query_map(
            &ToleranceConfig::new(vec![0.5], vec![5.0], 0.1, 1), no_substitution
        );
        // Both knobs are pure additions: every hash of the tight query survives
        for hash in tight.keys() {
            assert!(loose.contains_key(hash), "radius 2 lost a hash of radius 1");
            assert!(elastic.contains_key(hash), "elastic tolerance lost a hash");
        }
        assert!(loose.len() > tight.len());
        assert!(elastic.len() > tight.len());
    }

    #[test]
    fn substitutions_compose_with_geometric_tolerance() {
        // His at 207 as an alternative. Every hash of the unsubstituted query has to
        // stay, and the His variants have to appear at the tolerance neighbourhood
        // too and not only at the observed geometry.
        let tolerance = ToleranceConfig::new(vec![0.5], vec![5.0], 0.0, 1);
        let plain = zinc_finger_query_map(&tolerance, vec![None; 3]);
        let substituted = zinc_finger_query_map(
            &tolerance, vec![Some(vec![8]), None, None] // 8 = HIS
        );
        for hash in plain.keys() {
            assert!(substituted.contains_key(hash));
        }
        // 2 of the 6 ordered pairs involve residue 207, each gaining as many His
        // variants as it has geometric variants.
        assert!(
            substituted.len() > plain.len() + 2,
            "substitution did not compose with tolerance: {} vs {}",
            substituted.len(), plain.len()
        );
    }

    #[test]
    fn unknown_substitution_code_is_ignored() {
        // 255 is what an unrecognised one-letter code maps to; encoding it would
        // overflow the residue field of the hash.
        let tolerance = ToleranceConfig::default_query();
        let plain = zinc_finger_query_map(&tolerance, vec![None; 3]);
        let bogus = zinc_finger_query_map(&tolerance, vec![Some(vec![255]), None, None]);
        assert_eq!(plain.len(), bogus.len());
    }

    #[test]
    fn test_parse_query_string() {
        let query_string = "A250,B232,C269";
        let query_residues = parse_query_string(query_string, b'A');
        assert_eq!(query_residues, (vec![(b'A', 250), (b'B', 232), (b'C', 269)], vec![None, None, None]));
    }
    #[test]
    fn test_parse_query_string_with_space() {
        let query_string = "A250, A232, A269";
        let query_residues = parse_query_string(query_string, b'A');
        assert_eq!(query_residues, (vec![(b'A', 250), (b'A', 232), (b'A', 269)], vec![None, None, None]));
    }
    
    #[test]
    fn test_parse_query_string_with_space_and_no_chain() {
        let query_string = "250, 232, 269";
        let query_residues = parse_query_string(query_string, b'A');
        assert_eq!(query_residues, (vec![(b'A', 250), (b'A', 232), (b'A', 269)], vec![None, None, None]));
    }

    #[test]
    fn test_parse_query_string_with_aa_substitution() {
        let query_string = "A250:R,B232:K,C269:QK";
        let query_residues = parse_query_string(query_string, b'A');
        // R = 1, K = 11, Q = 5
        assert_eq!(query_residues, (vec![(b'A', 250), (b'B', 232), (b'C', 269)], vec![Some(vec![1]), Some(vec![11]), Some(vec![5, 11])]));
        let query_string = "250:R,232:K,269:QK";
        let query_residues = parse_query_string(query_string, b'A');
        // R = 1, K = 11, Q = 5
        assert_eq!(query_residues, (vec![(b'A', 250), (b'A', 232), (b'A', 269)], vec![Some(vec![1]), Some(vec![11]), Some(vec![5, 11])]));
    }
    #[test]
    fn test_parse_query_string_with_range() {
        let query_string = "A250-252,B232-234,C269:Q";
        let query_residues = parse_query_string(query_string, b'A');
        assert_eq!(query_residues, (vec![
            (b'A', 250), (b'A', 251), (b'A', 252), 
            (b'B', 232), (b'B', 233), (b'B', 234), 
            (b'C', 269),
        ], vec![None, None, None, None, None, None, Some(vec![5])]));
    }
}
