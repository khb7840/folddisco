// File: query.rs
// Created: 2023-12-22 17:00:50
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use rustc_hash::FxHashMap as HashMap;
use crate::geometry::core::{GeometricHash, HashType};
use crate::index::indextable::FolddiscoIndex;
use crate::structure::chain_id::{split_chain_and_rest, ChainId};
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
//
// `is_primary` marks the hash of the observed geometry, as opposed to one the tolerance
// expansion added. It has **no production consumer**: the rare-hash IDF filter that read
// it was removed once it measured worse than no filter at all, and the only remaining
// reader is `test_make_query_map`, which uses it to assert that the expansion produces
// at most one observed hash per ordered residue pair. It is kept rather than deleted
// because the tuple type `((usize, usize), bool, f32)` is spelled in ten places across
// three files with about seventeen edit sites, and the `f32` beside it *is* load-bearing
// (`retrieve::calculate_subgraph_idf` reads it), so removing one bool costs a wide
// refactor for six lines. Written down here so the next maintainer does not redo the grep.
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
    path: &String, query_residues: &Vec<(ChainId, u64)>, hash_type: HashType, 
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

/// Same as `make_query_map` for a structure already in memory.
fn make_query_map_from_structure(
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

/// Widest residue span a single `-q` range may expand to.
///
/// Ranges are expanded eagerly, one residue pushed per position, so the span is an
/// allocation. 100,000 sits above anything a real query can want - the default
/// `--residue` indexing limit is 50,000, i.e. the largest structure folddisco will index
/// unless told otherwise, and the longest human protein (titin) is about 35,000 residues
/// - so it cannot reject a range that any default-built index could match. Its purpose is
/// the other end of the scale: `A204-A2150000000`, one mistyped digit, asks for 2.15e9
/// residues and roughly 80 GB, and used to be killed by the OOM reaper with no message at
/// all. At this cap the two vectors stay a few megabytes and a typo gets a diagnostic.
const MAX_RESIDUE_RANGE_SPAN: u64 = 100_000;

/// Residue number at one end of a range, or a single position.
///
/// A chain letter may be repeated here - `F204-F215` means the same as `F204-215` - but
/// it has to agree with the chain the segment already declared, because a range cannot
/// span two chains.
fn parse_residue_number(token: &str, chain: u8, segment: &str) -> Result<u64, String> {
    let (token_chain, digits) = match token.chars().next() {
        Some(first) if first.is_ascii_alphabetic() => (Some(first as u8), &token[1..]),
        _ => (None, token),
    };
    if let Some(token_chain) = token_chain {
        if token_chain != chain {
            return Err(format!(
                "Query '{}' mixes chain '{}' and chain '{}'; a residue range stays in one chain",
                segment, chain as char, token_chain as char
            ));
        }
    }
    digits.parse::<u64>().map_err(|_| format!(
        "Query '{}' has '{}' where a residue number was expected", segment, token
    ))
}

/// Parse a query string for the CLI, exiting with a diagnostic rather than a panic when
/// it is malformed - a bad `-q` is user input, not a bug.
///
/// This **terminates the process** on a malformed query, which is right for a command
/// line and wrong for anything else. Library consumers should call
/// `parse_query_string_checked` and handle the `Err`.
pub fn parse_query_string(query_string: &str, default_chain: ChainId) -> (Vec<(ChainId, u64)>, Vec<Option<Vec<ChainId>>>) {
    match parse_query_string_checked(query_string, default_chain) {
        Ok((query_residues, amino_acid_substitutions)) => {
            warn_on_duplicate_residues(query_string, &query_residues);
            (query_residues, amino_acid_substitutions)
        }
        Err(err) => {
            print_log_msg(FAIL, &err);
            std::process::exit(1);
        }
    }
}

/// Warn when a query names the same residue more than once, as `B57,B57` or a pair of
/// overlapping ranges like `A1-A5,A3-A7` does.
///
/// Deduplicating would be the real fix, but the length of this list is the denominator of
/// `--covered-node-ratio`, `--max-node-ratio` and the novelty coverage, so dropping
/// repeats changes filtering on the default path - which is the one path this branch
/// keeps byte-identical to upstream. It waits for a benchmark run that can clear it.
fn warn_on_duplicate_residues(query_string: &str, query_residues: &[(ChainId, u64)]) {
    let mut distinct = query_residues.to_vec();
    distinct.sort_unstable();
    distinct.dedup();
    if distinct.len() < query_residues.len() {
        print_log_msg(WARN, &format!(
            "Query '{}' names {} residues but only {} distinct ones. The repeats still \
             count toward the query length, so coverage ratios and the novelty verdict are \
             computed against {} and read low",
            query_string, query_residues.len(), distinct.len(), query_residues.len()
        ));
    }
}

/// Parse a query string into its residues and their allowed substitutions, returning the
/// diagnostic instead of acting on it. This is the entry point for library use.
/// Parse a `-q` motif string into `(chain, residue)` pairs.
///
/// Both spellings of a residue are accepted, so that any `matching_residues`
/// field this tool prints can be pasted straight back in as a query:
///
///   `A250`     legacy, single-letter chain
///   `A_250`    same residue, explicit separator
///   `AA_250`   multi-character chain, separator required
///   `10_250`   numeric chain, separator required
///   `250`      no chain, `default_chain` applies
///
/// Ranges (`A250-252`) and substitutions (`A250:R`) work with either spelling.
pub fn parse_query_string_checked(
    query_string: &str, default_chain: ChainId
)-> Result<(Vec<(ChainId, u64)>, Vec<Option<Vec<ChainId>>>) , String> {
    let mut query_residues = Vec::new();
    let mut amino_acid_substitutions = Vec::new();

    if query_string.is_empty() {
        return Ok((query_residues, amino_acid_substitutions));
    }
    // A blank chain ID is no use as an implicit default; fall back to A as
    // before. A multi-character or numeric one is fine here, because the
    // default is supplied rather than parsed out of the query string.
    let default_chain = if default_chain.is_empty()
        || (default_chain.len() == 1 && !default_chain.first_byte().is_ascii_alphanumeric())
    {
        ChainId::from_byte(b'A')
    } else {
        default_chain
    };
    // Remove whitespace
    let query_string = query_string.replace(" ", "");
    for segment in query_string.split(',') {
        let (chain, rest) = match split_chain_and_rest(segment) {
            (Some(chain), rest) => (chain, rest),
            (None, rest) => (default_chain, rest),
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

        if let Some((start_str, end_str)) = range_part.split_once('-') {
            let start = parse_residue_number(start_str, chain, segment)?;
            let end = parse_residue_number(end_str, chain, segment)?;
            if end < start {
                return Err(format!(
                    "Query '{}' ends before it starts; a range runs from the lower residue",
                    segment
                ));
            }
            if end - start >= MAX_RESIDUE_RANGE_SPAN {
                return Err(format!(
                    "Query '{}' spans {} residues, past the {} a range may expand to. \
                     A range is expanded one residue at a time, so check for a mistyped \
                     digit before this becomes an out-of-memory kill",
                    segment, (end - start).saturating_add(1), MAX_RESIDUE_RANGE_SPAN
                ));
            }
            for r in start..=end {
                query_residues.push((chain, r));
                amino_acid_substitutions.push(subst_part.clone());
            }
        } else {
            let residue_num = parse_residue_number(range_part, chain, segment)?;
            query_residues.push((chain, residue_num));
            amino_acid_substitutions.push(subst_part);
        }
    }

    Ok((query_residues, amino_acid_substitutions))
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
        let path= String::from("query/1G2F.pdb");
        let query_residues = vec![
            (ChainId::from_byte(b'F'), 207), (ChainId::from_byte(b'F'), 212),
            (ChainId::from_byte(b'F'), 225)
        ];
        let amino_acid_substitutions = vec![None; query_residues.len()];
        // let path = String::from("data/serine_peptidases/1aq2.pdb");
        // let query_residues = vec![
        //     (b'A', 250), (b'A', 232), (b'A', 269)
        // ];
        let hash_type = HashType::PDBTrRosetta;
        let (hash_collection, _index_found, _observed_dist_map) = make_query_map(
            &path, &query_residues, hash_type, 16, 4, &None,
            &vec![0.0], &vec![0.0], &amino_acid_substitutions, 20.0, false,
            &None, 1000.0
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
            &ToleranceConfig::new(vec![0.5], vec![5.0], 1), no_substitution.clone()
        );
        let loose = zinc_finger_query_map(
            &ToleranceConfig::new(vec![0.5], vec![5.0], 2), no_substitution
        );
        // A wider radius is a pure addition: every hash of the tight query survives
        for hash in tight.keys() {
            assert!(loose.contains_key(hash), "radius 2 lost a hash of radius 1");
        }
        assert!(loose.len() > tight.len());
    }

    #[test]
    fn substitutions_compose_with_geometric_tolerance() {
        // His at 207 as an alternative. Every hash of the unsubstituted query has to
        // stay, and the His variants have to appear at the tolerance neighbourhood
        // too and not only at the observed geometry.
        let tolerance = ToleranceConfig::new(vec![0.5], vec![5.0], 1);
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

    fn chain(text: &str) -> ChainId {
        ChainId::from_str(text)
    }

    #[test]
    fn test_parse_query_string() {
        let query_string = "A250,B232,C269";
        let query_residues = parse_query_string(query_string, chain("A"));
        assert_eq!(query_residues, (vec![(chain("A"), 250), (chain("B"), 232), (chain("C"), 269)], vec![None, None, None]));
    }
    #[test]
    fn test_parse_query_string_with_space() {
        let query_string = "A250, A232, A269";
        let query_residues = parse_query_string(query_string, chain("A"));
        assert_eq!(query_residues, (vec![(chain("A"), 250), (chain("A"), 232), (chain("A"), 269)], vec![None, None, None]));
    }
    
    #[test]
    fn test_parse_query_string_with_space_and_no_chain() {
        let query_string = "250, 232, 269";
        let query_residues = parse_query_string(query_string, chain("A"));
        assert_eq!(query_residues, (vec![(chain("A"), 250), (chain("A"), 232), (chain("A"), 269)], vec![None, None, None]));
    }

    #[test]
    fn test_parse_query_string_with_aa_substitution() {
        let query_string = "A250:R,B232:K,C269:QK";
        let query_residues = parse_query_string(query_string, chain("A"));
        // R = 1, K = 11, Q = 5
        assert_eq!(query_residues, (vec![(chain("A"), 250), (chain("B"), 232), (chain("C"), 269)], vec![Some(vec![1]), Some(vec![11]), Some(vec![5, 11])]));
        let query_string = "250:R,232:K,269:QK";
        let query_residues = parse_query_string(query_string, chain("A"));
        // R = 1, K = 11, Q = 5
        assert_eq!(query_residues, (vec![(chain("A"), 250), (chain("A"), 232), (chain("A"), 269)], vec![Some(vec![1]), Some(vec![11]), Some(vec![5, 11])]));
    }
    #[test]
    fn range_end_may_repeat_the_chain() {
        // The author's benchmark commands write `F204-F215`; queries_used.tsv writes
        // `F204-215`. Both forms have to mean the same 23 residues, and the first one
        // used to panic on `"F215".parse::<u64>()`.
        let bare = parse_query_string_checked("F204-215,F222-232", b'A').unwrap();
        let prefixed = parse_query_string_checked("F204-F215,F222-F232", b'A').unwrap();
        assert_eq!(bare, prefixed);
        assert_eq!(bare.0.len(), 23);
        assert_eq!(bare.0[0], (b'F', 204));
        assert_eq!(*bare.0.last().unwrap(), (b'F', 232));
    }

    #[test]
    fn an_unbounded_range_is_a_diagnostic_not_an_allocation() {
        // `A204-A2150000000` asks for 2.15e9 residues, about 80 GB, and used to be killed
        // by the OOM reaper with no message. One mistyped digit is exactly the malformed
        // input this parser exists to name.
        let err = parse_query_string_checked("A204-A2150000000", b'A').unwrap_err();
        assert!(err.contains("2150000000") || err.contains("spans"), "{}", err);
        assert!(err.contains("100000"), "the message should name the limit: {}", err);
        // A range that a real index could match is not rejected: the default --residue
        // indexing limit is 50,000, and the cap sits above it
        let (residues, _) = parse_query_string_checked("A1-A50000", b'A').unwrap();
        assert_eq!(residues.len(), 50_000);
        // Exactly at the cap is allowed, one past it is not
        assert_eq!(parse_query_string_checked("A1-A100000", b'A').unwrap().0.len(), 100_000);
        assert!(parse_query_string_checked("A1-A100001", b'A').is_err());
    }

    #[test]
    fn duplicate_residues_are_kept_but_countable() {
        // Overlapping ranges repeat residues 3-5, and the length of this list is the
        // denominator of the coverage ratios, so a perfect 7-residue match would read
        // 7/10. Deduplicating changes filtering on the default path, so the parser keeps
        // them and the caller warns; this test pins the arithmetic the warning is about.
        let (residues, _) = parse_query_string_checked("A1-A5,A3-A7", b'A').unwrap();
        assert_eq!(residues.len(), 10);
        let mut distinct = residues.clone();
        distinct.sort_unstable();
        distinct.dedup();
        assert_eq!(distinct.len(), 7);
        // The plain repeated-residue form behaves the same way
        assert_eq!(parse_query_string_checked("B57,B57,B57", b'A').unwrap().0.len(), 3);
    }

    #[test]
    fn malformed_queries_are_diagnosed_not_panicked() {
        // A range cannot span two chains
        let err = parse_query_string_checked("F204-G215", b'A').unwrap_err();
        assert!(err.contains("chain 'F'") && err.contains("chain 'G'"), "{}", err);
        // Reversed ranges used to yield an empty residue list silently
        let err = parse_query_string_checked("F215-F204", b'A').unwrap_err();
        assert!(err.contains("ends before it starts"), "{}", err);
        // Non-numeric tokens name themselves
        let err = parse_query_string_checked("F20x", b'A').unwrap_err();
        assert!(err.contains("F20x") && err.contains("'20x'"), "{}", err);
        assert!(parse_query_string_checked("F204-", b'A').is_err());
        // The documented good forms keep working
        assert!(parse_query_string_checked("B57,B102,C195", b'A').is_ok());
        assert!(parse_query_string_checked("1-10,11:X", b'A').is_ok());
        assert!(parse_query_string_checked("164:H,195,221,247:ND", b'A').is_ok());
    }

    #[test]
    fn test_parse_query_string_with_range() {
        let query_string = "A250-252,B232-234,C269:Q";
        let query_residues = parse_query_string(query_string, chain("A"));
        assert_eq!(query_residues, (vec![
            (chain("A"), 250), (chain("A"), 251), (chain("A"), 252), 
            (chain("B"), 232), (chain("B"), 233), (chain("B"), 234), 
            (chain("C"), 269),
        ], vec![None, None, None, None, None, None, Some(vec![5])]));
    }

    #[test]
    fn test_parse_query_string_separator_is_optional_for_single_char_chains() {
        // The two spellings of the same motif have to agree, so that a printed
        // `matching_residues` field can be pasted back in as `-q` either way.
        assert_eq!(
            parse_query_string("A250,B232,C269", chain("A")),
            parse_query_string("A_250,B_232,C_269", chain("A"))
        );
        assert_eq!(
            parse_query_string("A250-252,B232-234,C269:Q", chain("A")),
            parse_query_string("A_250-252,B_232-234,C_269:Q", chain("A"))
        );
        assert_eq!(
            parse_query_string("A250:R,B232:K", chain("A")),
            parse_query_string("A_250:R,B_232:K", chain("A"))
        );
    }

    #[test]
    fn test_parse_query_string_multi_char_chain() {
        assert_eq!(
            parse_query_string("AA_250,AB_232,AC_269", chain("A")),
            (vec![(chain("AA"), 250), (chain("AB"), 232), (chain("AC"), 269)], vec![None, None, None])
        );
        // Ranges and substitutions too
        assert_eq!(
            parse_query_string("AA_250-252,AB_232:K", chain("A")),
            (
                vec![(chain("AA"), 250), (chain("AA"), 251), (chain("AA"), 252), (chain("AB"), 232)],
                vec![None, None, None, Some(vec![11])]
            )
        );
    }

    #[test]
    fn test_parse_query_string_numeric_chain() {
        // `10250` is unreadable, so a numeric chain has to be separated. The
        // separated spelling must not be mistaken for residue 10250 of the
        // default chain.
        assert_eq!(
            parse_query_string("10_250,10_252", chain("A")),
            (vec![(chain("10"), 250), (chain("10"), 252)], vec![None, None])
        );
        assert_eq!(
            parse_query_string("1_250", chain("A")),
            (vec![(chain("1"), 250)], vec![None])
        );
    }

    #[test]
    fn test_parse_query_string_default_chain_may_be_multi_char() {
        // A bare residue index belongs to the query structure's own first
        // chain, which can now be multi-character or numeric.
        assert_eq!(
            parse_query_string("250,252", chain("AA")),
            (vec![(chain("AA"), 250), (chain("AA"), 252)], vec![None, None])
        );
        assert_eq!(
            parse_query_string("250", chain("10")),
            (vec![(chain("10"), 250)], vec![None])
        );
        // A blank default still falls back to chain A, as before.
        assert_eq!(
            parse_query_string("250", ChainId::from_byte(b' ')),
            (vec![(chain("A"), 250)], vec![None])
        );
        assert_eq!(
            parse_query_string("250", ChainId::empty()),
            (vec![(chain("A"), 250)], vec![None])
        );
    }
    
    #[test]
    fn test_add_shifted_hashes() {
        let mut hash_collection = HashMap::default();
        let feature = vec![1.0, 2.0, 10.5, 15.2, 0.8, 1.2, 2.1, 0.0, 0.0]; // PDBTrRosetta feature
        let hash_type = HashType::PDBTrRosetta;
        let idf = 5.0; // Example IDF value
        
        // Test that shifted hashes are added for PDBTrRosetta
        _add_shifted_hashes(&feature, &mut hash_collection, 0, 1, hash_type, idf);
        
        // Should have added multiple shifted hash variants
        assert!(hash_collection.len() > 0);
        println!("Added {} shifted hash variants", hash_collection.len());
        
        // Test that no hashes are added for other hash types
        let mut hash_collection_other = HashMap::default();
        _add_shifted_hashes(&feature, &mut hash_collection_other, 0, 1, HashType::TrRosetta, idf);
        assert_eq!(hash_collection_other.len(), 0);
    }
}
