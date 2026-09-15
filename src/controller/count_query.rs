// First-pass scoring: count index hits per target structure before residue matching.

use rayon::prelude::*;

// use std::collections::HashMap;
use rustc_hash::FxHashMap as HashMap;

use crate::index::indextable::FolddiscoIndex;
use crate::index::lookup::LookupTable;
use crate::prelude::GeometricHash;

use super::result::StructureResult;


// Fixed-capacity bit set over structure ids
#[derive(Debug, Clone)]
struct BitVector {
    bits: Vec<u64>,
    capacity: usize,
}

impl BitVector {
    fn new(capacity: usize) -> Self {
        let words_needed = (capacity + 63) / 64;
        Self {
            bits: vec![0u64; words_needed],
            capacity,
        }
    }
    #[inline]
    fn set(&mut self, id: usize) {
        if id < self.capacity {
            let word_idx = id / 64;
            let bit_idx = id % 64;
            self.bits[word_idx] |= 1u64 << bit_idx;
        }
    }
    #[inline]
    fn is_set(&self, id: usize) -> bool {
        if id >= self.capacity {
            return false;
        }
        let word_idx = id / 64;
        let bit_idx = id % 64;
        (self.bits[word_idx] & (1u64 << bit_idx)) != 0
    }
    // #[inline]
    // fn count_ones(&self) -> u32 {
    //     self.bits.iter().map(|&word| word.count_ones()).sum()
    // }
    #[inline]
    fn clear(&mut self) {
        self.bits.fill(0);
    }
}


// Per-target counters, one dense vector per query node
#[derive(Debug, Clone)]
struct CompactEntry {
    node_count: u16,
    edge_count: u32,
    match_count: u32,
    idf_sum: f32,
    // idf_max_per_edge: f32,  // Track max IDF per edge and sum.
    initialized: bool,
}

impl Default for CompactEntry {
    fn default() -> Self {
        Self {
            node_count: 0,
            edge_count: 0,
            match_count: 0,
            idf_sum: 0.0,
            // idf_max_per_edge: 0.0,  // Initialize max IDF per edge
            initialized: false,
        }
    }
}

/// Score every structure hit by `queries`.
///
/// Hashes are grouped by the first residue of their query edge (one "node"); per
/// target, `node_count` is the number of nodes with a hit, `edge_count` the number
/// of residue pairs with a hit, and `idf` the summed log2(N / hash frequency)
/// scaled by `nres^-length_penalty` (default 0.5).
pub fn count_query<'a>(
    queries: &Vec<GeometricHash>, query_map: &HashMap<GeometricHash, ((usize, usize), bool, f32)>,
    index: &FolddiscoIndex,
    lookup: &'a LookupTable, 
    sampling_ratio: Option<f32>, sampling_count: Option<usize>,
    freq_filter: Option<f32>, length_penalty_power: Option<f32>,
) -> Vec<(usize, StructureResult<'a>)> {
    let queries_to_iter = sample_query(queries, index, sampling_ratio, sampling_count);
    let num_ids = lookup.len();
    let lp = length_penalty_power.unwrap_or(0.5);
    
    let node_grouped = build_node_groups(&queries_to_iter, query_map);
    
    // One node group per task; each owns a dense result vector
    let thread_results: Vec<Vec<CompactEntry>> = node_grouped
        .par_iter().map(|(_node, chunk)| {
            let mut local_results: Vec<CompactEntry> = vec![CompactEntry::default(); num_ids];
            let mut node_occupancy = BitVector::new(num_ids);
            let mut edge_occupancy = BitVector::new(num_ids);
            let mut prev_edge = None;

            for (_i, (e, query)) in chunk.iter().enumerate() {
                let need_edge_update = prev_edge.map_or(true, |prev| prev != *e);
                
                if need_edge_update {
                    // Close the previous edge
                    if prev_edge.is_some() {
                        for nid in 0..num_ids {
                            if edge_occupancy.is_set(nid) {
                                local_results[nid].edge_count += 1;
                            }
                        }
                    }
                    edge_occupancy.clear();
                    prev_edge = Some(*e);
                }
                
                let single_queried_values = index.get_entries(query.as_u32());
                let hash_count = single_queried_values.len();
                
                if let Some(freq_filter) = freq_filter {
                    if hash_count as f32 / lookup.len() as f32 > freq_filter {
                        continue;
                    }
                }

                let idf = if hash_count > 0 {
                    (lookup.len() as f32 / hash_count as f32).log2()
                } else {
                    continue;  // Hash absent from the index; nothing to count
                };

                for &value in single_queried_values.iter() {
                    if value >= lookup.len() {
                        continue;
                    }
                
                    let nid = lookup.records()[value].id as usize;
                    let entry = &mut local_results[nid];

                    if !entry.initialized {
                        entry.initialized = true;
                        node_occupancy.set(nid);
                    }
                    entry.match_count += 1;
                    entry.idf_sum += idf;
                    
                    edge_occupancy.set(nid);
                }
            }
            
            // Each node group contributes at most one node per target
            for nid in 0..num_ids {
                if node_occupancy.is_set(nid) {
                    local_results[nid].node_count = 1;
                }
            }
            
            // Close the last edge
            for nid in 0..num_ids {
                if edge_occupancy.is_set(nid) {
                    local_results[nid].edge_count += 1;
                }
            }
            local_results
        })
        .collect();

    // Merge node groups per target
    let results: Vec<(usize, StructureResult<'a>)> = (0..num_ids)
        .into_par_iter()
        .filter_map(|nid| {
            let mut merged_entry = CompactEntry::default();
            let mut found_data = false;
            
            for thread_array in &thread_results {
                let entry = &thread_array[nid];
                if entry.initialized {
                    if !found_data {
                        merged_entry = entry.clone();
                        found_data = true;
                    } else {
                        merged_entry.match_count += entry.match_count;
                        merged_entry.idf_sum += entry.idf_sum;
                        merged_entry.node_count += entry.node_count;
                        merged_entry.edge_count += entry.edge_count;
                    }
                }
            }
            
            if found_data && merged_entry.match_count > 0 {
                let lookup_entry = lookup.entry(nid);
                merged_entry.idf_sum *= (lookup_entry.nres as f32).powf(-lp);

                let sr = StructureResult::new(
                    lookup_entry.name,
                    nid,
                    merged_entry.match_count as usize,
                    merged_entry.node_count as usize,
                    merged_entry.edge_count as usize,
                    merged_entry.idf_sum,
                    lookup_entry.nres,
                    lookup_entry.plddt,
                    lookup_entry.db_key,
                );
                Some((nid, sr))
            } else {
                None
            }
        })
        .collect();

    results
}

/// Keep the rarest hashes: a fraction (`sampling_ratio`) or a count (`sampling_count`).
/// With neither or both set, all hashes are kept.
fn sample_query(
    queries: &Vec<GeometricHash>, 
    index: &FolddiscoIndex,
    sampling_ratio: Option<f32>, 
    sampling_count: Option<usize>,
) -> Vec<GeometricHash> {
    match (sampling_ratio, sampling_count) {
        (None, None) => queries.clone(),
        (Some(sampling_ratio), None) => {
            let mut sampled_queries = queries.par_iter().map(|query| {
                let single_queried_values = index.get_entries(query.as_u32());
                let hash_count = single_queried_values.len();
                (query, hash_count)
            }).collect::<Vec<_>>();
            sampled_queries.sort_by(|a, b| a.1.cmp(&b.1));
            let sample_query_size: usize = sampled_queries.len();
            sampled_queries.truncate((sampling_ratio * sample_query_size as f32).ceil() as usize);
            sampled_queries.into_iter().map(|(query, _)| *query).collect()
        },
        (None, Some(sampling_count)) => {
            let mut sampled_queries = queries.par_iter().map(|query| {
                let single_queried_values = index.get_entries(query.as_u32());
                let hash_count = single_queried_values.len();
                (query, hash_count)
            }).collect::<Vec<_>>();
            sampled_queries.sort_by(|a, b| a.1.cmp(&b.1));
            sampled_queries.truncate(sampling_count);
            sampled_queries.into_iter().map(|(query, _)| *query).collect()
        },
        (Some(_), Some(_)) => queries.clone(),
    }
}

/// Group hashes by the first residue of their edge, sorted by edge.
fn build_node_groups(
    sampled_queries: &[GeometricHash],
    query_map: &HashMap<GeometricHash, ((usize, usize), bool, f32)>,
) -> HashMap<usize, Vec<((usize, usize), GeometricHash)>> {
    let mut node_groups: HashMap<usize, Vec<((usize, usize), GeometricHash)>> = HashMap::default();

    for &query in sampled_queries {
        if let Some(((e0, e1), _, _)) = query_map.get(&query) {
            let edge = (*e0, *e1);
            node_groups.entry(edge.0).or_insert_with(Vec::new).push((edge, query));
        }
    }
    for (_, chunk) in node_groups.iter_mut() {
        chunk.sort_by_key(|(edge, _)| *edge);
    }
    node_groups
}