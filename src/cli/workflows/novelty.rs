// File: novelty.rs
// Description: `folddisco novelty` subcommand workflow.
//
// This workflow answers "is this structural motif novel?" by querying one or
// more reference-database indices and classifying the motif as:
//   NOVEL        – no structure covers the motif above the coverage threshold
//   PARTIAL_MATCH – at least one partial hit found but no full match
//   KNOWN        – at least one structure matches with sufficient coverage
//
// Features implemented
// --------------------
// Feature 1: dedicated `novelty` subcommand with structured TSV report
// Feature 3: sub-motif (pair-level) overlap summary via `--sub-motif`
// Feature 4: multi-database comparison via comma-separated `-i` values
// Feature 6: batch design screening via a query file passed to `-q`

use std::io::{BufRead, Write};

use rayon::prelude::*;

use crate::cli::config::read_index_config_from_file;
use crate::cli::*;
use crate::controller::count_query::count_query;
use crate::controller::filter::StructureFilter;
use crate::controller::io::{check_and_get_indices, get_lookup_and_type, read_compact_structure};
use crate::controller::novelty::{
    count_known_pairs_from_results, NoveltyReport,
};
use crate::controller::query::{make_query_map, parse_query_string, parse_threshold_string};
use crate::controller::result::StructureResult;
use crate::controller::retrieve::retrieval_wrapper;
use crate::index::indextable::load_folddisco_index;
use crate::index::lookup::load_lookup_from_file;
use crate::prelude::*;

#[cfg(feature = "foldcomp")]
use crate::controller::io::get_foldcomp_db_path_with_prefix;
#[cfg(feature = "foldcomp")]
use crate::controller::retrieve::retrieval_wrapper_for_foldcompdb;
#[cfg(feature = "foldcomp")]
use crate::structure::io::fcz::FoldcompDbReader;
#[cfg(feature = "foldcomp")]
use crate::structure::io::StructureFileFormat;

pub const HELP_NOVELTY: &str = "\
usage: folddisco novelty -p <i:PDB> -q <QUERY> -i <i:INDEX[,INDEX2,...]> [OPTIONS]

Classify a structural motif as NOVEL, PARTIAL_MATCH, or KNOWN by searching
one or more reference database indices.

input:
 -p, --pdb <PATH>                 Path of query PDB/mmCIF structure
 -q, --query <STR|FILE>           Comma-separated residue list (e.g. A10,A20,A30)
                                  or a TSV file: pdb_path<TAB>residues[<TAB>output_path]
 -i, --index <PATH[,PATH,...]>    Comma-separated index path(s) to search against [REQUIRED]
 -o, --output <PATH>              Output TSV file [stdout]. Headers are printed when --header is set.

search parameters:
 -t, --threads <INT>              Number of threads [1]
 -d, --distance <FLOAT>           Distance threshold in Å [0.5]
 -a, --angle <FLOAT>              Angle threshold in degrees [5.0]
 --ca-distance <FLOAT>            C-alpha distance threshold in matching [1.0]
 --skip-match                     Skip residue matching (faster, less precise coverage)
 --covered-node-ratio <FLOAT>     Pre-filter: minimum node coverage before matching [0.0]

novelty options:
 --coverage <FLOAT>               Minimum residue coverage fraction to call a motif KNOWN [0.8]
 --sub-motif                      Report sub-motif (residue-pair) overlap summary

general options:
 --header                         Print TSV header line
 -v, --verbose                    Print verbose messages
 -h, --help                       Print this help menu

output columns (TSV):
 design_id | query_residues | novelty_tier | query_residue_count |
 best_hit | best_hit_coverage | best_hit_rmsd | best_hit_idf | best_hit_evalue |
 sub_motif_summary | source_index

examples:
# Check a single designed motif against PDB
folddisco novelty -p design.pdb -q A10,A20,A30 -i index/pdb_folddisco -t 4

# Check against multiple indices (PDB + AFDB50)
folddisco novelty -p design.pdb -q A10,A20,A30 -i index/pdb_folddisco,index/afdb50_folddisco -t 4

# Batch: check a file of designs; write TSV report
folddisco novelty -q designs.txt -i index/afdb50_folddisco -o novelty_report.tsv --header -t 8

# Fast pre-filter only (skip residue matching)
folddisco novelty -p design.pdb -q A10,A20,A30 -i index/pdb_folddisco --skip-match

# Include sub-motif pair overlap summary
folddisco novelty -p design.pdb -q A10,A20,A30 -i index/pdb_folddisco --sub-motif
";

/// Entry point called from `main.rs`.
pub fn check_novelty(env: AppArgs) {
    match env {
        AppArgs::Novelty {
            pdb_path,
            query_string,
            index_paths,
            threads,
            dist_threshold,
            angle_threshold,
            ca_dist_threshold,
            skip_match,
            coverage_threshold,
            sub_motif,
            covered_node_ratio,
            output,
            header,
            verbose,
            help: _,
        } => {
            if verbose { print_logo(); }

            if index_paths.is_empty() {
                eprintln!("{}", HELP_NOVELTY);
                std::process::exit(1);
            }

            // Build thread pool
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();

            // Collect all individual index prefixes (comma-separated)
            let raw_index_list: Vec<String> = index_paths
                .split(',')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect();

            // Expand each comma-separated prefix into validated paths
            let all_indices: Vec<String> = raw_index_list.iter().flat_map(|prefix| {
                check_and_get_indices(Some(prefix.clone()), verbose)
            }).collect();

            if all_indices.is_empty() {
                print_log_msg(FAIL, "No valid index files found");
                std::process::exit(1);
            }

            // Build the list of (pdb_path, query_residues, output_path_per_query) tuples
            let queries: Vec<(String, String, String)> =
                if query_string.ends_with(".txt") || query_string.ends_with(".tsv") {
                    let file = std::fs::File::open(&query_string).expect(
                        &log_msg(FAIL, &format!("Failed to open query file: {}", &query_string))
                    );
                    let reader = std::io::BufReader::new(file);
                    reader.lines().filter_map(|l| l.ok()).map(|line| {
                        let mut parts = line.splitn(3, '\t');
                        let p = parts.next().unwrap_or("").to_string();
                        let q = parts.next().unwrap_or("").to_string();
                        let o = parts.next().unwrap_or("").to_string();
                        (p, q, o)
                    }).collect()
                } else {
                    vec![(pdb_path.clone(), query_string.clone(), output.clone())]
                };

            let dist_thresholds = parse_threshold_string(Some(dist_threshold.clone()));
            let angle_thresholds = parse_threshold_string(Some(angle_threshold.clone()));

            // Open a single global output writer when a global output path is given
            let global_output = if !output.is_empty() { Some(output.clone()) } else { None };

            if let Some(ref path) = global_output {
                if header {
                    let mut f = std::fs::File::create(path)
                        .expect(&log_msg(FAIL, &format!("Failed to create output: {}", path)));
                    writeln!(f, "{}", NoveltyReport::tsv_header())
                        .expect("Failed to write header");
                }
            } else if header {
                println!("{}", NoveltyReport::tsv_header());
            }

            // Process each query
            queries.into_par_iter().for_each(|(pdb_path, query_residues_str, query_output)| {
                let effective_output = if !query_output.is_empty() {
                    query_output.clone()
                } else {
                    global_output.clone().unwrap_or_default()
                };

                // Read the query structure
                let (query_structure, _) = match read_compact_structure(&pdb_path) {
                    Ok(s) => s,
                    Err(_) => {
                        print_log_msg(FAIL, &format!("Failed to read structure: {}", pdb_path));
                        return;
                    }
                };

                let default_chain = if query_structure.chains.is_empty() {
                    b'A'
                } else {
                    query_structure.chains[0]
                };
                let (query_res_parsed, aa_substitutions) =
                    parse_query_string(&query_residues_str, default_chain);

                let residue_count = if query_res_parsed.is_empty() {
                    query_structure.num_residues
                } else {
                    query_res_parsed.len()
                };

                let design_id = pdb_path.clone();

                // Query each index and keep the best report
                let mut best_report: Option<NoveltyReport> = None;

                for index_prefix in &all_indices {
                    let (index, _offset_mmap) =
                        measure_time!(load_folddisco_index(index_prefix), verbose);
                    let (lookup_path, hash_type_path) = get_lookup_and_type(index_prefix);
                    let config = read_index_config_from_file(&hash_type_path);
                    let lookup = measure_time!(load_lookup_from_file(&lookup_path), verbose);

                    let total_structures = lookup.len() as f32;
                    let hash_type = config.hash_type;
                    let num_bin_dist = config.num_bin_dist;
                    let num_bin_angle = config.num_bin_angle;
                    let dist_cutoff = config.grid_width;
                    let multiple_bin = &config.multiple_bin;

                    let (pdb_query_map, query_indices, aa_dist_map) =
                        measure_time!(make_query_map(
                            &pdb_path, &query_res_parsed, hash_type,
                            num_bin_dist, num_bin_angle, multiple_bin,
                            &dist_thresholds, &angle_thresholds,
                            &aa_substitutions, dist_cutoff, false,
                            &Some(&index), total_structures
                        ), verbose);

                    let pdb_query = pdb_query_map.keys().cloned().collect::<Vec<_>>();

                    let structure_filter = StructureFilter::new(
                        0, 0, covered_node_ratio, 0.0, 50000, 0.0,
                        0, 0.0, 0.0, residue_count,
                    );

                    let query_count_map = measure_time!(count_query(
                        &pdb_query, &pdb_query_map, &index, &lookup,
                        None, None, None, None
                    ), verbose);

                    let mut query_count_vec: Vec<(usize, StructureResult)> =
                        query_count_map.into_par_iter()
                            .filter(|(_, v)| structure_filter.filter_before_matching(v))
                            .collect();

                    if verbose {
                        print_log_msg(INFO, &format!(
                            "[{}] {} candidate structure(s)", index_prefix, query_count_vec.len()
                        ));
                    }

                    // Optionally retrieve residue matches
                    if !skip_match {
                        #[cfg(feature = "foldcomp")]
                        let using_foldcomp = config.foldcomp_db.is_some()
                            && config.input_format == StructureFileFormat::FCZDB;

                        #[cfg(feature = "foldcomp")]
                        let foldcomp_db_reader = match config.input_format {
                            StructureFileFormat::FCZDB => {
                                let mut foldcomp_db_path = config.foldcomp_db.clone().unwrap();
                                if !std::path::PathBuf::from(&foldcomp_db_path).is_file() {
                                    if let Some(local) = get_foldcomp_db_path_with_prefix(index_prefix) {
                                        foldcomp_db_path = local;
                                    }
                                }
                                measure_time!(FoldcompDbReader::new(foldcomp_db_path.as_str()), verbose)
                            }
                            _ => FoldcompDbReader::empty(),
                        };

                        query_count_vec.par_iter_mut().for_each(|(_, v)| {
                            #[cfg(not(feature = "foldcomp"))]
                            let retrieval_result = retrieval_wrapper(
                                &v.tid, crate::cli::workflows::query_pdb::MIN_CONNECTED_COMPONENT_SIZE,
                                &pdb_query, hash_type, num_bin_dist, num_bin_angle,
                                multiple_bin, dist_cutoff, &pdb_query_map,
                                &query_structure, &query_indices, &aa_dist_map,
                                ca_dist_threshold, false,
                            );
                            #[cfg(feature = "foldcomp")]
                            let retrieval_result = if using_foldcomp {
                                retrieval_wrapper_for_foldcompdb(
                                    v.db_key,
                                    crate::cli::workflows::query_pdb::MIN_CONNECTED_COMPONENT_SIZE,
                                    &pdb_query, hash_type, num_bin_dist, num_bin_angle,
                                    multiple_bin, dist_cutoff, &pdb_query_map,
                                    &query_structure, &query_indices, &aa_dist_map,
                                    ca_dist_threshold, false, &foldcomp_db_reader,
                                )
                            } else {
                                retrieval_wrapper(
                                    &v.tid,
                                    crate::cli::workflows::query_pdb::MIN_CONNECTED_COMPONENT_SIZE,
                                    &pdb_query, hash_type, num_bin_dist, num_bin_angle,
                                    multiple_bin, dist_cutoff, &pdb_query_map,
                                    &query_structure, &query_indices, &aa_dist_map,
                                    ca_dist_threshold, false,
                                )
                            };
                            v.matching_residues = retrieval_result.0;
                            v.matching_residues_processed = retrieval_result.1;
                            v.max_matching_node_count = retrieval_result.2;
                            v.min_rmsd_with_max_match = retrieval_result.3;
                        });
                    }

                    // Build report for this index
                    let report = if skip_match {
                        NoveltyReport::from_structure_results_skip_match(
                            design_id.clone(),
                            query_residues_str.clone(),
                            residue_count,
                            &query_count_vec,
                            coverage_threshold,
                            total_structures as usize,
                            Some(index_prefix.clone()),
                        )
                    } else {
                        NoveltyReport::from_structure_results(
                            design_id.clone(),
                            query_residues_str.clone(),
                            residue_count,
                            &query_count_vec,
                            coverage_threshold,
                            total_structures as usize,
                            Some(index_prefix.clone()),
                        )
                    };

                    // Attach sub-motif summary if requested
                    let report = if sub_motif {
                        let summary = count_known_pairs_from_results(&query_count_vec, residue_count);
                        report.with_sub_motif_summary(summary)
                    } else {
                        report
                    };

                    // Keep the report that gives the highest coverage (most informative)
                    best_report = Some(match best_report.take() {
                        None => report,
                        Some(prev) => {
                            if report.best_hit_coverage > prev.best_hit_coverage {
                                report
                            } else {
                                prev
                            }
                        }
                    });

                    drop(lookup);
                    drop(index);
                } // each index

                // Emit the result
                let line = match &best_report {
                    Some(r) => r.to_tsv(),
                    None => {
                        NoveltyReport {
                            design_id: design_id.clone(),
                            query_residues: query_residues_str.clone(),
                            novelty_tier: crate::controller::novelty::NoveltyTier::Novel,
                            query_residue_count: residue_count,
                            best_hit: None,
                            best_hit_coverage: 0.0,
                            best_hit_rmsd: f32::MAX,
                            best_hit_idf: 0.0,
                            best_hit_evalue: f64::MAX,
                            sub_motif_summary: None,
                            source_index: None,
                        }.to_tsv()
                    }
                };

                if effective_output.is_empty() {
                    println!("{}", line);
                } else {
                    let mut f = std::fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(&effective_output)
                        .expect(&log_msg(FAIL, &format!(
                            "Failed to open output file: {}", effective_output
                        )));
                    writeln!(f, "{}", line)
                        .expect("Failed to write novelty report line");
                }

                drop(query_structure);
            }); // queries
        }
        _ => {
            eprintln!("{}", HELP_NOVELTY);
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    #[ignore]
    fn test_novelty_workflow() {
        use super::*;
        let env = AppArgs::Novelty {
            pdb_path: "data/serine_peptidases/4cha.pdb".into(),
            query_string: "B57,B102,C195".into(),
            index_paths: "data/serine_peptidases_pdbtr_small".into(),
            threads: 1,
            dist_threshold: "0.5".into(),
            angle_threshold: "5.0".into(),
            ca_dist_threshold: 1.0,
            skip_match: false,
            coverage_threshold: 0.8,
            sub_motif: true,
            covered_node_ratio: 0.0,
            output: "".into(),
            header: true,
            verbose: true,
            help: false,
        };
        check_novelty(env);
    }
}
