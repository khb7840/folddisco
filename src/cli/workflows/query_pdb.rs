// File: query_pdb.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
// Description
// This file contains the workflow for querying PDB files
// When querying PDB files, we need index table and query file.

use std::io::BufRead;
use std::io::Write;

use rayon::prelude::*;

use crate::cli::config::read_index_config_from_file;
use crate::controller::filter::{MatchFilter, StructureFilter};
use crate::controller::mode::QueryMode;
use crate::controller::sort::{MatchSortStrategy, StructureSortStrategy};
use crate::cli::*;
use crate::controller::io::{
    read_compact_structure, 
    check_and_get_indices,
    get_lookup_and_type,
    resolve_tid_path_from_index_prefix,
};
use crate::controller::expand::ToleranceConfig;
use crate::controller::query::{make_query_map, parse_threshold_string};
use crate::controller::count_query::count_query;
use crate::controller::result::{
    convert_structure_query_result_to_match_query_results, 
    sort_and_print_match_query_result, sort_and_print_structure_query_result, StructureResult
};
use crate::controller::retrieve::retrieval_wrapper;
use crate::index::indextable::load_folddisco_index;
use crate::index::lookup::load_lookup_from_file;
use crate::prelude::*;

#[cfg(feature = "foldcomp")]
use crate::controller::retrieve::retrieval_wrapper_for_foldcompdb;
#[cfg(feature = "foldcomp")]
use crate::structure::io::fcz::FoldcompDbReader;
#[cfg(feature = "foldcomp")]
use crate::structure::io::StructureFileFormat;
#[cfg(feature = "foldcomp")]
use crate::controller::io::get_foldcomp_db_path_with_prefix;

pub const HELP_QUERY: &str = "\
usage: folddisco query -p <i:PDB> -q <QUERY> -i <i:INDEX> [OPTIONS] 

input/output:
 -p, --pdb <PATH>                 Path of PDB file to query
 -q, --query <STR>                Query string that specifies residues or a text file containing query
 -i, --index <PATH>               Path of index table to load [REQUIRED]
 -o, --output <PATH>              Output file path [stdout]
 
search parameters:
 -t, --threads <INT>              Number of threads [1]
 -d, --distance <FLOAT>           Distance tolerance in Angstroms. Widest value wins; it is
                                  sub-stepped so every bin it spans is searched [0.5]
 -a, --angle <FLOAT>              Angle tolerance in degrees. Widest value wins [5.0]
 --ca-distance <FLOAT>            C-alpha distance threshold in matching residues [1.0]
 --sampling-count <INT>           Number of sampled hashes to search [all]
 --sampling-ratio <FLOAT>         Sampling ratio for hashes used in searching. For long queries, smaller ratio is recommended [1.0]
 --freq-filter <FLOAT>            Skip queries with hash frequency higher than given ratio [0.0]
 --length-penalty <FLOAT>         Length penalty for searching. Zero means no penalty and higher value gives more penalty to longer structures [0.5]
 --skip-match                     Skip matching residues
 --serial-index                   Handle residue indices serially

non-rigid search:
 --nonrigid                       Preset for deformed motifs: --expand-radius 2. Raises recall
                                  and lowers precision, so pair it with --max-node <n_residues>,
                                  which is what turns the trade into a win. Measured on the
                                  human proteome, F1 over the whole list: 0.9421 -> 0.9641 on a
                                  matched 4-residue zinc motif, 0.9418 -> 0.9577 on a 3-residue
                                  one, 0.8831 -> 0.9160 on a Ser-His-Asp triad against an
                                  independent MEROPS set. Without --max-node the same 4-residue
                                  query is 0.9265 -> 0.9226, i.e. slightly worse. Long segment
                                  queries lose either way (23 residues: 0.9204 -> 0.9117), and
                                  residue matching runs ~11% longer. A looser explicit
                                  --expand-radius is kept

filtering options:
 --total-match <INT>              Filter out structures with less than total match count [0]
 --covered-node <INT>             Filter out structures not covered by given number of nodes with hashes [0]
 --covered-node-ratio <FLOAT>     Filter out structures not covered by given ratio of nodes with hashes [0.0]
 --max-node <INT>                 Filter out structures of maximum matching node size smaller than given value [0]
 --max-node-ratio <FLOAT>         Filter out structures of maximum matching node size smaller than given ratio [0.0]
 --score <FLOAT>                  IDF score cutoff [0.0]
 --connected-node <INT>           Filter out structures/matches with connected node count smaller than given value [0]
 --connected-node-ratio <FLOAT>   Filter out structures/matches with connected node count smaller than given ratio [0.0]
 --num-residue <INT>              Number of residues cutoff [50000]
 --plddt <FLOAT>                  pLDDT cutoff [0.0]
 --rmsd <FLOAT>                   Maximum RMSD cutoff [no limit]
 --tm-score <FLOAT>               Minimum TM-score cutoff [0.0]
 --gdt-ts <FLOAT>                 Minimum GDT-TS cutoff. Thresholds: 1.0Å, 2.0Å, 4.0Å, 8.0Å [0.0]
 --gdt-ha <FLOAT>                 Minimum GDT-HA cutoff. Thresholds: 0.5Å, 1.0Å, 2.0Å, 4.0Å [0.0]
 --chamfer <FLOAT>                Maximum Chamfer distance cutoff. Chamfer distance is mean of nearest neighbor distances between two point clouds [no limit]
 --hausdorff <FLOAT>              Maximum Hausdorff distance cutoff. Hausdorff distance is maximum of nearest neighbor distances between two point clouds [no limit]
 --drmsd <FLOAT>                  Maximum dRMSD cutoff. dRMSD compares the internal distances of the
                                  match instead of superposing it, so a motif bent on a hinge keeps
                                  a low dRMSD where its RMSD is large, so it is the more
                                  forgiving of the two cutoffs on a hinged motif. Sorting by it
                                  did not beat sorting by RMSD on the benchmark [no limit]
 --top <INT>                      Limit output to top N structures based on IDF score [all]

display options:
 --header                         Print header in output
 --web                            Print output for web
 --per-structure                  Print output per structure
 --per-match                      Print output per match. Not working with --skip-match
 --format-output <KEYS>           Comma-separated column names to output
                                  - Per-match: qid, tid, nid, db_key, node_count, idf, rmsd, matching_residues, u_matrix, t_vector,
                                    matching_coordinates, query_residues, tm_score, gdt_ts, gdt_ha, chamfer_distance, hausdorff_distance,
                                    drmsd, max_dist_deviation
                                  - Per-structure: qid, tid, nid, db_key, total_match_count, node_count, edge_count, idf, nres, plddt,
                                    max_node_cov, min_rmsd, min_drmsd, matching_residues, query_residues
                                  - Example: --format-output tid,idf,rmsd,tm_score
 --sort-by <KEYS>                 Comma-separated sort keys with optional :asc or :desc [default: node_count:desc,rmsd:asc]
                                  - Per-match: node_count, idf, rmsd, tm_score, gdt_ts, gdt_ha, chamfer_distance, hausdorff_distance,
                                    drmsd, max_dist_deviation
                                  - Per-structure: max_node_count, node_count, idf, min_rmsd, min_drmsd, total_match_count, edge_count, nres, plddt
                                  - Example: --sort-by tm_score,rmsd or --sort-by idf:desc
 --skip-ca-match                  Print matching residues before C-alpha distance check
 --partial-fit                    Superposition will find the best aligning substructure using LMS (Least Median of Squares)
 --superpose                      Print U, T, CA of matching residues

novelty options:
 --novelty-mode                   Replace the result listing with one verdict line per query:
                                  query_id, NOVEL/PARTIAL_MATCH/KNOWN, best hit, residue coverage,
                                  RMSD (NA with --skip-match), query residues. A query that
                                  hashes to nothing is reported as NO_HASHES, not as NOVEL.
                                  Composes with the sensitivity options above, though on 200
                                  measured comparisons --nonrigid changed exactly one verdict
 --novelty-coverage <FLOAT>       Residue coverage of the best hit needed to call a motif KNOWN.
                                  Coverage is covered/total residues, so for a motif under
                                  5 residues the default demands every residue [0.8]
 --novelty-rmsd <FLOAT>           Best-hit RMSD, in Angstroms, still allowed for KNOWN. A hit
                                  that covers every residue in a different arrangement is not
                                  the same motif: over 100 arbitrary motifs, a third of the
                                  KNOWN verdicts had a best hit worse than 2.0 A, up to 10.4 A.
                                  Covered but over the limit is PARTIAL_MATCH. With --skip-match
                                  no RMSD is computed, so the verdict falls back to coverage
                                  alone and is correspondingly weaker [2.0]

advanced options:
 --expand-radius <INT>            Geometric features of a residue pair allowed in a neighbouring
                                  bin at once. --nonrigid is exactly --expand-radius 2; use this
                                  only to search the observed bins alone (0) or to experiment [1]

general options:
 -v, --verbose                    Print verbose messages
 -h, --help                       Print this help menu

examples:
# Search with default settings. This will print out matching motifs sorted by node count then RMSD.
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6

# Print custom columns (tid, idf, RMSD, and TM-score only)
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --format-output tid,idf,rmsd,tm_score

# Print matches sorted by node count and TM-score
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --sort-by node_count,tm_score

# Print per-structure results sorted by max node count and IDF score
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --per-structure --sort-by max_node_count,idf

# Print per-structure results sorted by IDF only
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --per-structure --sort-by idf

# Query file given as separate text file
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 -d 0.5 -a 5

# Query with amino-acid substitutions and range. 
# Alternative amino acids can be given after colon. Range can be given with dash.
# This will query first 10 residues and 11th residue with subsitution to any amino acid.
folddisco query -p query/4CHA.pdb -q 1-10,11:X -i index/h_sapiens_folddisco -t 6 --serial-index

# Filtering
## Based on connected node and rmsd
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --connected-node 0.75 --rmsd 1.0

## Coverage based filtering & top N filtering without residue matching
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --covered-node 3 --top 1000 --per-structure --skip-match

# Non-rigid search for a deformed motif, ranked by superposition-free deformation
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --nonrigid \\
  --sort-by node_count,drmsd --format-output tid,node_count,idf,rmsd,drmsd,matching_residues

# Maximum recall, accepting a worse ranking: widen every tolerance and rescore by dRMSD
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 \\
  -d 1.0 -a 10 --expand-radius 2 --ca-distance 2.0 --sort-by drmsd

# Is a designed motif novel? One verdict line per design; grep NOVEL to keep the novel ones
folddisco query -p design.pdb -q A10,A20,A30 -i index/pdb_folddisco --novelty-mode --nonrigid

# Fast novelty screen of many designs without residue matching (RMSD prints as NA)
folddisco query -q designs.txt -i index/afdb50_folddisco -t 6 --skip-match \\
  --novelty-mode --novelty-coverage 0.9
";

pub const MIN_CONNECTED_COMPONENT_SIZE: usize = 2;
pub const MAX_NUM_LINES_FOR_WEB: usize = 1000;

/// Expansion radius `--nonrigid` raises the search to.
///
/// Radius 2 lets two features of a residue pair fall on the far side of their bin
/// boundary at the same time. Measured against the human proteome index with this
/// project's own benchmark protocol, F1 over the whole result list versus radius 1,
/// every delta holding its sign in 200 of 200 annotation-dropout replicates:
/// 0.9421 -> 0.9641 on the matched 4-residue zinc command, 0.9418 -> 0.9577 on the
/// 3-residue one, 0.8831 -> 0.9160 on the Ser-His-Asp triad against an independent
/// MEROPS set, against 0.9204 -> 0.9117 on a 23-residue two-segment query.
///
/// It needs the matching step to pay off: on the same 4-residue query without
/// --max-node it is 0.9265 -> 0.9226. Radius 3 measured no better than 2.
///
/// See feature_evaluation.md for the tables these come from.
const NONRIGID_EXPAND_RADIUS: usize = 2;

pub fn query_pdb(env: AppArgs) {
    match env {
        AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            skip_match,
            dist_threshold,
            angle_threshold,
            ca_dist_threshold,
            expand_radius,
            nonrigid,
            total_match_count,
            covered_node_count,
            covered_node_ratio,
            max_matching_node_count,
            max_matching_node_ratio,
            idf_score_cutoff,
            connected_node_count,
            connected_node_ratio,
            num_res_cutoff,
            plddt_cutoff,
            rmsd_cutoff,
            tm_score_cutoff,
            gdt_ts_cutoff,
            gdt_ha_cutoff,
            chamfer_distance_cutoff,
            hausdorff_distance_cutoff,
            drmsd_cutoff,
            top_n,
            web_mode,
            sampling_count,
            sampling_ratio,
            freq_filter,
            length_penalty,
            sort_by,
            format_output,
            output_per_structure,
            output_per_match,
            output_with_superpose,
            skip_ca_match,
            partial_fit,
            header,
            serial_query,
            output,
            novelty_mode,
            novelty_coverage_threshold,
            novelty_rmsd_threshold,
            verbose,
            help: _,
        } => {
            if verbose { print_logo(); }
            // help is already handled in main.rs
            // Check if arguments are valid
            if index_path.is_none() {
                eprintln!("{}", HELP_QUERY);
                std::process::exit(1);
            }
            
            // Determine query mode first to decide which sorting strategy to use
            let query_mode = QueryMode::from_flags(
                skip_match, web_mode, output_per_structure, output_per_match
            );
            
            // Error handling
            if query_mode == QueryMode::ContradictoryPrintError {
                print_log_msg(FAIL, 
                    "Cannot print output per structure and per match at the same time. Use either --per-structure or --per-match"
                );
                std::process::exit(1);
            }
            
            
            // Parse the appropriate sorting strategy based on query mode
            // For per-structure output modes, use StructureSortStrategy
            // For per-match output modes, use MatchSortStrategy
            let use_structure_sort = matches!(
                query_mode,
                QueryMode::PerStructure | QueryMode::SkipMatch
            );
            
            let match_sort_strategy = if !use_structure_sort {
                // If sort_by is given, parse it
                if sort_by.is_empty() {
                    MatchSortStrategy::default()
                } else {
                    MatchSortStrategy::from_str(&sort_by)
                        .unwrap_or_else(|e| {
                            print_log_msg(FAIL, &format!("Error parsing --sort-by: {}", e));
                            std::process::exit(1);
                        })
                }
            } else {
                MatchSortStrategy::default()
            };
            
            let structure_sort_strategy = if use_structure_sort {
                // If sort_by is given, parse it
                if sort_by.is_empty() {
                    StructureSortStrategy::default()
                } else {
                    StructureSortStrategy::from_str(&sort_by)
                        .unwrap_or_else(|e| {
                            print_log_msg(FAIL, &format!("Error parsing --sort-by: {}", e));
                            std::process::exit(1);
                        })
                }
            } else {
                StructureSortStrategy::default()
            };
            
            // Parse format_output if provided
            let parsed_columns: Option<Vec<String>> = format_output.map(|cols| {
                cols.split(',').map(|s| s.trim().to_string()).collect()
            });
            let column_refs: Option<Vec<&str>> = parsed_columns.as_ref().map(|cols| {
                cols.iter().map(|s| s.as_str()).collect()
            });

            // Print query mode and sorting strategy
            if verbose  {
                if use_structure_sort {
                    print_log_msg(INFO, &format!("Printing results {} sorting with {}", query_mode, structure_sort_strategy));
                } else {
                    print_log_msg(INFO, &format!("Printing results {} sorting with {}", query_mode, match_sort_strategy));
                }
            }
            
            // Print query information
            if verbose {
                // If pdb_path is empty
                if pdb_path.is_empty() {
                    print_log_msg(INFO, &format!("Querying {} to {}", &query_string, &index_path.clone().unwrap()));
                } else {
                    print_log_msg(INFO, &format!("Querying {}:{} to {}", &pdb_path, &query_string, &index_path.clone().unwrap()));
                }
                // NOTE: If needed, print filter information
            }
            
            // Set thread pool if there's no global thread pool yet
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
            
            // Get index paths
            let index_paths = check_and_get_indices(index_path.clone(), verbose);
            if verbose {
                print_log_msg(INFO, &format!("Found {} index file(s). Querying with {} threads", index_paths.len(), threads));
            }
            
            // Always use big index
            let index_prefix = index_paths[0].clone();
            let (index, offset_mmap) = measure_time!(load_folddisco_index(&index_prefix), verbose);
            
            // Load lookup and config
            let (lookup_path, hash_type_path) = get_lookup_and_type(&index_prefix);
            let config = read_index_config_from_file(&hash_type_path);
            let lookup = measure_time!(load_lookup_from_file(&lookup_path), verbose);

            
            let queries = if query_string.ends_with(".txt") || query_string.ends_with(".tsv") {
                // Read file and get path, query, output by line
                let mut queries: Vec<(String, String, String)> = Vec::new();
                let file = std::fs::File::open(&query_string).expect(
                    &log_msg(FAIL, &format!("Failed to open query file: {}", &query_string))
                );
                let reader = std::io::BufReader::new(file);
                for line in reader.lines() {
                    let line = line.expect("Failed to read line");
                    let mut split = line.split('\t');
                    let pdb_path = split.next().expect("Failed to get pdb path").to_string();
                    let query_string = split.next().unwrap_or("").to_string();
                    let output_path = split.next().unwrap_or("").to_string();
                    queries.push((pdb_path, query_string, output_path));
                }
                queries
            } else {
                vec![(pdb_path.clone(), query_string.clone(), output.clone())]
            };

            // Every other output mode truncates its file (`File::create` in result.rs).
            // Novelty verdicts are appended, so that a batch of queries sharing one
            // output file accumulates instead of racing to overwrite each other; that
            // makes a rerun of the same command double the file unless each distinct
            // path is emptied once, here, before the queries start.
            if novelty_mode {
                let mut truncated: Vec<&str> = Vec::new();
                for (_, _, output_path) in queries.iter() {
                    if output_path.is_empty() || truncated.contains(&output_path.as_str()) {
                        continue;
                    }
                    std::fs::File::create(output_path).expect(
                        &log_msg(FAIL, &format!("Failed to create file: {}", output_path))
                    );
                    truncated.push(output_path);
                }
            }

            let dist_thresholds = parse_threshold_string(Some(dist_threshold.clone()));
            let angle_thresholds = parse_threshold_string(Some(angle_threshold.clone()));

            // `--nonrigid` only raises the expansion radius; it never lowers an
            // explicit --expand-radius.
            let expand_radius = if nonrigid { expand_radius.max(NONRIGID_EXPAND_RADIUS) } else { expand_radius };
            let tolerance = ToleranceConfig::new(
                dist_thresholds, angle_thresholds, expand_radius
            );
            if verbose {
                print_log_msg(INFO, &format!(
                    "Tolerance: distance {} A, angle {} deg, expansion radius {}",
                    &dist_threshold, &angle_threshold, expand_radius
                ));
            }

            // Load foldcomp db 
            #[cfg(feature = "foldcomp")]
            let using_foldcomp = config.foldcomp_db.is_some() && config.input_format == StructureFileFormat::FCZDB;

            #[cfg(feature = "foldcomp")]
            let foldcomp_db_reader = match config.input_format {
                StructureFileFormat::FCZDB => {
                    if !skip_match {
                        let mut foldcomp_db_path = config.foldcomp_db.clone().unwrap();
                        // If foldcomp_db_path is not a valid path, check foldcomp db with index prefix
                        if !std::path::PathBuf::from(&foldcomp_db_path).is_file() {
                            let local_foldcomp_db_path = get_foldcomp_db_path_with_prefix(&index_prefix);
                            if local_foldcomp_db_path.is_some() {
                                foldcomp_db_path = local_foldcomp_db_path.unwrap();
                            }
                        }
                        measure_time!(FoldcompDbReader::new(foldcomp_db_path.as_str()), verbose)
                    } else {
                        FoldcompDbReader::empty()
                    }
                },
                _ => FoldcompDbReader::empty(),
            };

            // #[cfg(not(feature = "foldcomp"))]
            // let using_foldcomp = false;

            // Iterate over queries
            queries.into_par_iter().for_each(|(pdb_path, query_string, output_path)| {
                let (query_structure, _) = read_compact_structure(&pdb_path).expect(
                    &log_msg(FAIL, &format!("Failed to read structure: {}", &pdb_path))
                );
                
                let (query_residues, aa_substitutions) = parse_query_string(&query_string, query_structure.chains[0]);
                
                let residue_count = if query_residues.is_empty() {
                    query_structure.num_residues
                } else {
                    query_residues.len()
                };
                let query_string = if query_residues.is_empty() {
                    query_string
                } else {
                    let query_residues = query_residues.clone();
                    // query_residues.sort();
                    res_chain_to_string(&query_residues)
                };

                // Get query map for the index
                let hash_type = config.hash_type;
                let num_bin_dist = config.num_bin_dist;
                let num_bin_angle = config.num_bin_angle;
                let dist_cutoff = config.grid_width;
                let multiple_bin = &config.multiple_bin;
                let total_structures = lookup.len() as f32;
                        
                let (pdb_query_map, query_indices, aa_dist_map ) = measure_time!(make_query_map(
                    &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, multiple_bin,
                    &tolerance, &aa_substitutions, dist_cutoff, serial_query,
                    &Some(&index), total_structures
                ), verbose);

                let pdb_query = pdb_query_map.keys().cloned().collect::<Vec<_>>();
                if verbose {
                    print_log_msg(INFO, &format!("Expanded query into {} hashes", pdb_query.len()));
                }
                // A query whose residues are all further apart than the index cutoff
                // produces no hashes at all, so it matches nothing. That is an
                // unrepresentable query, not a result: saying so is the difference
                // between "this motif is new" and "this motif cannot be searched".
                if pdb_query.is_empty() {
                    print_log_msg(WARN, &format!(
                        "{}:{} produced no hashes; its residues are further apart than the index \
                         distance cutoff, so nothing can match it",
                        &pdb_path, &query_string
                    ));
                    if novelty_mode {
                        print_novelty_verdict(
                            &novelty_no_hash_line(&pdb_path, &query_string), &output_path
                        );
                        return;
                    }
                }
                // Make filters out of filtering parameters
                let structure_filter = StructureFilter::new(
                    total_match_count, covered_node_count, covered_node_ratio,
                    idf_score_cutoff, num_res_cutoff, plddt_cutoff,
                    max_matching_node_count, max_matching_node_ratio, rmsd_cutoff,
                    drmsd_cutoff, residue_count,
                );

                let query_count_map = measure_time!(count_query(
                    &pdb_query, &pdb_query_map, &index, &lookup,
                    sampling_ratio, sampling_count, freq_filter, length_penalty
                ), verbose);
                let mut query_count_vec: Vec<(usize, StructureResult)> = query_count_map.into_par_iter().filter(|(_k, v)| {
                    structure_filter.filter_before_matching(v)
                }).collect();

                if verbose {
                    print_log_msg(INFO, &format!("Found {} structures from inverted index", query_count_vec.len()));
                }

                        // 
                measure_time!(query_count_vec.par_sort_by(|a, b| b.1.idf.partial_cmp(&a.1.idf).unwrap()), verbose);
                // Apply top N filter if top_n is not usize::MAX
                if top_n != usize::MAX {
                    if verbose {
                        print_log_msg(INFO, &format!("Limiting result to top {} structures", top_n));
                    }
                    query_count_vec.truncate(top_n);
                }
                        
                // IF retrieve is true, retrieve matching residues
                if !skip_match {
                    measure_time!(query_count_vec.par_iter_mut().for_each(|(_, v)| {
                        let joined_path;
                        let resolved_tid: &str = if std::path::Path::new(&v.tid).is_file() {
                            &v.tid
                        } else {
                            joined_path = resolve_tid_path_from_index_prefix(&v.tid, &index_prefix);
                            &joined_path
                        };
                        
                        #[cfg(not(feature = "foldcomp"))]
                        let retrieval_result = retrieval_wrapper(
                            resolved_tid, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                            hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                            &pdb_query_map, &query_structure, &query_indices,
                            &aa_dist_map, ca_dist_threshold, partial_fit
                        );
                        #[cfg(feature = "foldcomp")]
                        let retrieval_result = if using_foldcomp {
                            retrieval_wrapper_for_foldcompdb(
                                v.db_key, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                &pdb_query_map, &query_structure, &query_indices,
                                &aa_dist_map, ca_dist_threshold, partial_fit,
                                &foldcomp_db_reader
                            )
                        } else {
                            retrieval_wrapper(
                                resolved_tid, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                &pdb_query_map, &query_structure, &query_indices,
                                &aa_dist_map, ca_dist_threshold, partial_fit,
                            )
                        };
                        v.matching_residues = retrieval_result.0;
                        v.matching_residues_processed = retrieval_result.1;
                        v.max_matching_node_count = retrieval_result.2;
                        v.min_rmsd_with_max_match = retrieval_result.3;
                        v.min_drmsd_with_max_match = retrieval_result.4;
                    }), verbose);

                    // Filter query_count_vec with reasonable retrieval results
                    query_count_vec.retain(|(_, v)| structure_filter.filter_after_matching(v));
                }
                let mut queried_from_indices = query_count_vec;
                drop(query_residues);
                let evalue_cutoff = f64::MAX; // Currently not used
                let match_filter= MatchFilter::new(
                    connected_node_count, connected_node_ratio, idf_score_cutoff, evalue_cutoff,
                    rmsd_cutoff, tm_score_cutoff, gdt_ts_cutoff, gdt_ha_cutoff,
                    chamfer_distance_cutoff, hausdorff_distance_cutoff, drmsd_cutoff,
                    residue_count,
                );

                match query_mode {
                    QueryMode::PerMatch => {
                        let mut match_results = convert_structure_query_result_to_match_query_results(
                            &queried_from_indices, skip_ca_match, 
                            total_structures as usize, residue_count
                        );
                        match_results.retain(|(_, v)| match_filter.filter(v));
                        if novelty_mode {
                            // Highest residue coverage wins, lower RMSD breaks ties
                            let best = match_results.iter().map(|(_, v)| v).max_by(
                                |a, b| a.node_count.cmp(&b.node_count).then_with(
                                    || b.rmsd.partial_cmp(&a.rmsd).unwrap_or(std::cmp::Ordering::Equal)
                                )
                            ).map(|v| (v.tid, v.node_count, Some(v.rmsd)));
                            print_novelty_verdict(&novelty_verdict_line(
                                &pdb_path, &query_string, residue_count,
                                novelty_coverage_threshold, novelty_rmsd_threshold, best,
                            ), &output_path);
                        } else {
                            sort_and_print_match_query_result(
                                &mut match_results, top_n, 
                                &output_path, &pdb_path, &query_string, 
                                column_refs.as_deref(), output_with_superpose, header, verbose,
                                match_sort_strategy.clone(),
                            );
                        }
                    }
                    QueryMode::Web => {
                        let mut match_results = convert_structure_query_result_to_match_query_results(
                            &queried_from_indices, skip_ca_match, 
                            total_structures as usize, residue_count
                        );
                        match_results.retain(|(_, v)| match_filter.filter(v));
                        // If web, set superpose to true.
                        sort_and_print_match_query_result(
                            &mut match_results, MAX_NUM_LINES_FOR_WEB,
                            &output_path, &pdb_path, &query_string, 
                            column_refs.as_deref(), true, header, verbose,
                            match_sort_strategy.clone(),
                        );
                    }
                    QueryMode::PerStructure | QueryMode::SkipMatch => {
                        if novelty_mode {
                            // Without matching, the hashes covering a node are all we know,
                            // and the IDF score has to break ties instead of the RMSD
                            let best = queried_from_indices.iter().map(|(_, v)| v).max_by(|a, b| if skip_match {
                                a.node_count.cmp(&b.node_count).then_with(
                                    || a.idf.partial_cmp(&b.idf).unwrap_or(std::cmp::Ordering::Equal)
                                )
                            } else {
                                a.max_matching_node_count.cmp(&b.max_matching_node_count).then_with(
                                    || b.min_rmsd_with_max_match.partial_cmp(&a.min_rmsd_with_max_match)
                                        .unwrap_or(std::cmp::Ordering::Equal)
                                )
                            }).map(|v| if skip_match {
                                (v.tid, v.node_count, None)
                            } else {
                                (v.tid, v.max_matching_node_count, Some(v.min_rmsd_with_max_match))
                            });
                            print_novelty_verdict(&novelty_verdict_line(
                                &pdb_path, &query_string, residue_count,
                                novelty_coverage_threshold, novelty_rmsd_threshold, best,
                            ), &output_path);
                        } else {
                            sort_and_print_structure_query_result(
                                &mut queried_from_indices, &output_path, 
                                &pdb_path, &query_string, column_refs.as_deref(), header, verbose, structure_sort_strategy.clone()
                            );
                        }
                    }
                    QueryMode::ContradictoryPrintError => {
                        // This should have been caught earlier, but handle it just in case
                        print_log_msg(FAIL, "Invalid query mode");
                        std::process::exit(1);
                    }
                }
                drop(queried_from_indices);
                drop(query_structure);
            }); // queries
            drop(lookup);
            drop(offset_mmap);
            drop(index);
        }, // AppArgs::Query
        _ => {
            eprintln!("{}", HELP_QUERY);
            std::process::exit(1);
        }
    }
}

/// Verdict for a query that produced no hashes, and so cannot match anything. Kept
/// distinct from NOVEL: an unrepresentable query is not a discovery.
fn novelty_no_hash_line(query_id: &str, query_residues: &str) -> String {
    format!("{}\tNO_HASHES\tNA\t0.0000\tNA\t{}", query_id, query_residues)
}

/// One-line novelty verdict for a query, tab separated:
/// `query_id, NOVEL|PARTIAL_MATCH|KNOWN, best hit or NA, coverage, RMSD or NA, query residues`
///
/// KNOWN needs both enough coverage and a close enough best hit; anything covered but
/// failing either is PARTIAL_MATCH. (A query that produced no hashes at all is reported
/// as NO_HASHES by the caller, and never reaches this function.)
///
/// `best` is the highest-coverage hit as (target id, covered residues, RMSD). Its
/// RMSD is `None` when residue matching was skipped, printed as `NA` instead of the
/// 0.0 it was never computed into. No hit at all is the NOVEL verdict.
fn novelty_verdict_line(
    query_id: &str, query_residues: &str, query_residue_count: usize,
    coverage_threshold: f32, rmsd_threshold: f32, best: Option<(&str, usize, Option<f32>)>,
) -> String {
    let (tid, covered_nodes, rmsd) = match best {
        Some(best) => best,
        None => return format!("{}\tNOVEL\tNA\t0.0000\tNA\t{}", query_id, query_residues),
    };
    // A query without residues has nothing to cover, so it is never KNOWN
    let coverage = if query_residue_count > 0 {
        covered_nodes as f32 / query_residue_count as f32
    } else {
        0.0
    };
    // Coverage alone does not make a motif known: the same residues in a different
    // arrangement cover everything and are still a different motif, so the RMSD the
    // line already prints has to gate the verdict too. With --skip-match there is no
    // RMSD to gate on and coverage is all there is.
    let close_enough = rmsd.map_or(true, |rmsd| rmsd <= rmsd_threshold);
    let tier = if query_residue_count > 0 && coverage >= coverage_threshold && close_enough {
        "KNOWN"
    } else if covered_nodes > 0 {
        "PARTIAL_MATCH"
    } else {
        "NOVEL"
    };
    let rmsd = rmsd.map_or("NA".to_string(), |rmsd| format!("{:.4}", rmsd));
    format!("{}\t{}\t{}\t{:.4}\t{}\t{}", query_id, tier, tid, coverage, rmsd, query_residues)
}

/// Append one verdict line to the output file, or print it to stdout. Verdicts are
/// appended because a batch of queries usually collects into a single output file; the
/// caller empties each distinct path once before the queries start, so a rerun replaces
/// the file rather than growing it. One `write_all` per line keeps lines from
/// interleaving on a local filesystem (NFS makes no such promise).
fn print_novelty_verdict(line: &str, output_path: &str) {
    if output_path.is_empty() {
        println!("{}", line);
        return;
    }
    let mut file = std::fs::OpenOptions::new().create(true).append(true).open(output_path).expect(
        &log_msg(FAIL, &format!("Failed to open file: {}", output_path))
    );
    file.write_all(format!("{}\n", line).as_bytes()).expect(
        &log_msg(FAIL, &format!("Failed to write to file: {}", output_path))
    );
}

pub fn res_chain_to_string(res_chain: &Vec<(u8, u64)>) -> String {
    let mut output = String::new();
    for (i, (chain, res)) in res_chain.iter().enumerate() {
        output.push_str(&format!("{}{}", *chain as char, res));
        if i < res_chain.len() - 1 {
            output.push(',');
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    #[ignore]
    fn test_query_pdb_workflow() {
        let pdb_path = String::from("data/serine_peptidase/4cha.pdb");
        let query_string = String::from("B57,B102,C195");
        let threads = 1;
        let index_path = Some(String::from("data/serine_peptidases_pdbtr_small"));
        let env = AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            skip_match: false,
            dist_threshold: String::from("0.5"),
            angle_threshold: String::from("5.0"),
            ca_dist_threshold: 1.0,
            expand_radius: 1,
            nonrigid: false,
            total_match_count: 0,
            covered_node_count: 0,
            covered_node_ratio: 0.0,
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            idf_score_cutoff: 0.0,
            connected_node_count: 0,
            connected_node_ratio: 0.0,
            num_res_cutoff: 3000,
            plddt_cutoff: 0.0,
            rmsd_cutoff: 1.0,
            tm_score_cutoff: 0.0,
            gdt_ts_cutoff: 0.0,
            gdt_ha_cutoff: 0.0,
            chamfer_distance_cutoff: 0.0,
            hausdorff_distance_cutoff: 0.0,
            drmsd_cutoff: 0.0,
            top_n: 1000,
            web_mode: false,
            sampling_count: None,
            sampling_ratio: None,
            freq_filter: None,
            length_penalty: None,
            sort_by: String::from("node_count,rmsd"),
            format_output: None,
            output_per_structure: false,
            output_per_match: true,
            output_with_superpose: false,
            skip_ca_match: false,
            partial_fit: false,
            header: true,
            serial_query: false,
            output: String::from(""),
            novelty_mode: false,
            novelty_coverage_threshold: 0.8,
            novelty_rmsd_threshold: 2.0,
            verbose: true,
            help: false,
        };
        query_pdb(env);
    }
    #[test]
    #[ignore]
    fn test_query_with_foldcompdb() {
        #[cfg(feature = "foldcomp")] {
            let pdb_path = String::from("data/foldcomp/example_db:d1asha_");
            let query_string = String::from("1,2,3,4");
            let threads = 1;
            let index_path = Some(String::from("data/example_db_folddisco_db"));
            let env = AppArgs::Query {
                pdb_path,
                query_string,
                threads,
                index_path,
                skip_match: false,
                dist_threshold: String::from("0.5"),
                angle_threshold: String::from("5.0"),
                ca_dist_threshold: 1.0,
                expand_radius: 1,
                nonrigid: false,
                total_match_count: 0,
                covered_node_count: 0,
                covered_node_ratio: 0.0,
                idf_score_cutoff: 0.0,
                connected_node_count: 0,
                connected_node_ratio: 0.0,
                max_matching_node_count: 0,
                max_matching_node_ratio: 0.0,
                num_res_cutoff: 3000,
                plddt_cutoff: 0.0,
                rmsd_cutoff: 1.0,
                tm_score_cutoff: 0.0,
                gdt_ts_cutoff: 0.0,
                gdt_ha_cutoff: 0.0,
                chamfer_distance_cutoff: 0.0,
                hausdorff_distance_cutoff: 0.0,
                drmsd_cutoff: 0.0,
                top_n: 1000,
                web_mode: false,
                sampling_count: None,
                sampling_ratio: None,
                freq_filter: None,
                length_penalty: None,
                sort_by: String::from("node_count,rmsd"),
                format_output: None,
                output_per_structure: true,
                output_per_match: false,
                output_with_superpose: true,
                skip_ca_match: false,
                header: true,
                serial_query: false,
                output: String::from(""),
                novelty_mode: false,
                novelty_coverage_threshold: 0.8,
                novelty_rmsd_threshold: 2.0,
                verbose: true,
                partial_fit: false,
                help: false,
            };
            query_pdb(env);
        }
    }
    #[test]
    #[ignore]
    fn test_query_pdb_with_file() {
        let pdb_path = String::from("");
        let query_string = String::from("data/query.tsv");
        let threads = 4;
        let index_path = Some(String::from("analysis/e_coli/test"));
        let env = AppArgs::Query {
            pdb_path,
            query_string,
            threads, 
            index_path,
            skip_match: true,
            dist_threshold: String::from("0.5"),
            angle_threshold: String::from("5.0"),
            ca_dist_threshold: 1.0,
            expand_radius: 1,
            nonrigid: false,
            total_match_count: 0,
            covered_node_count: 0,
            covered_node_ratio: 0.0,
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            idf_score_cutoff: 0.0,
            connected_node_count: 0,
            connected_node_ratio: 0.0,
            num_res_cutoff: 3000,
            plddt_cutoff: 0.0,
            rmsd_cutoff: 1.0,
            tm_score_cutoff: 0.0,
            gdt_ts_cutoff: 0.0,
            gdt_ha_cutoff: 0.0,
            chamfer_distance_cutoff: 0.0,
            hausdorff_distance_cutoff: 0.0,
            drmsd_cutoff: 0.0,
            top_n: 1000,
            web_mode: false,
            sampling_count: None,
            sampling_ratio: None,
            freq_filter: None,
            length_penalty: None,
            sort_by: String::from("node_count,rmsd"),
            format_output: None,
            output_per_structure: true,
            output_per_match: false,
            output_with_superpose: true,
            skip_ca_match: false,
            partial_fit: false,
            header: true,
            serial_query: false,
            output: String::from(""),
            novelty_mode: false,
            novelty_coverage_threshold: 0.8,
            novelty_rmsd_threshold: 2.0,
            verbose: true,
            help: false,
        };
        query_pdb(env);
    }

    #[test]
    fn test_novelty_verdict_line() {
        let residues = "A10,A20,A30";
        // No hit at all is the NOVEL verdict, and it is still emitted
        assert_eq!(
            novelty_verdict_line("design.pdb", residues, 3, 0.8, 2.0, None),
            "design.pdb\tNOVEL\tNA\t0.0000\tNA\tA10,A20,A30"
        );
        // Fully covered motif with a matched RMSD
        assert_eq!(
            novelty_verdict_line("design.pdb", residues, 3, 0.8, 2.0, Some(("1abc", 3, Some(0.5)))),
            "design.pdb\tKNOWN\t1abc\t1.0000\t0.5000\tA10,A20,A30"
        );
        // Covered below the threshold
        assert_eq!(
            novelty_verdict_line("design.pdb", residues, 3, 0.8, 2.0, Some(("1abc", 2, Some(1.25)))),
            "design.pdb\tPARTIAL_MATCH\t1abc\t0.6667\t1.2500\tA10,A20,A30"
        );
        // --skip-match: the RMSD was never computed, so it is NA and not 0.0000
        assert_eq!(
            novelty_verdict_line("design.pdb", residues, 3, 0.8, 2.0, Some(("1abc", 3, None))),
            "design.pdb\tKNOWN\t1abc\t1.0000\tNA\tA10,A20,A30"
        );
        // Zero residues must not divide by zero, and is never KNOWN
        assert_eq!(
            novelty_verdict_line("design.pdb", "", 0, 0.0, 2.0, Some(("1abc", 0, None))),
            "design.pdb\tNOVEL\t1abc\t0.0000\tNA\t"
        );
    }

    #[test]
    fn test_novelty_verdict_needs_a_close_hit_not_just_a_covering_one() {
        // The measured failure mode: an 8-residue motif covered 1.0 by a hit 6.15 A
        // away was called KNOWN. Covering every residue in a different arrangement is
        // a different motif.
        let residues = "A10,A20,A30,A40,A50,A60,A70,A80";
        assert_eq!(
            novelty_verdict_line("design.pdb", residues, 8, 0.8, 2.0, Some(("1abc", 8, Some(6.1487)))),
            format!("design.pdb\tPARTIAL_MATCH\t1abc\t1.0000\t6.1487\t{}", residues)
        );
        // The genuine hits of the same measurement stay KNOWN, even at 0.5 A
        for (rmsd, threshold, tier) in [
            (0.0599f32, 2.0f32, "KNOWN"), (0.1669, 2.0, "KNOWN"),
            (0.0599, 0.5, "KNOWN"), (4.3167, 2.0, "PARTIAL_MATCH"),
        ] {
            let line = novelty_verdict_line("d.pdb", "A1,A2,A3", 3, 0.8, threshold,
                Some(("1abc", 3, Some(rmsd))));
            assert_eq!(line.split('\t').nth(1).unwrap(), tier, "{} A at limit {}", rmsd, threshold);
        }
        // Without a computed RMSD (--skip-match) the gate cannot apply
        let line = novelty_verdict_line("d.pdb", "A1,A2,A3", 3, 0.8, 0.5, Some(("1abc", 3, None)));
        assert_eq!(line.split('\t').nth(1).unwrap(), "KNOWN");
    }

    #[test]
    fn test_no_hash_query_is_not_reported_as_novel() {
        // Residues further apart than the index cutoff hash to nothing. Reporting that
        // as NOVEL turns unsearchable input into a discovery.
        assert_eq!(
            novelty_no_hash_line("query/1SU6.pdb", "A4,A39,A70"),
            "query/1SU6.pdb\tNO_HASHES\tNA\t0.0000\tNA\tA4,A39,A70"
        );
    }
}
