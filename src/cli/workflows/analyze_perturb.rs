// File: analyze_perturb.rs
// Created: 2026-04-21
// Author: Hyunbin Kim (khb7840@gmail.com)

use std::error::Error;
use std::path::{Path, PathBuf};

use crate::cli::{print_logo, AppArgs};
use crate::controller::perturbation::{
    collect_query_directory_feature_drifts, generate_discretized_distribution_plot,
    generate_feature_drift_vs_transform_plots, generate_raw_distribution_plots, FeatureDrift,
    Perturbator, ReportGenerator,
};
use crate::prelude::*;

pub const HELP_ANALYZE_PERTURB: &str = "\
usage: folddisco analyze-perturb --query-dir <DIR> [options]

input:
    --query-dir <DIR>              Query directory containing structure files and/or motif TXT files

output:
    -o, --output <DIR>             Output directory prefix [default: <query-dir>/perturbation_analysis]
    --dataset <STR>                Dataset to export (whole|motif|both) [default: both]

analysis options:
    --max-pairs <INT>              Max residue pairs per structure [default: 200]
    -d, --distance <INT>           Distance bins for discretized encoding [default: 0]
    -a, --angle <INT>              Angle bins for discretized encoding [default: 0]
    --translation-min <FLOAT>      Minimum translation perturbation in Å [default: 0.1]
    --translation-max <FLOAT>      Maximum translation perturbation in Å [default: 1.5]
    --rotation-min <FLOAT>         Minimum rotation perturbation in degrees [default: 1.0]
    --rotation-max <FLOAT>         Maximum rotation perturbation in degrees [default: 45.0]

general options:
    -v, --verbose                  Print verbose messages
    -h, --help                     Print this help menu
";

pub fn analyze_perturb(env: AppArgs) {
    match env {
        AppArgs::AnalyzePerturb {
            query_dir,
            output,
            max_pairs,
            num_bin_dist,
            num_bin_angle,
            translation_min,
            translation_max,
            rotation_min,
            rotation_max,
            dataset,
            verbose,
            help: _,
        } => {
            if verbose {
                print_logo();
            }
            if query_dir.is_empty() {
                eprintln!("{}", HELP_ANALYZE_PERTURB);
                std::process::exit(1);
            }
            if translation_min > translation_max || rotation_min > rotation_max {
                print_log_msg(
                    FAIL,
                    "Invalid perturbation range: min value cannot be larger than max value",
                );
                std::process::exit(1);
            }

            let query_path = Path::new(&query_dir);
            if !query_path.is_dir() {
                print_log_msg(
                    FAIL,
                    &format!("Query directory does not exist: {}", query_dir),
                );
                std::process::exit(1);
            }

            let selection = dataset.trim().to_ascii_lowercase();
            if !matches!(selection.as_str(), "whole" | "motif" | "both") {
                print_log_msg(
                    FAIL,
                    &format!(
                        "Invalid --dataset value: {} (allowed: whole, motif, both)",
                        dataset
                    ),
                );
                std::process::exit(1);
            }

            let output_root = if output.is_empty() {
                query_path.join("perturbation_analysis")
            } else {
                PathBuf::from(output)
            };
            std::fs::create_dir_all(&output_root).expect(&log_msg(
                FAIL,
                &format!(
                    "Failed to create output directory: {}",
                    output_root.display()
                ),
            ));
            if verbose {
                print_log_msg(
                    INFO,
                    &format!("Analyzing perturbation drift from {}", query_path.display()),
                );
            }

            let perturbator = Perturbator::new(
                (translation_min, translation_max),
                (rotation_min, rotation_max),
            );
            let datasets = collect_query_directory_feature_drifts(
                query_path,
                &perturbator,
                max_pairs,
                num_bin_dist,
                num_bin_angle,
            )
            .expect(&log_msg(
                FAIL,
                "Failed to collect perturbation drift datasets",
            ));

            if selection == "whole" || selection == "both" {
                write_dataset_outputs(
                    "whole_structure",
                    &datasets.whole_structure,
                    &output_root,
                    verbose,
                )
                .expect(&log_msg(
                    FAIL,
                    "Failed to generate perturbation outputs for whole_structure",
                ));
            }
            if selection == "motif" || selection == "both" {
                write_dataset_outputs("motif_only", &datasets.motif_only, &output_root, verbose)
                    .expect(&log_msg(
                        FAIL,
                        "Failed to generate perturbation outputs for motif_only",
                    ));
            }
            if verbose {
                print_log_msg(
                    DONE,
                    &format!(
                        "Perturbation analysis completed at {}",
                        output_root.display()
                    ),
                );
            }
        }
        _ => {
            eprintln!("{}", HELP_ANALYZE_PERTURB);
            std::process::exit(1);
        }
    }
}

fn write_dataset_outputs(
    dataset_name: &str,
    records: &[FeatureDrift],
    output_root: &Path,
    verbose: bool,
) -> Result<(), Box<dyn Error>> {
    let dataset_dir = output_root.join(dataset_name);
    std::fs::create_dir_all(&dataset_dir)?;
    let plots_dir = dataset_dir.join("plots");
    std::fs::create_dir_all(&plots_dir)?;
    let mut plot_paths = Vec::new();

    if !records.is_empty() {
        match generate_feature_drift_vs_transform_plots(records, &plots_dir) {
            Ok(paths) => plot_paths.extend(paths),
            Err(e) => print_log_msg(
                WARN,
                &format!(
                    "Failed to generate drift-vs-transform plots for {}: {}",
                    dataset_name, e
                ),
            ),
        }
        match generate_raw_distribution_plots(records, &plots_dir) {
            Ok(paths) => plot_paths.extend(paths),
            Err(e) => print_log_msg(
                WARN,
                &format!(
                    "Failed to generate raw distribution plots for {}: {}",
                    dataset_name, e
                ),
            ),
        }
        let discretized_plot = plots_dir.join("discretized_distribution.png");
        match generate_discretized_distribution_plot(records, &discretized_plot) {
            Ok(()) => plot_paths.push(discretized_plot),
            Err(e) => print_log_msg(
                WARN,
                &format!(
                    "Failed to generate discretized distribution plot for {}: {}",
                    dataset_name, e
                ),
            ),
        }
    }

    let report_path = dataset_dir.join("feature_drift_report.md");
    ReportGenerator::generate_markdown_report(&report_path, records, &plot_paths)?;
    if verbose {
        print_log_msg(
            INFO,
            &format!(
                "Saved {} records for {} to {}",
                records.len(),
                dataset_name,
                dataset_dir.display()
            ),
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_temp_dir(test_name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "folddisco_{}_{}_{}",
            test_name,
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn test_write_dataset_outputs_generates_report_for_empty_records() {
        let out_dir = make_temp_dir("analyze_perturb_empty");
        let records = vec![];
        write_dataset_outputs("motif_only", &records, &out_dir, false).unwrap();
        assert!(out_dir
            .join("motif_only")
            .join("feature_drift_report.md")
            .exists());
        assert!(!out_dir
            .join("motif_only")
            .join("plots")
            .join("discretized_distribution.png")
            .exists());
        let report =
            std::fs::read_to_string(out_dir.join("motif_only").join("feature_drift_report.md"))
                .unwrap();
        assert!(report.contains("Total records: 0"));
    }
}
