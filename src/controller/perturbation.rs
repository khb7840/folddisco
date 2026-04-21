use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};

use nalgebra::{Isometry3, Point3, Translation3, Unit, UnitQuaternion, Vector3};
use plotters::prelude::*;
use rand::Rng;

use crate::controller::feature::get_single_feature;
use crate::controller::io::read_structure_from_path;
use crate::controller::query::parse_query_string;
use crate::geometry::core::{GeometricHash, HashType};
use crate::structure::coordinate::Coordinate;
use crate::structure::core::CompactStructure;

const FEATURE_NAMES: [&str; 7] = [
    "ca_distance",
    "cb_distance",
    "ca_cb_angle",
    "theta1",
    "theta2",
    "phi1",
    "phi2",
];
const MAX_DISCRETIZED_ENCODING_BINS: usize = 30;
const HISTOGRAM_BIN_COUNT: usize = 30;

#[derive(Debug, Clone, Copy)]
pub struct ResidueBackbone {
    pub index: usize,
    pub chain: u8,
    pub residue_serial: u64,
    pub residue_name: [u8; 3],
    pub n: Coordinate,
    pub ca: Coordinate,
    pub cb: Coordinate,
}

impl ResidueBackbone {
    pub fn from_compact(structure: &CompactStructure, index: usize) -> Option<Self> {
        Some(Self {
            index,
            chain: structure.chain_per_residue.get(index).copied()?,
            residue_serial: *structure.residue_serial.get(index)?,
            residue_name: *structure.residue_name.get(index)?,
            n: structure.get_n(index)?,
            ca: structure.get_ca(index)?,
            cb: structure.get_cb(index)?,
        })
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PerturbationRecord {
    pub original_residue: ResidueBackbone,
    pub perturbed_residue: ResidueBackbone,
    pub translation_magnitude_angstrom: f32,
    pub rotation_magnitude_degrees: f32,
    pub translation_vector: [f32; 3],
    pub rotation_vector: [f32; 3],
}

pub struct Perturbator {
    pub translation_range_angstrom: (f32, f32),
    pub rotation_range_degrees: (f32, f32),
}

impl Default for Perturbator {
    fn default() -> Self {
        Self {
            translation_range_angstrom: (0.1, 1.5),
            rotation_range_degrees: (1.0, 45.0),
        }
    }
}

impl Perturbator {
    pub fn new(translation_range_angstrom: (f32, f32), rotation_range_degrees: (f32, f32)) -> Self {
        Self {
            translation_range_angstrom,
            rotation_range_degrees,
        }
    }

    pub fn perturb_residue(&self, residue: &ResidueBackbone) -> PerturbationRecord {
        let mut rng = rand::thread_rng();
        let translation_magnitude =
            rng.gen_range(self.translation_range_angstrom.0..=self.translation_range_angstrom.1);
        let rotation_magnitude =
            rng.gen_range(self.rotation_range_degrees.0..=self.rotation_range_degrees.1);
        let translation_dir = random_unit_vector(&mut rng);
        let rotation_axis = random_unit_vector(&mut rng);
        let translation = translation_dir * translation_magnitude;
        self.perturb_with_transform(residue, translation, rotation_axis, rotation_magnitude)
    }

    pub fn perturb_with_transform(
        &self,
        residue: &ResidueBackbone,
        translation: Vector3<f32>,
        rotation_axis: Vector3<f32>,
        rotation_degrees: f32,
    ) -> PerturbationRecord {
        let axis_unit = Unit::new_normalize(if rotation_axis.norm_squared() == 0.0 {
            Vector3::new(1.0, 0.0, 0.0)
        } else {
            rotation_axis
        });
        let rotation = UnitQuaternion::from_axis_angle(&axis_unit, rotation_degrees.to_radians());
        let isometry = Isometry3::from_parts(
            Translation3::new(translation.x, translation.y, translation.z),
            rotation,
        );
        let perturbed = transform_residue(residue, &isometry);
        let rotation_vec = axis_unit.into_inner() * rotation_degrees;
        PerturbationRecord {
            original_residue: *residue,
            perturbed_residue: perturbed,
            translation_magnitude_angstrom: translation.norm(),
            rotation_magnitude_degrees: rotation_degrees.abs(),
            translation_vector: [translation.x, translation.y, translation.z],
            rotation_vector: [rotation_vec.x, rotation_vec.y, rotation_vec.z],
        }
    }
}

#[derive(Debug, Clone)]
pub struct FeatureDrift {
    pub residue_pair: (usize, usize),
    pub translation_distances_angstrom: [f32; 2],
    pub rotation_angles_degrees: [f32; 2],
    pub initial_raw_features: [f32; 7],
    pub perturbed_raw_features: [f32; 7],
    pub initial_discretized_encoding: u32,
    pub perturbed_discretized_encoding: u32,
    pub absolute_delta: [f32; 7],
}

pub fn compute_feature_drift(
    structure: &CompactStructure,
    residue_pair: (usize, usize),
    perturbed_first: &PerturbationRecord,
    perturbed_second: &PerturbationRecord,
    nbin_dist: usize,
    nbin_angle: usize,
) -> Option<FeatureDrift> {
    let mut initial_feature = vec![0.0; 9];
    if !get_single_feature(
        residue_pair.0,
        residue_pair.1,
        structure,
        HashType::Hybrid,
        1000.0,
        &mut initial_feature,
    ) {
        return None;
    }

    let mut perturbed_structure = structure.clone();
    apply_record_to_structure(&mut perturbed_structure, perturbed_first);
    apply_record_to_structure(&mut perturbed_structure, perturbed_second);

    let mut perturbed_feature = vec![0.0; 9];
    if !get_single_feature(
        residue_pair.0,
        residue_pair.1,
        &perturbed_structure,
        HashType::Hybrid,
        1000.0,
        &mut perturbed_feature,
    ) {
        return None;
    }

    let initial_raw = feature7_from_full(&initial_feature);
    let perturbed_raw = feature7_from_full(&perturbed_feature);
    let mut absolute_delta = [0.0; 7];
    for i in 0..7 {
        absolute_delta[i] = (perturbed_raw[i] - initial_raw[i]).abs();
    }

    let initial_encoding = if nbin_dist == 0 || nbin_angle == 0 {
        GeometricHash::perfect_hash_default_as_u32(&initial_feature, HashType::Hybrid)
    } else {
        GeometricHash::perfect_hash_as_u32(
            &initial_feature,
            HashType::Hybrid,
            nbin_dist,
            nbin_angle,
        )
    };
    let perturbed_encoding = if nbin_dist == 0 || nbin_angle == 0 {
        GeometricHash::perfect_hash_default_as_u32(&perturbed_feature, HashType::Hybrid)
    } else {
        GeometricHash::perfect_hash_as_u32(
            &perturbed_feature,
            HashType::Hybrid,
            nbin_dist,
            nbin_angle,
        )
    };

    Some(FeatureDrift {
        residue_pair,
        translation_distances_angstrom: [
            perturbed_first.translation_magnitude_angstrom,
            perturbed_second.translation_magnitude_angstrom,
        ],
        rotation_angles_degrees: [
            perturbed_first.rotation_magnitude_degrees,
            perturbed_second.rotation_magnitude_degrees,
        ],
        initial_raw_features: initial_raw,
        perturbed_raw_features: perturbed_raw,
        initial_discretized_encoding: initial_encoding,
        perturbed_discretized_encoding: perturbed_encoding,
        absolute_delta,
    })
}

pub fn collect_feature_drift_for_structure(
    structure: &CompactStructure,
    residue_indices: &[usize],
    perturbator: &Perturbator,
    max_pairs: usize,
    nbin_dist: usize,
    nbin_angle: usize,
) -> Vec<FeatureDrift> {
    let mut output = Vec::new();
    let mut tried = 0usize;
    for i in 0..residue_indices.len() {
        for j in (i + 1)..residue_indices.len() {
            if tried >= max_pairs {
                return output;
            }
            let idx1 = residue_indices[i];
            let idx2 = residue_indices[j];
            tried += 1;
            let res1 = if let Some(v) = ResidueBackbone::from_compact(structure, idx1) {
                v
            } else {
                continue;
            };
            let res2 = if let Some(v) = ResidueBackbone::from_compact(structure, idx2) {
                v
            } else {
                continue;
            };
            let record1 = perturbator.perturb_residue(&res1);
            let record2 = perturbator.perturb_residue(&res2);
            if let Some(drift) = compute_feature_drift(
                structure,
                (idx1, idx2),
                &record1,
                &record2,
                nbin_dist,
                nbin_angle,
            ) {
                output.push(drift);
            }
        }
    }
    output
}

#[derive(Debug, Default, Clone)]
pub struct QueryDriftDatasets {
    pub whole_structure: Vec<FeatureDrift>,
    pub motif_only: Vec<FeatureDrift>,
}

pub fn collect_query_directory_feature_drifts(
    query_dir: &Path,
    perturbator: &Perturbator,
    max_pairs_per_structure: usize,
    nbin_dist: usize,
    nbin_angle: usize,
) -> Result<QueryDriftDatasets, Box<dyn Error>> {
    let mut datasets = QueryDriftDatasets::default();
    let entries = fs::read_dir(query_dir)?;

    for entry in entries {
        let path = entry?.path();
        if is_structure_path(&path) {
            let Some(path_str) = path.to_str() else {
                continue;
            };
            let Some(structure) = read_structure_from_path(path_str) else {
                continue;
            };
            let compact = structure.to_compact();
            if compact.num_residues < 5 {
                continue;
            }
            let whole_indices: Vec<usize> = (1..(compact.num_residues - 1)).collect();
            datasets
                .whole_structure
                .extend(collect_feature_drift_for_structure(
                    &compact,
                    &whole_indices,
                    perturbator,
                    max_pairs_per_structure,
                    nbin_dist,
                    nbin_angle,
                ));
        } else if path.extension().and_then(|s| s.to_str()) == Some("txt") {
            let text = fs::read_to_string(&path)?;
            for line in text.lines() {
                let cols: Vec<&str> = line.split('\t').collect();
                if cols.len() < 2 || cols[0].trim().is_empty() {
                    continue;
                }
                let structure_path = resolve_structure_path(query_dir, cols[0].trim());
                let Some(structure_path_str) = structure_path.to_str() else {
                    continue;
                };
                let Some(structure) = read_structure_from_path(structure_path_str) else {
                    continue;
                };
                let compact = structure.to_compact();
                if compact.num_residues < 5 {
                    continue;
                }
                let (query_residues, _amino_acid_substitutions) =
                    parse_query_string(cols[1].trim(), b'A');
                let motif_indices: Vec<usize> = query_residues
                    .iter()
                    .filter_map(|(chain, serial)| compact.get_index(chain, serial))
                    .filter(|i| *i > 0 && *i < compact.num_residues - 1)
                    .collect();
                if motif_indices.len() < 2 {
                    continue;
                }
                datasets
                    .motif_only
                    .extend(collect_feature_drift_for_structure(
                        &compact,
                        &motif_indices,
                        perturbator,
                        max_pairs_per_structure,
                        nbin_dist,
                        nbin_angle,
                    ));
            }
        }
    }
    Ok(datasets)
}

pub fn generate_feature_drift_vs_transform_plots(
    records: &[FeatureDrift],
    output_dir: &Path,
) -> Result<Vec<PathBuf>, Box<dyn Error>> {
    fs::create_dir_all(output_dir)?;
    let mut outputs = Vec::new();
    if records.is_empty() {
        return Ok(outputs);
    }
    for feature_idx in 0..7 {
        let path = output_dir.join(format!(
            "drift_vs_translation_{}.png",
            FEATURE_NAMES[feature_idx]
        ));
        scatter_plot(
            &path,
            &format!("Δ{} vs Translation", FEATURE_NAMES[feature_idx]),
            "Translation (Å)",
            records
                .iter()
                .map(|r| {
                    (
                        (r.translation_distances_angstrom[0] + r.translation_distances_angstrom[1])
                            * 0.5,
                        r.absolute_delta[feature_idx],
                    )
                })
                .collect::<Vec<_>>()
                .as_slice(),
        )?;
        outputs.push(path);
    }
    Ok(outputs)
}

pub fn generate_raw_distribution_plots(
    records: &[FeatureDrift],
    output_dir: &Path,
) -> Result<Vec<PathBuf>, Box<dyn Error>> {
    fs::create_dir_all(output_dir)?;
    let mut outputs = Vec::new();
    if records.is_empty() {
        return Ok(outputs);
    }
    for feature_idx in 0..7 {
        let mut init = Vec::with_capacity(records.len());
        let mut pert = Vec::with_capacity(records.len());
        for record in records {
            init.push(record.initial_raw_features[feature_idx]);
            pert.push(record.perturbed_raw_features[feature_idx]);
        }
        let path = output_dir.join(format!(
            "raw_distribution_{}.png",
            FEATURE_NAMES[feature_idx]
        ));
        histogram_overlay_plot(
            &path,
            &format!("Raw Distribution {}", FEATURE_NAMES[feature_idx]),
            &init,
            &pert,
        )?;
        outputs.push(path);
    }
    Ok(outputs)
}

pub fn generate_discretized_distribution_plot(
    records: &[FeatureDrift],
    output_path: &Path,
) -> Result<(), Box<dyn Error>> {
    let mut init_counts = HashMap::<u32, u32>::new();
    let mut pert_counts = HashMap::<u32, u32>::new();
    for r in records {
        *init_counts
            .entry(r.initial_discretized_encoding)
            .or_insert(0) += 1;
        *pert_counts
            .entry(r.perturbed_discretized_encoding)
            .or_insert(0) += 1;
    }
    let mut keys: Vec<u32> = init_counts
        .keys()
        .chain(pert_counts.keys())
        .copied()
        .collect();
    keys.sort_unstable();
    keys.dedup();
    if keys.len() > MAX_DISCRETIZED_ENCODING_BINS {
        keys.truncate(MAX_DISCRETIZED_ENCODING_BINS);
    }
    let max_count = keys
        .iter()
        .map(|k| {
            let a = *init_counts.get(k).unwrap_or(&0);
            let b = *pert_counts.get(k).unwrap_or(&0);
            a.max(b)
        })
        .max()
        .unwrap_or(1);

    let root = BitMapBackend::new(output_path, (1400, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption("Discretized Encoding Distribution", ("sans-serif", 30))
        .x_label_area_size(50)
        .y_label_area_size(50)
        .build_cartesian_2d(0f32..(keys.len() as f32), 0u32..(max_count + 1))?;
    chart
        .configure_mesh()
        .x_desc("Encoding rank")
        .y_desc("Count")
        .draw()?;
    for (i, key) in keys.iter().enumerate() {
        let x = i as f32;
        let a = *init_counts.get(key).unwrap_or(&0);
        let b = *pert_counts.get(key).unwrap_or(&0);
        chart.draw_series(std::iter::once(Rectangle::new(
            [(x - 0.30, 0), (x - 0.05, a)],
            BLUE.mix(0.5).filled(),
        )))?;
        chart.draw_series(std::iter::once(Rectangle::new(
            [(x + 0.05, 0), (x + 0.30, b)],
            RED.mix(0.5).filled(),
        )))?;
    }
    root.present()?;
    Ok(())
}

pub struct ReportGenerator;

impl ReportGenerator {
    pub fn generate_markdown_report(
        output_path: &Path,
        records: &[FeatureDrift],
        plot_paths: &[PathBuf],
    ) -> Result<(), Box<dyn Error>> {
        let mut report = String::new();
        report.push_str("# Feature Drift Report\n\n");
        report.push_str(&format!("Total records: {}\n\n", records.len()));
        report.push_str("## Overall Drift Summary\n\n");
        report.push_str("| Feature | Mean Δ | Max Δ |\n|---|---:|---:|\n");
        for (idx, name) in FEATURE_NAMES.iter().enumerate() {
            let (mean, max) = summarize_feature_delta(records, idx);
            report.push_str(&format!("| {} | {:.6} | {:.6} |\n", name, mean, max));
        }

        report.push_str("\n## Grouped by Translation Magnitude (Å)\n\n");
        for (label, filtered) in grouped_by_translation(records) {
            report.push_str(&format!("### {}\n\n", label));
            report.push_str("| Feature | Mean Δ | Max Δ |\n|---|---:|---:|\n");
            for (idx, name) in FEATURE_NAMES.iter().enumerate() {
                let (mean, max) = summarize_feature_delta(&filtered, idx);
                report.push_str(&format!("| {} | {:.6} | {:.6} |\n", name, mean, max));
            }
            report.push('\n');
        }

        report.push_str("## Grouped by Rotation Magnitude (degrees)\n\n");
        for (label, filtered) in grouped_by_rotation(records) {
            report.push_str(&format!("### {}\n\n", label));
            report.push_str("| Feature | Mean Δ | Max Δ |\n|---|---:|---:|\n");
            for (idx, name) in FEATURE_NAMES.iter().enumerate() {
                let (mean, max) = summarize_feature_delta(&filtered, idx);
                report.push_str(&format!("| {} | {:.6} | {:.6} |\n", name, mean, max));
            }
            report.push('\n');
        }

        report.push_str("## Generated Plots\n\n");
        for path in plot_paths {
            report.push_str(&format!("- {}\n", path.display()));
        }
        fs::write(output_path, report)?;
        Ok(())
    }
}

fn grouped_by_translation(records: &[FeatureDrift]) -> Vec<(String, Vec<FeatureDrift>)> {
    let bins = [(0.0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, f32::INFINITY)];
    bins.iter()
        .map(|(lo, hi)| {
            let label = if hi.is_infinite() {
                format!("Translation [{:.1}, inf)", lo)
            } else {
                format!("Translation [{:.1}, {:.1})", lo, hi)
            };
            let selected = records
                .iter()
                .filter(|r| {
                    let t = (r.translation_distances_angstrom[0]
                        + r.translation_distances_angstrom[1])
                        * 0.5;
                    t >= *lo && t < *hi
                })
                .cloned()
                .collect();
            (label, selected)
        })
        .collect()
}

fn grouped_by_rotation(records: &[FeatureDrift]) -> Vec<(String, Vec<FeatureDrift>)> {
    let bins = [
        (0.0, 10.0),
        (10.0, 20.0),
        (20.0, 30.0),
        (30.0, 45.0),
        (45.0, f32::INFINITY),
    ];
    bins.iter()
        .map(|(lo, hi)| {
            let label = if hi.is_infinite() {
                format!("Rotation [{:.1}, inf)", lo)
            } else {
                format!("Rotation [{:.1}, {:.1})", lo, hi)
            };
            let selected = records
                .iter()
                .filter(|r| {
                    let rot = (r.rotation_angles_degrees[0] + r.rotation_angles_degrees[1]) * 0.5;
                    rot >= *lo && rot < *hi
                })
                .cloned()
                .collect();
            (label, selected)
        })
        .collect()
}

fn summarize_feature_delta(records: &[FeatureDrift], idx: usize) -> (f32, f32) {
    if records.is_empty() {
        return (0.0, 0.0);
    }
    let mut sum = 0.0f32;
    let mut max = 0.0f32;
    for record in records {
        let val = record.absolute_delta[idx];
        sum += val;
        if val > max {
            max = val;
        }
    }
    (sum / records.len() as f32, max)
}

fn feature7_from_full(feature: &[f32]) -> [f32; 7] {
    [
        feature[2], feature[3], feature[4], feature[5], feature[6], feature[7], feature[8],
    ]
}

fn apply_record_to_structure(structure: &mut CompactStructure, record: &PerturbationRecord) {
    let idx = record.perturbed_residue.index;
    let n = record.perturbed_residue.n;
    let ca = record.perturbed_residue.ca;
    let cb = record.perturbed_residue.cb;
    if idx < structure.n_vector.x.len() {
        structure.n_vector.x[idx] = Some(n.x);
        structure.n_vector.y[idx] = Some(n.y);
        structure.n_vector.z[idx] = Some(n.z);
        structure.ca_vector.x[idx] = Some(ca.x);
        structure.ca_vector.y[idx] = Some(ca.y);
        structure.ca_vector.z[idx] = Some(ca.z);
        structure.cb_vector.x[idx] = Some(cb.x);
        structure.cb_vector.y[idx] = Some(cb.y);
        structure.cb_vector.z[idx] = Some(cb.z);
    }
}

fn transform_residue(residue: &ResidueBackbone, isometry: &Isometry3<f32>) -> ResidueBackbone {
    let n = transform_coordinate(residue.n, isometry);
    let ca = transform_coordinate(residue.ca, isometry);
    let cb = transform_coordinate(residue.cb, isometry);
    ResidueBackbone {
        n,
        ca,
        cb,
        ..*residue
    }
}

fn transform_coordinate(c: Coordinate, iso: &Isometry3<f32>) -> Coordinate {
    let p = iso.transform_point(&Point3::new(c.x, c.y, c.z));
    Coordinate::new(p.x, p.y, p.z)
}

fn random_unit_vector<R: Rng + ?Sized>(rng: &mut R) -> Vector3<f32> {
    loop {
        let v = Vector3::new(
            rng.gen_range(-1.0..=1.0),
            rng.gen_range(-1.0..=1.0),
            rng.gen_range(-1.0..=1.0),
        );
        if v.norm_squared() > 1e-8 {
            return v.normalize();
        }
    }
}

fn is_structure_path(path: &Path) -> bool {
    let ext = path.extension().and_then(|s| s.to_str()).unwrap_or("");
    matches!(ext, "pdb" | "ent" | "cif")
}

fn resolve_structure_path(query_dir: &Path, rel_or_abs: &str) -> PathBuf {
    let candidate = PathBuf::from(rel_or_abs);
    if candidate.is_absolute() {
        candidate
    } else {
        let direct = query_dir.join(rel_or_abs);
        if direct.exists() {
            direct
        } else {
            query_dir
                .parent()
                .map(|p| p.join(rel_or_abs))
                .unwrap_or(direct)
        }
    }
}

fn scatter_plot(
    output_path: &Path,
    title: &str,
    x_desc: &str,
    points: &[(f32, f32)],
) -> Result<(), Box<dyn Error>> {
    if points.is_empty() {
        return Ok(());
    }
    let x_max = points
        .iter()
        .map(|(x, _)| *x)
        .fold(0.0_f32, f32::max)
        .max(1.0);
    let y_max = points
        .iter()
        .map(|(_, y)| *y)
        .fold(0.0_f32, f32::max)
        .max(1.0);
    let root = BitMapBackend::new(output_path, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption(title, ("sans-serif", 30))
        .x_label_area_size(45)
        .y_label_area_size(60)
        .build_cartesian_2d(0f32..x_max, 0f32..y_max)?;
    chart
        .configure_mesh()
        .x_desc(x_desc)
        .y_desc("Absolute Δ")
        .draw()?;
    chart.draw_series(
        points
            .iter()
            .map(|(x, y)| Circle::new((*x, *y), 3, RED.filled())),
    )?;
    root.present()?;
    Ok(())
}

fn histogram_overlay_plot(
    output_path: &Path,
    title: &str,
    initial: &[f32],
    perturbed: &[f32],
) -> Result<(), Box<dyn Error>> {
    if initial.is_empty() && perturbed.is_empty() {
        return Ok(());
    }
    let mut all = Vec::with_capacity(initial.len() + perturbed.len());
    all.extend_from_slice(initial);
    all.extend_from_slice(perturbed);
    let min_v = all
        .iter()
        .copied()
        .fold(f32::INFINITY, |a, b| if b < a { b } else { a });
    let max_v = all
        .iter()
        .copied()
        .fold(f32::NEG_INFINITY, |a, b| if b > a { b } else { a });
    let bins = HISTOGRAM_BIN_COUNT;
    let span = (max_v - min_v).max(1e-6);
    let bin_width = span / bins as f32;
    let mut init_counts = vec![0u32; bins];
    let mut pert_counts = vec![0u32; bins];
    for &v in initial {
        let idx = (((v - min_v) / bin_width).floor() as usize).min(bins - 1);
        init_counts[idx] += 1;
    }
    for &v in perturbed {
        let idx = (((v - min_v) / bin_width).floor() as usize).min(bins - 1);
        pert_counts[idx] += 1;
    }
    let max_count = init_counts
        .iter()
        .chain(pert_counts.iter())
        .copied()
        .max()
        .unwrap_or(1);

    let root = BitMapBackend::new(output_path, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption(title, ("sans-serif", 30))
        .x_label_area_size(45)
        .y_label_area_size(60)
        .build_cartesian_2d(min_v..max_v, 0u32..(max_count + 1))?;
    chart
        .configure_mesh()
        .x_desc("Feature value")
        .y_desc("Count")
        .draw()?;

    chart.draw_series((0..bins).map(|i| {
        let x0 = min_v + i as f32 * bin_width;
        let x1 = x0 + bin_width * 0.45;
        Rectangle::new([(x0, 0), (x1, init_counts[i])], BLUE.mix(0.5).filled())
    }))?;
    chart.draw_series((0..bins).map(|i| {
        let x0 = min_v + i as f32 * bin_width + (bin_width * 0.5);
        let x1 = x0 + bin_width * 0.45;
        Rectangle::new([(x0, 0), (x1, pert_counts[i])], RED.mix(0.5).filled())
    }))?;
    root.present()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::coordinate::CarbonCoordinateVector;

    #[test]
    fn test_rigid_body_transformation_preserves_internal_distances() {
        let residue = ResidueBackbone {
            index: 1,
            chain: b'A',
            residue_serial: 10,
            residue_name: *b"ALA",
            n: Coordinate::new(0.0, 0.0, 0.0),
            ca: Coordinate::new(1.0, 0.0, 0.0),
            cb: Coordinate::new(1.0, 1.0, 0.0),
        };
        let perturbator = Perturbator::default();
        let record = perturbator.perturb_with_transform(
            &residue,
            Vector3::new(1.2, -0.8, 0.5),
            Vector3::new(0.0, 1.0, 0.0),
            27.0,
        );
        let eps = 1e-4f32;
        let d1_before = residue.n.calc_distance(&residue.ca);
        let d2_before = residue.ca.calc_distance(&residue.cb);
        let d3_before = residue.n.calc_distance(&residue.cb);
        let d1_after = record
            .perturbed_residue
            .n
            .calc_distance(&record.perturbed_residue.ca);
        let d2_after = record
            .perturbed_residue
            .ca
            .calc_distance(&record.perturbed_residue.cb);
        let d3_after = record
            .perturbed_residue
            .n
            .calc_distance(&record.perturbed_residue.cb);
        assert!((d1_before - d1_after).abs() < eps);
        assert!((d2_before - d2_after).abs() < eps);
        assert!((d3_before - d3_after).abs() < eps);
    }

    #[test]
    fn test_known_translation_vector_delta_expectation() {
        let structure = build_test_structure();
        let perturbator = Perturbator::default();
        let r1 = ResidueBackbone::from_compact(&structure, 1).unwrap();
        let r3 = ResidueBackbone::from_compact(&structure, 3).unwrap();
        let translation = Vector3::new(1.0, 0.0, 0.0);
        let rec1 =
            perturbator.perturb_with_transform(&r1, translation, Vector3::new(1.0, 0.0, 0.0), 0.0);
        let rec2 =
            perturbator.perturb_with_transform(&r3, translation, Vector3::new(1.0, 0.0, 0.0), 0.0);
        let drift = compute_feature_drift(&structure, (1, 3), &rec1, &rec2, 0, 0).unwrap();
        for delta in drift.absolute_delta {
            assert!(delta.abs() < 1e-4);
        }
    }

    fn build_test_structure() -> CompactStructure {
        let num_residues = 5usize;
        let mut n = CarbonCoordinateVector::with_capacity(num_residues);
        let mut ca = CarbonCoordinateVector::with_capacity(num_residues);
        let mut cb = CarbonCoordinateVector::with_capacity(num_residues);
        for i in 0..num_residues {
            let x = i as f32 * 3.8;
            n.push(&Coordinate::new(x - 0.5, 0.0, 0.0));
            ca.push(&Coordinate::new(x, 0.5, 0.0));
            cb.push(&Coordinate::new(x + 0.2, 1.4, 0.6));
        }
        CompactStructure {
            num_chains: 1,
            chains: vec![b'A'],
            chain_per_residue: vec![b'A'; num_residues],
            num_residues,
            residue_serial: (1..=num_residues as u64).collect(),
            residue_name: vec![*b"ALA"; num_residues],
            n_vector: n,
            ca_vector: ca,
            cb_vector: cb,
            b_factors: vec![50.0; num_residues],
        }
    }
}
