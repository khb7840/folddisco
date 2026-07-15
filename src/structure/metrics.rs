// File: metrics.rs
// Created: 2025-10-29
// Description: Structure similarity metrics implementation
//   - TM-score: Template Modeling score
//   - GDT-TS: Global Distance Test - Total Score.
//   - GDT-HA: Global Distance Test - High Accuracy
//   - Chamfer Distance: Average of nearest neighbor distances
//   - Hausdorff Distance: Maximum of minimum distances. Outlier sensitive.
//   - RMSD: Root Mean Square Deviation. Outlier sensitive.
//
// Performance optimizations:
//   - Uses f64 for internal arithmetic to maintain precision
//   - PrecomputedDistances struct for efficient batch metric calculation
//   - Avoids redundant distance calculations across metrics
//
// Usage example:
// ```rust
// use crate::structure::metrics::*;
//
// // Method 1: Calculate all metrics (simple but slower)
// let metrics = StructureMetrics::calculate_all(&ref_coords, &model_coords);
// println!("TM-score: {:.4}", metrics.tm_score);
//
// // Method 2: Calculate all metrics fast (2x faster, recommended)
// let metrics = StructureMetrics::calculate_all_fast(&ref_coords, &model_coords);
// metrics.print();
//
// // Method 3: Pre-compute distances and calculate individual metrics
// let distances = PrecomputedDistances::new(&ref_coords, &model_coords);
// let tm = tm_score_fast(&distances, None);
// let gdt = gdt_ts_fast(&distances);
// let lddt = lddt_default_fast(&distances);
// let rmsd = rmsd_fast(&distances);
// ```

use core::fmt;

/// Pre-computed distances for efficient metric calculation
/// 
/// This struct stores pre-calculated distances to avoid redundant computations
/// when calculating multiple metrics on the same structure pair.
/// 
/// Memory usage: ~4MB for 1000-residue protein
/// Speedup: 2-3x faster when calculating all metrics together
pub struct PrecomputedDistances {
    /// Squared distances between corresponding atoms: dist_sq(ref[i], model[i])
    /// Used by: TM-score, GDT, RMSD
    pub pairwise_dist: Vec<f32>,
    /// Number of atoms
    pub n: usize,
}

impl PrecomputedDistances {
    /// Pre-calculate all distances between two structures
    /// 
    /// # Arguments
    /// * `reference_coords` - Reference structure coordinates
    /// * `coords` - Model structure coordinates
    /// 
    /// # Returns
    /// PrecomputedDistances struct containing all distance calculations
    pub fn new(reference_coords: &[[f32; 3]], coords: &[[f32; 3]]) -> Self {
        let n = reference_coords.len();
        // Currently only supports equal-length structures. If empty or different lengths, return empty.
        if n == 0 || n != coords.len() {
            return Self {
                pairwise_dist: Vec::new(),
                n: 0,
            };
        }

        // Pre-calculate pairwise squared distances (for TM/GDT/RMSD)
        let mut pairwise_dist: Vec<f32> = Vec::with_capacity(n * n);

        for c in coords {
            for r in reference_coords {
                pairwise_dist.push(dist(*r, *c)); // This is after sqrt
            }
        }
        // Pairwise distances can be retrieved with pairwise_dist[i * n + j] 
        // instead of distances[i][j]
        
        Self {
            pairwise_dist,
            n,
        }
    }

    #[inline(always)]
    pub fn get_distance(&self, i: usize, j: usize) -> f32 {
        self.pairwise_dist[i * self.n + j]
    }

    
}

/// Calculate squared Euclidean distance between two 3D points.
/// This function uses f64 for intermediate calculations to improve precision.
#[inline(always)]
fn dist_sq_as_f64(a: [f32; 3], b: [f32; 3]) -> f64 {
    let dx = a[0] as f64 - b[0] as f64;
    let dy = a[1] as f64 - b[1] as f64;
    let dz = a[2] as f64 - b[2] as f64;
    (dx * dx + dy * dy + dz * dz) as f64
}

/// Calculate Euclidean distance between two 3D points
#[inline(always)]
fn dist(a: [f32; 3], b: [f32; 3]) -> f32 {
    dist_sq_as_f64(a, b).sqrt() as f32
}

/// Distance tolerance (Å) for the Distance Matrix Score (DMS).
pub const D_TOL: f32 = 2.0;

/// TM-score normalization function for final evaluation
/// d0(L) = 1.24 * (L - 15)^(1/3) - 1.8 for L > 19
/// d0(L) = 0.5 for L <= 21
#[inline]
fn d0_scale(length: usize) -> f32 {
    if length > 21 {
        1.24 * ((length as f32 - 15.0).powf(1.0 / 3.0)) - 1.8
    } else {
        0.5
    }
}

/// Fast TM-score using precomputed distances
/// 
/// # Arguments
/// * `distances` - Precomputed distance data
/// * `d0` - Optional normalization parameter
/// 
/// # Returns
/// TM-score in range [0, 1]
pub fn tm_score(distances: &PrecomputedDistances, d0: Option<f32>) -> f32 {
    if distances.n == 0 {
        return 0.0;
    }
    
    let d0 = d0.unwrap_or_else(|| d0_scale(distances.n));
    let d0_sq = (d0 * d0) as f64;
    
    let sum: f64 = (0..distances.n)
        .map(|i| {
            let d_sq = distances.get_distance(i, i) as f64;
            1.0 / (1.0 + d_sq / d0_sq)
        })
        .sum();
    
    
    (sum / distances.n as f64) as f32
}

/// Fast GDT using precomputed distances
fn gdt_generic(distances: &PrecomputedDistances, cutoffs: &[f64]) -> f32 {
    if distances.n == 0 || cutoffs.is_empty() {
        return 0.0;
    }
    
    let mut sum = 0.0_f64;
    
    for &cutoff in cutoffs {
        let cutoff_sq = cutoff * cutoff;
        let count = (0..distances.n)
            .filter(|&i| {
                let d_sq = distances.get_distance(i, i) as f64;
                d_sq <= cutoff_sq
            })
            .count();
        
        sum += (count as f64) / (distances.n as f64);
    }
    
    (sum / cutoffs.len() as f64) as f32
}

/// Fast GDT-TS using precomputed distances
pub fn gdt_ts(distances: &PrecomputedDistances) -> f32 {
    const CUTOFFS: [f64; 4] = [1.0, 2.0, 4.0, 8.0];
    gdt_generic(distances, &CUTOFFS)
}

/// Fast GDT-HA using precomputed distances
pub fn gdt_ha(distances: &PrecomputedDistances) -> f32 {
    const CUTOFFS: [f64; 4] = [0.5, 1.0, 2.0, 4.0];
    gdt_generic(distances, &CUTOFFS)
}

/// Fast GDT-strict using precomputed distances
// pub fn gdt_strict(distances: &PrecomputedDistances) -> f32 {
//     const CUTOFFS: [f64; 4] = [0.25, 0.5, 1.0, 2.0];
//     gdt_generic(distances, &CUTOFFS)
// }


/// Calculate Chamfer Distance between two point sets
/// 
/// Chamfer Distance is the mean of:
/// - Average nearest neighbor distance from coords to reference_coords
///
/// # Arguments
/// * `distance` - Precomputed distance data
///
/// # Returns
/// Chamfer distance (lower is better, 0 is perfect match)
pub fn chamfer_distance(distance: &PrecomputedDistances) -> f32 {
    if distance.n == 0 {
        return f32::INFINITY;
    }

    // Use f64 for accumulation to maintain precision
    // Average min distance from coords to reference
    let sum_coords_to_ref: f64 = (0..distance.n)
        .map(|i| {
            (0..distance.n)
                .map(|j| distance.get_distance(i, j) as f64)
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap()
        })
        .sum();
    
    (sum_coords_to_ref / distance.n as f64) as f32
}


/// Calculate Hausdorff Distance between two point sets
/// 
/// Hausdorff Distance is the maximum of:
/// - Maximum nearest neighbor distance from coords to reference_coords
///
/// # Arguments
/// * `distance` - Precomputed distance data
///
/// # Returns
/// Hausdorff distance (lower is better, 0 is perfect match)
pub fn hausdorff_distance(distance: &PrecomputedDistances) -> f32 {
    if distance.n == 0 {
        return f32::INFINITY;
    }

    // Max min distance from coords to reference
    let max_coords_to_ref = (0..distance.n)
        .map(|i| {
            (0..distance.n)
                .map(|j| distance.get_distance(i, j))
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap()
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();

    max_coords_to_ref
}


/// Calculate RMSD (Root Mean Square Deviation) between aligned structures
/// 
/// # Arguments
/// * `distance` - Precomputed distance data
///
/// # Returns
/// RMSD value in Ångströms
/// 
pub fn rmsd(distances: &PrecomputedDistances) -> f32 {
    if distances.n == 0 {
        return 0.0;
    }
    
    // Use f64 for accumulation to maintain precision
    let sum_sq: f64 = (0..distances.n)
        .map(|i| (distances.get_distance(i, i) as f64).powi(2))
        .sum();

    ((sum_sq / distances.n as f64).sqrt()) as f32
}

/// Distance Matrix Score (DMS) on Cα coordinates.
///
/// Returns 1.0 for N < 2, treating underspecified motifs as a neutral/perfect
/// agreement case because no pairwise distance disagreement can be measured.
pub fn distance_matrix_score(reference_coords: &[[f32; 3]], coords: &[[f32; 3]]) -> f32 {
    let n = reference_coords.len().min(coords.len());
    if n < 2 {
        return 1.0;
    }
    let total_pairs = n * (n - 1) / 2;
    let mut sum = 0.0_f32;
    for i in 0..n {
        for j in (i + 1)..n {
            let d_ref = dist(reference_coords[i], reference_coords[j]);
            let d_model = dist(coords[i], coords[j]);
            let diff = (d_ref - d_model).abs();
            sum += (1.0 - diff / D_TOL).max(0.0);
        }
    }
    sum / total_pairs as f32
}

#[inline]
fn pseudo_bond_angle(a: [f32; 3], b: [f32; 3], c: [f32; 3]) -> f32 {
    let v1 = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let v2 = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    let len1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
    let len2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
    if len1 < 1e-6 || len2 < 1e-6 {
        return 0.0;
    }
    (dot / (len1 * len2)).clamp(-1.0, 1.0).acos()
}

/// Pseudo-Bond Angle Score (PAS) on Cα coordinates.
///
/// Returns 1.0 for N < 3, treating underspecified motifs as a neutral/perfect
/// agreement case because no pseudo-bond angle disagreement can be measured.
pub fn pseudo_bond_angle_score(reference_coords: &[[f32; 3]], coords: &[[f32; 3]]) -> f32 {
    let n = reference_coords.len().min(coords.len());
    if n < 3 {
        return 1.0;
    }
    let count = n - 2;
    let sum: f32 = (0..count)
        .map(|k| {
            let theta_ref = pseudo_bond_angle(reference_coords[k], reference_coords[k + 1], reference_coords[k + 2]);
            let theta_model = pseudo_bond_angle(coords[k], coords[k + 1], coords[k + 2]);
            let delta = theta_ref - theta_model;
            (1.0 + delta.cos()) / 2.0
        })
        .sum();
    sum / count as f32
}

#[inline]
fn cb_direction(ca: [f32; 3], cb: [f32; 3]) -> Option<[f32; 3]> {
    let dx = cb[0] - ca[0];
    let dy = cb[1] - ca[1];
    let dz = cb[2] - ca[2];
    let len = (dx * dx + dy * dy + dz * dz).sqrt();
    if len < 1e-6 {
        None
    } else {
        Some([dx / len, dy / len, dz / len])
    }
}

/// Side-chain Orientation Score (SOS) using Cα/Cβ vectors.
///
/// Returns 1.0 when no valid residues exist for orientation comparison.
pub fn side_chain_orientation_score(
    reference_ca: &[[f32; 3]],
    coords_ca: &[[f32; 3]],
    reference_cb: &[[f32; 3]],
    coords_cb: &[[f32; 3]],
) -> f32 {
    let n = reference_ca
        .len()
        .min(coords_ca.len())
        .min(reference_cb.len())
        .min(coords_cb.len());
    if n == 0 {
        return 1.0;
    }
    let mut sum = 0.0_f32;
    let mut count = 0usize;
    for i in 0..n {
        if let (Some(d_ref), Some(d_model)) = (
            cb_direction(reference_ca[i], reference_cb[i]),
            cb_direction(coords_ca[i], coords_cb[i]),
        ) {
            let cos_sim = d_ref[0] * d_model[0] + d_ref[1] * d_model[1] + d_ref[2] * d_model[2];
            sum += (1.0 + cos_sim.clamp(-1.0, 1.0)) / 2.0;
            count += 1;
        }
    }
    if count == 0 {
        1.0
    } else {
        sum / count as f32
    }
}


/// Structure similarity metrics calculator
#[derive(Debug, Clone, Default, PartialEq, Copy)]
pub struct StructureSimilarityMetrics {
    pub tm_score: f32,
    pub gdt_ts: f32,
    pub gdt_ha: f32,
    pub dms: f32,
    pub pas: f32,
    pub sos: f32,
    pub chamfer_distance: f32,
    pub hausdorff_distance: f32,
}

impl StructureSimilarityMetrics {
    
    /// Calculate all metrics efficiently using precomputed distances
    /// 
    /// This is 2-3x faster than calculate_all() because distances are computed only once.
    /// Recommended for batch processing or when calculating multiple metrics.
    /// 
    /// # Arguments
    /// * `reference_coords` - Reference structure coordinates
    /// * `coords` - Model structure coordinates (should be pre-aligned)
    ///
    /// # Returns
    /// StructureMetrics containing all calculated metrics
    pub fn new() -> Self {      
        Self {
            tm_score: 0.0,
            gdt_ts: 0.0,
            gdt_ha: 0.0,
            dms: 0.0,
            pas: 0.0,
            sos: 0.0,
            chamfer_distance: 0.0,
            hausdorff_distance: 0.0,
        }
    }

    pub fn calculate_tm_score(&self, precomputed: &PrecomputedDistances) -> f32 {
        tm_score(precomputed, None)
    }

    pub fn calculate_gdt_ts(&self, precomputed: &PrecomputedDistances) -> f32 {
        gdt_ts(precomputed)
    }

    pub fn calculate_gdt_ha(&self, precomputed: &PrecomputedDistances) -> f32 {
        gdt_ha(precomputed)
    }

    pub fn calculate_chamfer_distance(&self, precomputed: &PrecomputedDistances) -> f32 {
        chamfer_distance(precomputed)
    }

    pub fn calculate_hausdorff_distance(&self, precomputed: &PrecomputedDistances) -> f32 {
        hausdorff_distance(precomputed)
    }

    pub fn calculate_dms(&self, reference_ca: &[[f32; 3]], coords_ca: &[[f32; 3]]) -> f32 {
        distance_matrix_score(reference_ca, coords_ca)
    }

    pub fn calculate_pas(&self, reference_ca: &[[f32; 3]], coords_ca: &[[f32; 3]]) -> f32 {
        pseudo_bond_angle_score(reference_ca, coords_ca)
    }

    pub fn calculate_sos(
        &self,
        reference_ca: &[[f32; 3]],
        coords_ca: &[[f32; 3]],
        reference_cb: &[[f32; 3]],
        coords_cb: &[[f32; 3]],
    ) -> f32 {
        side_chain_orientation_score(reference_ca, coords_ca, reference_cb, coords_cb)
    }

    pub fn calculate_all(&mut self, precomputed: &PrecomputedDistances) {
        self.tm_score = self.calculate_tm_score(precomputed);
        self.gdt_ts = self.calculate_gdt_ts(precomputed);
        self.gdt_ha = self.calculate_gdt_ha(precomputed);
        self.chamfer_distance = self.calculate_chamfer_distance(precomputed);
        self.hausdorff_distance = self.calculate_hausdorff_distance(precomputed);
    }

    /// Calculate DMS/PAS/SOS from interleaved [CA, CB, CA, CB, ...] coordinates.
    pub fn calculate_dms_pas_sos_from_interleaved_ca_cb(
        &mut self,
        reference_coords: &[[f32; 3]],
        coords: &[[f32; 3]],
    ) {
        if reference_coords.is_empty()
            || reference_coords.len() != coords.len()
            || reference_coords.len() % 2 != 0
        {
            self.dms = 0.0;
            self.pas = 0.0;
            self.sos = 0.0;
            return;
        }

        let mut reference_ca = Vec::with_capacity(reference_coords.len() / 2);
        let mut reference_cb = Vec::with_capacity(reference_coords.len() / 2);
        let mut coords_ca = Vec::with_capacity(coords.len() / 2);
        let mut coords_cb = Vec::with_capacity(coords.len() / 2);

        for i in (0..reference_coords.len()).step_by(2) {
            reference_ca.push(reference_coords[i]);
            reference_cb.push(reference_coords[i + 1]);
            coords_ca.push(coords[i]);
            coords_cb.push(coords[i + 1]);
        }

        self.dms = self.calculate_dms(&reference_ca, &coords_ca);
        self.pas = self.calculate_pas(&reference_ca, &coords_ca);
        self.sos = self.calculate_sos(&reference_ca, &coords_ca, &reference_cb, &coords_cb);
    }
    
    /// Print metrics in a formatted way
    pub fn print_in_a_formatted_way(&self) {
        println!("Structure Similarity Metrics:");
        println!("  TM-score:           {:.4}", self.tm_score);
        println!("  GDT-TS:             {:.4}", self.gdt_ts);
        println!("  GDT-HA:             {:.4}", self.gdt_ha);
        println!("  DMS:                {:.4}", self.dms);
        println!("  PAS:                {:.4}", self.pas);
        println!("  SOS:                {:.4}", self.sos);
        println!("  Chamfer Distance:   {:.4} Å", self.chamfer_distance);
        println!("  Hausdorff Distance: {:.4} Å", self.hausdorff_distance);
    }
}

impl fmt::Display for StructureSimilarityMetrics {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Print all metrics in a tab-separated format with 4 decimal places
        write!(
            f,
            "{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
            self.tm_score,
            self.gdt_ts,
            self.gdt_ha,
            self.dms,
            self.pas,
            self.sos,
            self.chamfer_distance,
            self.hausdorff_distance
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::structure::{kabsch::KabschSuperimposer, lms_qcp::LmsQcpSuperimposer};

    use super::*;

    #[inline]
    fn rotate_x(v: [f32; 3], angle: f32) -> [f32; 3] {
        let (s, c) = angle.sin_cos();
        [v[0], c * v[1] - s * v[2], s * v[1] + c * v[2]]
    }

    #[inline]
    fn rotate_z(v: [f32; 3], angle: f32) -> [f32; 3] {
        let (s, c) = angle.sin_cos();
        [c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]]
    }

    fn mock_per_residue_rt_false_positive(
        reference_interleaved: &[[f32; 3]],
        base_angle: f32,
        base_shift: f32,
        cb_extra_angle: f32,
    ) -> Vec<[f32; 3]> {
        let mut transformed = Vec::with_capacity(reference_interleaved.len());
        for i in (0..reference_interleaved.len()).step_by(2) {
            let residue_idx = (i / 2) as f32 + 1.0;
            let ca = reference_interleaved[i];
            let cb = reference_interleaved[i + 1];
            let translation = [
                base_shift * residue_idx * 0.31,
                -base_shift * residue_idx * 0.17,
                base_shift * residue_idx * 0.11,
            ];

            let ca_shifted = [
                ca[0] + translation[0],
                ca[1] + translation[1],
                ca[2] + translation[2],
            ];

            let mut ca_cb = [cb[0] - ca[0], cb[1] - ca[1], cb[2] - ca[2]];
            let angle = base_angle * residue_idx;
            ca_cb = rotate_z(ca_cb, angle);
            ca_cb = rotate_x(ca_cb, angle * 0.5 + cb_extra_angle);

            let cb_shifted = [
                ca_shifted[0] + ca_cb[0],
                ca_shifted[1] + ca_cb[1],
                ca_shifted[2] + ca_cb[2],
            ];

            transformed.push(ca_shifted);
            transformed.push(cb_shifted);
        }
        transformed
    }

    fn dms_pas_sos_average(
        reference_interleaved: &[[f32; 3]],
        model_interleaved: &[[f32; 3]],
    ) -> (f32, f32, f32, f32) {
        let mut metrics = StructureSimilarityMetrics::new();
        metrics.calculate_dms_pas_sos_from_interleaved_ca_cb(reference_interleaved, model_interleaved);
        let avg = (metrics.dms + metrics.pas + metrics.sos) / 3.0;
        (metrics.dms, metrics.pas, metrics.sos, avg)
    }

    #[test]
    fn test_real_motif_and_mocked_false_positive_validity() {
        use crate::prelude::PDBReader;

        let query_reader = PDBReader::from_file("query/1G2F.pdb").unwrap();
        let query_structure = query_reader.read_structure().unwrap().to_compact();
        let reference_indices = vec![
            query_structure.get_index(&b'F', &207).unwrap(),
            query_structure.get_index(&b'F', &212).unwrap(),
            query_structure.get_index(&b'F', &225).unwrap(),
            query_structure.get_index(&b'F', &229).unwrap(),
        ];

        let reference_interleaved = vec![
            query_structure.get_ca(reference_indices[0]).unwrap().to_array(),
            query_structure.get_cb(reference_indices[0]).unwrap().to_array(),
            query_structure.get_ca(reference_indices[1]).unwrap().to_array(),
            query_structure.get_cb(reference_indices[1]).unwrap().to_array(),
            query_structure.get_ca(reference_indices[2]).unwrap().to_array(),
            query_structure.get_cb(reference_indices[2]).unwrap().to_array(),
            query_structure.get_ca(reference_indices[3]).unwrap().to_array(),
            query_structure.get_cb(reference_indices[3]).unwrap().to_array(),
        ];

        // Composition-matched motif labels (same multiset)
        let aa_exact = [b'F', b'F', b'F', b'F'];
        let aa_composition_matched_different = [b'F', b'F', b'F', b'F'];
        assert_eq!(aa_exact, aa_composition_matched_different);

        let exact = reference_interleaved.clone();
        // Mild per-residue perturbation: small local rotation + sub-Å residue shifts.
        let slight_deviation = mock_per_residue_rt_false_positive(&reference_interleaved, 0.06, 0.12, 0.0);
        // Strong perturbation while preserving composition labels:
        // large local rotations and larger per-residue shifts.
        let composition_matched_different = mock_per_residue_rt_false_positive(&reference_interleaved, 0.9, 1.2, 0.75);

        let (_, _, _, exact_avg) = dms_pas_sos_average(&reference_interleaved, &exact);
        let (_, _, _, slight_avg) = dms_pas_sos_average(&reference_interleaved, &slight_deviation);
        let (_, _, _, different_avg) =
            dms_pas_sos_average(&reference_interleaved, &composition_matched_different);

        assert!((exact_avg - 1.0).abs() < 1e-6);
        assert!(exact_avg > slight_avg);
        assert!(slight_avg > different_avg);
    }

    #[test]
    fn test_dms_pas_sos_identical() {
        let ca = vec![
            [0.0_f32, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [3.0, 1.0, 1.0],
        ];
        let cb = vec![
            [0.0_f32, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 2.0, 0.0],
            [3.0, 2.0, 1.0],
        ];

        assert!((distance_matrix_score(&ca, &ca) - 1.0).abs() < 1e-6);
        assert!((pseudo_bond_angle_score(&ca, &ca) - 1.0).abs() < 1e-6);
        assert!((side_chain_orientation_score(&ca, &ca, &cb, &cb) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_dms_pas_sos_ranges() {
        let ca_ref = vec![
            [0.0_f32, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let ca_model = vec![
            [0.0_f32, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [1.5, -1.0, 0.0],
            [3.5, 0.0, 0.5],
        ];
        let cb_ref = vec![
            [0.0_f32, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let cb_model = vec![
            [0.0_f32, -1.0, 0.0],
            [0.5, 2.0, 0.0],
            [1.5, -2.0, 0.0],
            [3.5, -1.0, 0.5],
        ];

        let dms = distance_matrix_score(&ca_ref, &ca_model);
        let pas = pseudo_bond_angle_score(&ca_ref, &ca_model);
        let sos = side_chain_orientation_score(&ca_ref, &ca_model, &cb_ref, &cb_model);

        assert!((0.0..=1.0).contains(&dms));
        assert!((0.0..=1.0).contains(&pas));
        assert!((0.0..=1.0).contains(&sos));
    }

    #[test]
    fn test_dms_pas_sos_interleaved_input() {
        let interleaved_ref = vec![
            [0.0_f32, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 2.0, 0.0],
        ];
        let interleaved_model = interleaved_ref.clone();

        let mut metrics = StructureSimilarityMetrics::new();
        metrics.calculate_dms_pas_sos_from_interleaved_ca_cb(&interleaved_ref, &interleaved_model);

        assert!((metrics.dms - 1.0).abs() < 1e-6);
        assert!((metrics.pas - 1.0).abs() < 1e-6);
        assert!((metrics.sos - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_dms_pas_sos_interleaved_invalid_input() {
        let mut metrics = StructureSimilarityMetrics::new();
        metrics.calculate_dms_pas_sos_from_interleaved_ca_cb(&[[0.0, 0.0, 0.0]], &[[0.0, 0.0, 0.0]]);

        assert_eq!(metrics.dms, 0.0);
        assert_eq!(metrics.pas, 0.0);
        assert_eq!(metrics.sos, 0.0);
    }

    #[test]
    fn test_dms_pas_sos_interleaved_odd_length_input() {
        let mut metrics = StructureSimilarityMetrics::new();
        metrics.calculate_dms_pas_sos_from_interleaved_ca_cb(
            &[
                [0.0_f32, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
            ],
            &[
                [0.0_f32, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
            ],
        );

        assert_eq!(metrics.dms, 0.0);
        assert_eq!(metrics.pas, 0.0);
        assert_eq!(metrics.sos, 0.0);
    }

    #[test]
    fn test_metrics_calculate_all_with_identical() {
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let precomputed = PrecomputedDistances::new(&coords, &coords);
        let mut metrics = StructureSimilarityMetrics::new();
        metrics.calculate_all(&precomputed);

        assert!((metrics.tm_score - 1.0).abs() < 1e-6);
        assert!((metrics.gdt_ts - 1.0).abs() < 1e-6);
        assert!((metrics.gdt_ha - 1.0).abs() < 1e-6);
        assert!(metrics.chamfer_distance.abs() < 1e-6);
        assert!(metrics.hausdorff_distance.abs() < 1e-6);
        
        metrics.print_in_a_formatted_way();
    }
    
    
    #[test]
    fn test_with_real_coordinates() {
        // Read coordinates from PDB files
        use crate::prelude::PDBReader;
        let query_reader = PDBReader::from_file("query/1G2F.pdb").unwrap();
        let query_zinc_structure = query_reader.read_structure().unwrap().to_compact();
        let target_reader = PDBReader::from_file("data/zinc/AF-P36508-F1-model_v6.pdb").unwrap();
        let target_zinc_structure = target_reader.read_structure().unwrap().to_compact();
        // Get reference coordinates: F207,F212,F225,F229
        let reference_indices = vec![
            query_zinc_structure.get_index(&b'F', &207).unwrap(),
            query_zinc_structure.get_index(&b'F', &212).unwrap(),
            query_zinc_structure.get_index(&b'F', &225).unwrap(),
            query_zinc_structure.get_index(&b'F', &229).unwrap(),
        ];
        println!("Reference indices: {:?}", reference_indices);
        let reference_coords = vec![
            query_zinc_structure.get_ca(reference_indices[0]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[0]).unwrap(),
            query_zinc_structure.get_ca(reference_indices[1]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[1]).unwrap(),
            query_zinc_structure.get_ca(reference_indices[2]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[2]).unwrap(),
            query_zinc_structure.get_ca(reference_indices[3]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[3]).unwrap(),
        ];
        println!("Reference coords: {:?}", reference_coords);
        
        // Get target coordinates: A257,A262,A275,A279
        let target_indices = vec![
            target_zinc_structure.get_index(&b'A', &257).unwrap(),
            target_zinc_structure.get_index(&b'A', &262).unwrap(),
            target_zinc_structure.get_index(&b'A', &275).unwrap(),
            target_zinc_structure.get_index(&b'A', &279).unwrap(),
        ];
        println!("Target indices: {:?}", target_indices);
        let target_coords = vec![
            target_zinc_structure.get_ca(target_indices[0]).unwrap(),
            target_zinc_structure.get_cb(target_indices[0]).unwrap(),
            target_zinc_structure.get_ca(target_indices[1]).unwrap(),
            target_zinc_structure.get_cb(target_indices[1]).unwrap(),
            target_zinc_structure.get_ca(target_indices[2]).unwrap(),
            target_zinc_structure.get_cb(target_indices[2]).unwrap(),
            target_zinc_structure.get_ca(target_indices[3]).unwrap(),
            target_zinc_structure.get_cb(target_indices[3]).unwrap(),
        ];
        println!("Target coords: {:?}", target_coords);
        
        let mut kabsch = KabschSuperimposer::new();
        // let mut kabsch = LmsQcpSuperimposer::new();

        kabsch.set_atoms(&reference_coords, &target_coords);
        kabsch.run();

        let precomputed = PrecomputedDistances::new(
            &kabsch.reference_coords.unwrap(), &kabsch.transformed_coords.unwrap()
        );
        let mut metrics = StructureSimilarityMetrics::new();
        metrics.calculate_all(&precomputed);

        metrics.print_in_a_formatted_way();
        
    }

    
    #[test]
    fn test_with_long_coordinates() {
        // Read coordinates from PDB files
        use crate::prelude::PDBReader;
        let query_reader = PDBReader::from_file("query/1G2F.pdb").unwrap();
        let query_zinc_structure = query_reader.read_structure().unwrap().to_compact();
        let target_reader = PDBReader::from_file("data/zinc/AF-P36508-F1-model_v6.pdb").unwrap();
        let target_zinc_structure = target_reader.read_structure().unwrap().to_compact();
        // Get reference coordinates: F205-214,F223-232
        let mut reference_indices = (207..213).map(|res_num| {
            query_zinc_structure.get_index(&b'F', &res_num).unwrap()
        }).collect::<Vec<usize>>();
        reference_indices.extend((225..230).map(|res_num| {
            query_zinc_structure.get_index(&b'F', &res_num).unwrap()
        }));
        println!("Reference indices: {:?}", reference_indices);
        let reference_coords = reference_indices.iter().flat_map(|&idx| {
            vec![
                query_zinc_structure.get_n(idx).unwrap(),
                query_zinc_structure.get_ca(idx).unwrap(),
                query_zinc_structure.get_cb(idx).unwrap(),
            ]
        }).collect::<Vec<_>>();
        println!("Reference coords: {:?}", reference_coords);
        
        // Get target coordinates: A255-260,A273-282
        let mut target_indices = (256..262).map(|res_num| {
            target_zinc_structure.get_index(&b'A', &res_num).unwrap()
        }).collect::<Vec<usize>>();
        target_indices.extend((275..280).map(|res_num| {
            target_zinc_structure.get_index(&b'A', &res_num).unwrap()
        }));
        println!("Target indices: {:?}", target_indices);
        let target_coords = target_indices.iter().flat_map(|&idx| {
            vec![
                target_zinc_structure.get_n(idx).unwrap(),
                target_zinc_structure.get_ca(idx).unwrap(),
                target_zinc_structure.get_cb(idx).unwrap(),
            ]
        }).collect::<Vec<_>>();
        println!("Target coords: {:?}", target_coords);
        
        // let mut kabsch = KabschSuperimposer::new();
        let mut kabsch = LmsQcpSuperimposer::new();

        kabsch.set_atoms(&reference_coords, &target_coords);
        kabsch.run();

        let precomputed = PrecomputedDistances::new(
            &kabsch.reference_coords.unwrap(), &kabsch.transformed_coords.unwrap()
        );
        let mut metrics = StructureSimilarityMetrics::new();
        metrics.calculate_all(&precomputed);

        metrics.print_in_a_formatted_way();
        
    }
    
    #[test]
    fn test_with_outlier_coordinates() {
        // Read coordinates from PDB files
        use crate::prelude::PDBReader;
        let query_reader = PDBReader::from_file("query/1G2F.pdb").unwrap();
        let query_zinc_structure = query_reader.read_structure().unwrap().to_compact();
        let target_reader = PDBReader::from_file("data/zinc/AF-P36508-F1-model_v6.pdb").unwrap();
        let target_zinc_structure = target_reader.read_structure().unwrap().to_compact();
        // Get reference coordinates: F205,F212,F225,F229
        let reference_indices = vec![
            query_zinc_structure.get_index(&b'F', &205).unwrap(), // Outlier
            query_zinc_structure.get_index(&b'F', &212).unwrap(),
            query_zinc_structure.get_index(&b'F', &225).unwrap(),
            query_zinc_structure.get_index(&b'F', &229).unwrap(),
        ];
        println!("Reference indices: {:?}", reference_indices);
        let reference_coords = vec![
            query_zinc_structure.get_ca(reference_indices[0]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[0]).unwrap(),
            query_zinc_structure.get_ca(reference_indices[1]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[1]).unwrap(),
            query_zinc_structure.get_ca(reference_indices[2]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[2]).unwrap(),
            query_zinc_structure.get_ca(reference_indices[3]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[3]).unwrap(),
        ];
        println!("Reference coords: {:?}", reference_coords);
        
        // Get target coordinates: A257,A262,A275,A279
        let target_indices = vec![
            target_zinc_structure.get_index(&b'A', &257).unwrap(),
            target_zinc_structure.get_index(&b'A', &262).unwrap(),
            target_zinc_structure.get_index(&b'A', &275).unwrap(),
            target_zinc_structure.get_index(&b'A', &279).unwrap(),
        ];
        println!("Target indices: {:?}", target_indices);
        let target_coords = vec![
            target_zinc_structure.get_ca(target_indices[0]).unwrap(),
            target_zinc_structure.get_cb(target_indices[0]).unwrap(),
            target_zinc_structure.get_ca(target_indices[1]).unwrap(),
            target_zinc_structure.get_cb(target_indices[1]).unwrap(),
            target_zinc_structure.get_ca(target_indices[2]).unwrap(),
            target_zinc_structure.get_cb(target_indices[2]).unwrap(),
            target_zinc_structure.get_ca(target_indices[3]).unwrap(),
            target_zinc_structure.get_cb(target_indices[3]).unwrap(),
        ];
        println!("Target coords: {:?}", target_coords);
        
        // let mut kabsch = KabschSuperimposer::new();
        let mut kabsch = LmsQcpSuperimposer::new();
        kabsch.set_atoms(&reference_coords, &target_coords);
        kabsch.run();

        let precomputed = PrecomputedDistances::new(
            &kabsch.reference_coords.unwrap(), &kabsch.transformed_coords.unwrap()
        );
        let mut metrics = StructureSimilarityMetrics::new();
        metrics.calculate_all(&precomputed);

        metrics.print_in_a_formatted_way();
        
    }

}
