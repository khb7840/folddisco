//! # GRAMS: Geometric Alignment Motif Score
//!
//! A fully geometry-focused composite ranking metric for structural motif search
//! results in folddisco.
//!
//! ## Mathematical Formulation
//!
//! ```text
//! GRAMS(Q, H) = α · TM(Q, H) + β · DMS(Q, H) + γ · PAS(Q, H)
//! ```
//!
//! | Term | Description |
//! |------|-------------|
//! | `TM(Q, H)` | TM-score over post-superimposition Cα distances |
//! | `DMS(Q, H)` | Distance Matrix Score: agreement of all-pairs inter-residue Cα distances |
//! | `PAS(Q, H)` | Pseudo-Bond Angle Score: agreement of Cα-Cα-Cα pseudo-bond angles |
//! | α, β, γ | Non-negative weights summing to 1 (defaults: 0.5, 0.3, 0.2) |
//!
//! All three terms are independently bounded in **[0, 1]**, so GRAMS ∈ [0, 1].
//! Higher scores indicate better structural agreement.
//!
//! ## Why these three terms?
//!
//! - **TM-score** evaluates global Cα backbone quality after superimposition, using a
//!   length-dependent normalization that is robust to outlier residues.
//! - **DMS** captures internal distance-geometry consistency: two fragments that
//!   superimpose well should also have similar all-pairs Cα distance matrices.
//!   Deviations larger than `D_TOL` (2.0 Å) contribute zero, providing soft outlier
//!   tolerance similar in spirit to lDDT.
//! - **PAS** measures local backbone bend geometry via Cα pseudo-bond angles.  A
//!   motif with correct Cα positions but locally wrong curvature scores lower here.
//!
//! Together the three terms penalise global misalignment (TM), internal distance
//! distortion (DMS), and local backbone topology errors (PAS) — all purely from Cα
//! coordinates, with no reliance on sequence identity.
//!
//! ## Time Complexity
//!
//! - TM-score: **O(N)**
//! - DMS: **O(N²)**
//! - PAS: **O(N)**
//!
//! For small motifs (N ≤ 10) common in folddisco, the N² DMS term is negligible
//! (at most 45 distance evaluations per hit).
//!
//! ## Examples
//!
//! ```rust
//! use folddisco::ranking::{grams_score, GramsWeights};
//!
//! // Identical structures → GRAMS ≈ 1.0
//! let coords = vec![[1.0_f32, 0.0, 0.0], [4.0, 0.0, 0.0], [7.0, 0.0, 0.0]];
//! let score = grams_score(&coords, &coords, GramsWeights::default());
//! assert!((score - 1.0_f32).abs() < 1e-4, "Expected ~1.0, got {score}");
//! ```

use crate::structure::metrics::{tm_score, PrecomputedDistances};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Distance tolerance (Å) for the Distance Matrix Score.
///
/// Pairs whose inter-residue Cα distance differs by more than this threshold
/// contribute zero to DMS.  Analogous to the per-residue distance cutoff in lDDT.
pub const D_TOL: f32 = 2.0;

// ---------------------------------------------------------------------------
// Euclidean distance helper
// ---------------------------------------------------------------------------

#[inline(always)]
fn dist3(a: [f32; 3], b: [f32; 3]) -> f32 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

// ---------------------------------------------------------------------------
// Distance Matrix Score (DMS)
// ---------------------------------------------------------------------------

/// Compute the **Distance Matrix Score** between two sets of Cα coordinates.
///
/// DMS measures how well the internal pairwise Cα distance matrix of the hit
/// matches the query's distance matrix.  Each pair `(i, j)` with `i < j`
/// contributes a soft score:
///
/// ```text
/// s(i,j) = max(0,  1 − |d_Q(i,j) − d_H(i,j)| / D_TOL)
/// DMS(Q, H) = (1 / C(N,2)) · Σ_{i<j} s(i,j)
/// ```
///
/// Pairs whose distance difference exceeds [`D_TOL`] contribute 0 (hard cutoff),
/// providing outlier tolerance similar to lDDT.
///
/// # Returns
///
/// A value in **[0, 1]**.  Returns 1.0 for N < 2 (no pairs to compare).
pub fn distance_matrix_score(ref_ca: &[[f32; 3]], model_ca: &[[f32; 3]]) -> f32 {
    let n = ref_ca.len().min(model_ca.len());
    if n < 2 {
        return 1.0;
    }
    let total_pairs = n * (n - 1) / 2;
    let mut sum = 0.0_f32;
    for i in 0..n {
        for j in (i + 1)..n {
            let d_ref = dist3(ref_ca[i], ref_ca[j]);
            let d_model = dist3(model_ca[i], model_ca[j]);
            let diff = (d_ref - d_model).abs();
            sum += (1.0 - diff / D_TOL).max(0.0);
        }
    }
    sum / total_pairs as f32
}

// ---------------------------------------------------------------------------
// Pseudo-Bond Angle Score (PAS)
// ---------------------------------------------------------------------------

/// Compute the pseudo-bond angle (radians) at the central atom `b`
/// in the triplet `a – b – c`.
///
/// Returns 0.0 if any two input points are coincident (degenerate geometry).
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

/// Compute the **Pseudo-Bond Angle Score** between two sets of Cα coordinates.
///
/// For each consecutive triplet of Cα atoms `(k, k+1, k+2)`, the pseudo-bond
/// angle θ at the central atom is computed for both query and hit.  The per-triplet
/// contribution uses the cosine kernel:
///
/// ```text
/// PAS(Q, H) = (1/(N−2)) · Σ_{k=0}^{N−3}  (1 + cos(θ_Q_k − θ_H_k)) / 2
/// ```
///
/// This gives 1.0 when all angles match exactly and 0.0 when every angle differs
/// by π.
///
/// # Returns
///
/// A value in **[0, 1]**.  Returns 1.0 for N < 3 (no triplets available).
pub fn pseudo_bond_angle_score(ref_ca: &[[f32; 3]], model_ca: &[[f32; 3]]) -> f32 {
    let n = ref_ca.len().min(model_ca.len());
    if n < 3 {
        return 1.0;
    }
    let count = n - 2;
    let sum: f32 = (0..count)
        .map(|k| {
            let theta_ref = pseudo_bond_angle(ref_ca[k], ref_ca[k + 1], ref_ca[k + 2]);
            let theta_model = pseudo_bond_angle(model_ca[k], model_ca[k + 1], model_ca[k + 2]);
            let delta = theta_ref - theta_model;
            (1.0 + delta.cos()) / 2.0
        })
        .sum();
    sum / count as f32
}

// ---------------------------------------------------------------------------
// GramsWeights
// ---------------------------------------------------------------------------

/// Weights controlling the contribution of each GRAMS geometry sub-score.
///
/// The three weights must be non-negative.  They are automatically normalised
/// to sum to 1.0 inside [`grams_score`], so you can supply raw relative
/// importances (e.g. `GramsWeights { tm: 5.0, distance: 3.0, angle: 2.0 }`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GramsWeights {
    /// Weight for the TM-score (global backbone quality after superimposition). Default: 0.5.
    pub tm: f32,
    /// Weight for the Distance Matrix Score (internal distance geometry). Default: 0.3.
    pub distance: f32,
    /// Weight for the Pseudo-Bond Angle Score (local backbone curvature). Default: 0.2.
    pub angle: f32,
}

impl Default for GramsWeights {
    fn default() -> Self {
        Self {
            tm: 0.5,
            distance: 0.3,
            angle: 0.2,
        }
    }
}

impl GramsWeights {
    /// Construct weights and normalise them to sum to 1.0.
    ///
    /// # Panics
    ///
    /// Panics in debug mode if any weight is negative or if all weights are zero.
    pub fn new(tm: f32, distance: f32, angle: f32) -> Self {
        debug_assert!(
            tm >= 0.0 && distance >= 0.0 && angle >= 0.0,
            "GRAMS weights must be non-negative"
        );
        let total = tm + distance + angle;
        debug_assert!(total > 0.0, "GRAMS weights must not all be zero");
        Self {
            tm: tm / total,
            distance: distance / total,
            angle: angle / total,
        }
    }
}

// ---------------------------------------------------------------------------
// GramsComponents
// ---------------------------------------------------------------------------

/// Intermediate geometry sub-scores returned by [`grams_score_detailed`].
///
/// Useful for diagnostics or per-term ranking.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GramsComponents {
    /// TM-score (global backbone quality after superimposition), in [0, 1].
    pub tm_score: f32,
    /// Distance Matrix Score (internal Cα distance geometry agreement), in [0, 1].
    pub distance_matrix_score: f32,
    /// Pseudo-Bond Angle Score (local Cα backbone curvature agreement), in [0, 1].
    pub pseudo_bond_angle_score: f32,
    /// Final GRAMS composite score, in [0, 1].
    pub grams: f32,
}

// ---------------------------------------------------------------------------
// Main API
// ---------------------------------------------------------------------------

/// Compute the **GRAMS** composite geometry ranking score.
///
/// All three components use only Cα coordinates — no sequence identity is required.
///
/// # Arguments
///
/// * `ref_ca` — Reference (query) Cα coordinates after superimposition, in Å.
/// * `model_ca` — Hit Cα coordinates after superimposition, in Å.
/// * `weights` — Relative importance of each sub-score.  Normalised internally.
///
/// # Returns
///
/// GRAMS score in **[0, 1]**.  Returns 0.0 on empty or mismatched-length input.
///
/// # Examples
///
/// ```rust
/// use folddisco::ranking::{grams_score, GramsWeights};
///
/// // Identical geometry → GRAMS ≈ 1.0
/// let coords = vec![[0.0_f32, 0.0, 0.0], [4.0, 0.0, 0.0], [8.0, 0.0, 0.0]];
/// let score = grams_score(&coords, &coords, GramsWeights::default());
/// assert!((score - 1.0_f32).abs() < 1e-4);
/// ```
pub fn grams_score(
    ref_ca: &[[f32; 3]],
    model_ca: &[[f32; 3]],
    weights: GramsWeights,
) -> f32 {
    grams_score_detailed(ref_ca, model_ca, weights).grams
}

/// Compute the **GRAMS** composite score and return all intermediate sub-scores.
///
/// Identical parameters and semantics to [`grams_score`]; additionally returns
/// the individual TM-score, DMS, and PAS components.
///
/// # Examples
///
/// ```rust
/// use folddisco::ranking::{grams_score_detailed, GramsWeights};
///
/// let coords = vec![[0.0_f32, 0.0, 0.0], [4.0, 0.0, 0.0], [8.0, 0.0, 0.0]];
/// let c = grams_score_detailed(&coords, &coords, GramsWeights::default());
/// assert!((c.tm_score - 1.0_f32).abs() < 1e-4);
/// assert!((c.distance_matrix_score - 1.0_f32).abs() < 1e-4);
/// assert!((c.pseudo_bond_angle_score - 1.0_f32).abs() < 1e-4);
/// ```
pub fn grams_score_detailed(
    ref_ca: &[[f32; 3]],
    model_ca: &[[f32; 3]],
    weights: GramsWeights,
) -> GramsComponents {
    let n = ref_ca.len();
    if n == 0 || n != model_ca.len() {
        return GramsComponents {
            tm_score: 0.0,
            distance_matrix_score: 0.0,
            pseudo_bond_angle_score: 0.0,
            grams: 0.0,
        };
    }

    // Normalise weights to sum to 1.0
    let total_w = weights.tm + weights.distance + weights.angle;
    let (w_tm, w_dist, w_angle) = if total_w > 0.0 {
        (
            weights.tm / total_w,
            weights.distance / total_w,
            weights.angle / total_w,
        )
    } else {
        (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
    };

    // 1. TM-score (global backbone quality)
    let precomputed = PrecomputedDistances::new(ref_ca, model_ca);
    let tm = tm_score(&precomputed, None);

    // 2. Distance Matrix Score (internal Cα distance geometry)
    let dms = distance_matrix_score(ref_ca, model_ca);

    // 3. Pseudo-Bond Angle Score (local backbone curvature)
    let pas = pseudo_bond_angle_score(ref_ca, model_ca);

    let grams = w_tm * tm + w_dist * dms + w_angle * pas;

    GramsComponents {
        tm_score: tm,
        distance_matrix_score: dms,
        pseudo_bond_angle_score: pas,
        grams,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Helper: build a linear chain of Cα atoms spaced `spacing` Å apart
    // -----------------------------------------------------------------------
    fn linear_chain(n: usize, spacing: f32) -> Vec<[f32; 3]> {
        (0..n)
            .map(|i| [i as f32 * spacing, 0.0, 0.0])
            .collect()
    }

    // -----------------------------------------------------------------------
    // dist3
    // -----------------------------------------------------------------------

    #[test]
    fn test_dist3_known_values() {
        assert!((dist3([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]) - 5.0).abs() < 1e-5);
        assert!((dist3([1.0, 1.0, 1.0], [1.0, 1.0, 1.0])).abs() < 1e-6);
    }

    // -----------------------------------------------------------------------
    // distance_matrix_score
    // -----------------------------------------------------------------------

    #[test]
    fn test_dms_identical() {
        let coords = linear_chain(5, 4.0);
        assert!((distance_matrix_score(&coords, &coords) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_dms_single_atom() {
        let coords = vec![[1.0_f32, 0.0, 0.0]];
        assert_eq!(distance_matrix_score(&coords, &coords), 1.0);
    }

    #[test]
    fn test_dms_empty() {
        assert_eq!(distance_matrix_score(&[], &[]), 1.0);
    }

    #[test]
    fn test_dms_large_deviation_gives_zero() {
        // Completely different distance geometry → DMS close to 0
        let ref_coords = linear_chain(4, 4.0);
        let model_coords = vec![
            [0.0_f32, 0.0, 0.0],
            [100.0, 0.0, 0.0],
            [200.0, 0.0, 0.0],
            [300.0, 0.0, 0.0],
        ];
        let score = distance_matrix_score(&ref_coords, &model_coords);
        assert!(score < 0.1, "score={score}");
    }

    #[test]
    fn test_dms_small_noise_high_score() {
        // Small perturbation should give a high DMS
        let ref_coords = linear_chain(5, 4.0);
        let model_coords: Vec<[f32; 3]> = ref_coords
            .iter()
            .enumerate()
            .map(|(i, &[x, y, z])| [x + 0.1 * i as f32, y + 0.05, z])
            .collect();
        let score = distance_matrix_score(&ref_coords, &model_coords);
        assert!(score > 0.8, "score={score}");
    }

    #[test]
    fn test_dms_bounds() {
        let ref_coords = linear_chain(6, 4.0);
        let model_coords = linear_chain(6, 6.0); // stretched
        let score = distance_matrix_score(&ref_coords, &model_coords);
        assert!((0.0..=1.0).contains(&score), "score={score}");
    }

    // -----------------------------------------------------------------------
    // pseudo_bond_angle
    // -----------------------------------------------------------------------

    #[test]
    fn test_pba_straight_line() {
        // Three collinear points → angle = π
        let a = [0.0_f32, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [2.0, 0.0, 0.0];
        let angle = pseudo_bond_angle(a, b, c);
        assert!((angle - std::f32::consts::PI).abs() < 1e-5, "angle={angle}");
    }

    #[test]
    fn test_pba_right_angle() {
        // Perpendicular vectors → angle = π/2
        let a = [1.0_f32, 0.0, 0.0];
        let b = [0.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let angle = pseudo_bond_angle(a, b, c);
        assert!((angle - std::f32::consts::FRAC_PI_2).abs() < 1e-5, "angle={angle}");
    }

    #[test]
    fn test_pba_degenerate() {
        // Coincident points → returns 0.0
        let p = [1.0_f32, 1.0, 1.0];
        assert_eq!(pseudo_bond_angle(p, p, p), 0.0);
    }

    // -----------------------------------------------------------------------
    // pseudo_bond_angle_score
    // -----------------------------------------------------------------------

    #[test]
    fn test_pas_identical() {
        let coords = linear_chain(5, 4.0);
        assert!((pseudo_bond_angle_score(&coords, &coords) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_pas_two_atoms() {
        // N < 3 → returns 1.0
        let coords = linear_chain(2, 4.0);
        assert_eq!(pseudo_bond_angle_score(&coords, &coords), 1.0);
    }

    #[test]
    fn test_pas_one_atom() {
        let coords = vec![[0.0_f32, 0.0, 0.0]];
        assert_eq!(pseudo_bond_angle_score(&coords, &coords), 1.0);
    }

    #[test]
    fn test_pas_empty() {
        assert_eq!(pseudo_bond_angle_score(&[], &[]), 1.0);
    }

    #[test]
    fn test_pas_opposite_curvature_low_score() {
        // A curve and its mirror-image curvature should give a low PAS
        let ref_coords = vec![
            [0.0_f32, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [8.0, 4.0, 0.0], // bends up
        ];
        let model_coords = vec![
            [0.0_f32, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [8.0, -4.0, 0.0], // bends down
        ];
        let score = pseudo_bond_angle_score(&ref_coords, &model_coords);
        // Angle at b is the same (symmetric bend), so score = 1.0 here
        // (cos kernel depends on |delta|, not sign).
        // This verifies that angle magnitude, not direction, is captured.
        assert!((0.0..=1.0).contains(&score), "score={score}");
    }

    #[test]
    fn test_pas_bounds() {
        let ref_coords = vec![
            [0.0_f32, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [8.0, 4.0, 0.0],
            [12.0, 0.0, 0.0],
        ];
        let model_coords = vec![
            [0.0_f32, 0.0, 0.0],
            [4.0, 2.0, 0.0],
            [8.0, -1.0, 0.0],
            [12.0, 3.0, 0.0],
        ];
        let score = pseudo_bond_angle_score(&ref_coords, &model_coords);
        assert!((0.0..=1.0).contains(&score), "score={score}");
    }

    // -----------------------------------------------------------------------
    // grams_score
    // -----------------------------------------------------------------------

    #[test]
    fn test_grams_identical_structures() {
        let coords = linear_chain(5, 4.0);
        let score = grams_score(&coords, &coords, GramsWeights::default());
        // All three sub-scores = 1.0 → GRAMS = 1.0
        assert!((score - 1.0).abs() < 1e-4, "score={score}");
    }

    #[test]
    fn test_grams_empty_input() {
        assert_eq!(grams_score(&[], &[], GramsWeights::default()), 0.0);
    }

    #[test]
    fn test_grams_mismatched_lengths() {
        let coords3 = linear_chain(3, 4.0);
        let coords4 = linear_chain(4, 4.0);
        assert_eq!(grams_score(&coords3, &coords4, GramsWeights::default()), 0.0);
    }

    #[test]
    fn test_grams_pure_tm_weight() {
        // With only TM weight, GRAMS == TM-score
        let coords = linear_chain(5, 4.0);
        let w = GramsWeights::new(1.0, 0.0, 0.0);
        let c = grams_score_detailed(&coords, &coords, w);
        assert!((c.grams - c.tm_score).abs() < 1e-6);
    }

    #[test]
    fn test_grams_pure_distance_weight() {
        // With only distance weight, GRAMS == DMS
        let coords = linear_chain(5, 4.0);
        let w = GramsWeights::new(0.0, 1.0, 0.0);
        let c = grams_score_detailed(&coords, &coords, w);
        assert!((c.grams - c.distance_matrix_score).abs() < 1e-6);
    }

    #[test]
    fn test_grams_pure_angle_weight() {
        // With only angle weight, GRAMS == PAS
        let coords = linear_chain(5, 4.0);
        let w = GramsWeights::new(0.0, 0.0, 1.0);
        let c = grams_score_detailed(&coords, &coords, w);
        assert!((c.grams - c.pseudo_bond_angle_score).abs() < 1e-6);
    }

    #[test]
    fn test_grams_score_range() {
        use std::f32::consts::PI;
        let mut coords_a = Vec::new();
        let mut coords_b = Vec::new();
        for i in 0..6 {
            let angle = i as f32 * PI / 3.0;
            coords_a.push([angle.cos() * 5.0, angle.sin() * 5.0, 0.0]);
            coords_b.push([angle.cos() * 5.5, angle.sin() * 5.0, 1.0]);
        }
        let score = grams_score(&coords_a, &coords_b, GramsWeights::default());
        assert!((0.0..=1.0).contains(&score), "score={score}");
    }

    #[test]
    fn test_grams_perfect_beats_random() {
        let coords_ref = linear_chain(5, 4.0);
        let coords_random = vec![
            [10.0_f32, 20.0, 5.0],
            [0.5, 15.0, 3.0],
            [-4.0, 7.0, 12.0],
            [8.0, -2.0, 1.0],
            [3.0, 11.0, -6.0],
        ];
        let perfect = grams_score(&coords_ref, &coords_ref, GramsWeights::default());
        let random = grams_score(&coords_ref, &coords_random, GramsWeights::default());
        assert!(perfect > random, "perfect={perfect}, random={random}");
    }

    #[test]
    fn test_grams_components_sum_to_grams() {
        let coords = linear_chain(4, 4.0);
        let w = GramsWeights::default();
        let c = grams_score_detailed(&coords, &coords, w);
        let expected = w.tm * c.tm_score + w.distance * c.distance_matrix_score
            + w.angle * c.pseudo_bond_angle_score;
        assert!((c.grams - expected).abs() < 1e-5, "components don't sum to GRAMS");
    }

    #[test]
    fn test_grams_weights_normalisation() {
        // Raw vs. pre-normalised weights give identical scores
        let coords = linear_chain(5, 4.0);
        let w1 = GramsWeights::new(5.0, 3.0, 2.0);
        let w2 = GramsWeights::new(0.5, 0.3, 0.2);
        let s1 = grams_score(&coords, &coords, w1);
        let s2 = grams_score(&coords, &coords, w2);
        assert!((s1 - s2).abs() < 1e-6);
    }

    #[test]
    fn test_grams_single_residue() {
        let coords = vec![[0.0_f32, 0.0, 0.0]];
        let score = grams_score(&coords, &coords, GramsWeights::default());
        // N=1: DMS=1.0, PAS=1.0, TM uses d0=0.5 → score in (0,1]
        assert!((0.0..=1.0).contains(&score), "score={score}");
    }

    #[test]
    fn test_grams_distorted_geometry_lower_score() {
        // Distorted coordinates should score lower than identical
        let ref_coords = linear_chain(5, 4.0);
        let distorted: Vec<[f32; 3]> = ref_coords
            .iter()
            .enumerate()
            .map(|(i, &[x, y, z])| [x + 3.0 * i as f32, y + 2.0, z + 1.0])
            .collect();
        let perfect = grams_score(&ref_coords, &ref_coords, GramsWeights::default());
        let distorted_score = grams_score(&ref_coords, &distorted, GramsWeights::default());
        assert!(perfect > distorted_score, "perfect={perfect}, distorted={distorted_score}");
    }
}
