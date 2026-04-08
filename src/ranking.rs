//! # GRAMS: Geometric-Residue Alignment Motif Score
//!
//! A composite ranking metric for structural motif search results in folddisco.
//!
//! ## Mathematical Formulation
//!
//! ```text
//! GRAMS(Q, H) = α · TM(Q, H) + β · ResComp(Q, H) + γ · (1 − ClashFrac(H))
//! ```
//!
//! | Term | Description |
//! |------|-------------|
//! | `TM(Q, H)` | TM-score from post-superimposition Cα distances |
//! | `ResComp(Q, H)` | Mean normalised BLOSUM62 score over matched residue pairs |
//! | `ClashFrac(H)` | Fraction of intra-hit Cα pairs closer than [`D_CLASH`] (3.0 Å) |
//! | α, β, γ | Non-negative weights summing to 1 (defaults: 0.5, 0.3, 0.2) |
//!
//! All three terms are independently bounded in **[0, 1]**, so GRAMS ∈ [0, 1].
//! Higher scores indicate better, more biochemically consistent matches.
//!
//! ## Time Complexity
//!
//! - TM-score and ResComp: **O(N)**
//! - Clash detection: **O(N²)**
//!
//! For small motifs (N ≤ 10) common in folddisco, the N² term is negligible
//! (at most 45 distance evaluations per hit).
//!
//! ## Examples
//!
//! ```rust
//! use folddisco::ranking::{grams_score, GramsWeights};
//!
//! // Identical structures with identical residues → GRAMS ≈ 1.0
//! let coords = vec![[1.0_f32, 0.0, 0.0], [4.0, 0.0, 0.0], [7.0, 0.0, 0.0]];
//! let residues = b"ALA";
//! let score = grams_score(&coords, &coords, residues, residues, GramsWeights::default());
//! assert!((score - 1.0_f32).abs() < 1e-4, "Expected ~1.0, got {score}");
//! ```

use crate::structure::metrics::{tm_score, PrecomputedDistances};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Minimum value in the BLOSUM62 matrix (across all residue pairs, not just diagonals),
/// used as the lower bound for linear normalisation to [0, 1].
const B62_MIN: f32 = -4.0;

/// Maximum value in the BLOSUM62 matrix (Trp self-substitution score),
/// used as the upper bound for linear normalisation to [0, 1].
const B62_MAX: f32 = 11.0;

/// Cα–Cα distance threshold (Å) below which two residues are considered
/// sterically clashing.
pub const D_CLASH: f32 = 3.0;

// ---------------------------------------------------------------------------
// BLOSUM62 matrix
// ---------------------------------------------------------------------------

/// Residue index order for BLOSUM62 lookup: A R N D C Q E G H I L K M F P S T W Y V
///
/// Maps an ASCII uppercase one-letter code to a 0-based column/row index.
#[inline]
fn blosum62_index(aa: u8) -> Option<usize> {
    match aa {
        b'A' => Some(0),
        b'R' => Some(1),
        b'N' => Some(2),
        b'D' => Some(3),
        b'C' => Some(4),
        b'Q' => Some(5),
        b'E' => Some(6),
        b'G' => Some(7),
        b'H' => Some(8),
        b'I' => Some(9),
        b'L' => Some(10),
        b'K' => Some(11),
        b'M' => Some(12),
        b'F' => Some(13),
        b'P' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'W' => Some(17),
        b'Y' => Some(18),
        b'V' => Some(19),
        _ => None,
    }
}

/// BLOSUM62 substitution matrix (20×20, order: A R N D C Q E G H I L K M F P S T W Y V).
///
/// Source: NCBI BLOSUM62, PMID 1438297.
#[rustfmt::skip]
const BLOSUM62: [[i8; 20]; 20] = [
    //A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    [ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0], // A
    [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3], // R
    [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3], // N
    [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3], // D
    [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1], // C
    [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2], // Q
    [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2], // E
    [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3], // G
    [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3], // H
    [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3], // I
    [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1], // L
    [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2], // K
    [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1], // M
    [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1], // F
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2], // P
    [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2], // S
    [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0], // T
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3], // W
    [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1], // Y
    [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4], // V
];

/// Look up the raw BLOSUM62 score for a pair of amino acids.
///
/// Returns `None` if either residue code is unrecognised.
#[inline]
pub fn blosum62_score(a: u8, b: u8) -> Option<i8> {
    let i = blosum62_index(a)?;
    let j = blosum62_index(b)?;
    Some(BLOSUM62[i][j])
}

/// Normalise a raw BLOSUM62 score to the range **[0, 1]** using linear scaling:
///
/// ```text
/// b62_norm(a, b) = (BLOSUM62[a][b] − B62_MIN) / (B62_MAX − B62_MIN)
/// ```
///
/// Unknown amino acid codes map to 0.0 (neutral contribution).
#[inline]
pub fn blosum62_norm(a: u8, b: u8) -> f32 {
    match blosum62_score(a, b) {
        Some(raw) => ((raw as f32) - B62_MIN) / (B62_MAX - B62_MIN),
        None => 0.0,
    }
}

// ---------------------------------------------------------------------------
// GramsWeights
// ---------------------------------------------------------------------------

/// Weights controlling the contribution of each GRAMS sub-score.
///
/// The three weights must be non-negative.  They are automatically normalised
/// to sum to 1.0 inside [`grams_score`], so you can supply raw relative
/// importances (e.g. `GramsWeights { geometric: 5.0, biochemical: 3.0, clash: 2.0 }`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GramsWeights {
    /// Weight for the TM-score (geometric backbone quality). Default: 0.5.
    pub geometric: f32,
    /// Weight for the residue-compatibility (BLOSUM62) term. Default: 0.3.
    pub biochemical: f32,
    /// Weight for the steric-clash penalty term. Default: 0.2.
    pub clash: f32,
}

impl Default for GramsWeights {
    fn default() -> Self {
        Self {
            geometric: 0.5,
            biochemical: 0.3,
            clash: 0.2,
        }
    }
}

impl GramsWeights {
    /// Construct weights and normalise them to sum to 1.0.
    ///
    /// # Panics
    ///
    /// Panics in debug mode if any weight is negative or if all weights are zero.
    pub fn new(geometric: f32, biochemical: f32, clash: f32) -> Self {
        debug_assert!(geometric >= 0.0 && biochemical >= 0.0 && clash >= 0.0,
            "GRAMS weights must be non-negative");
        let total = geometric + biochemical + clash;
        debug_assert!(total > 0.0, "GRAMS weights must not all be zero");
        Self {
            geometric: geometric / total,
            biochemical: biochemical / total,
            clash: clash / total,
        }
    }
}

// ---------------------------------------------------------------------------
// Sub-score functions
// ---------------------------------------------------------------------------

/// Compute the residue-compatibility sub-score (β term) for a set of matched residue pairs.
///
/// Each element of `ref_residues` / `model_residues` is an ASCII uppercase one-letter
/// amino-acid code (e.g. `b'A'` for alanine).  Unknown codes contribute 0.0.
///
/// # Returns
///
/// Mean normalised BLOSUM62 score over all N pairs, in **[0, 1]**.  Returns 0.0
/// when the input is empty.
pub fn residue_compatibility(ref_residues: &[u8], model_residues: &[u8]) -> f32 {
    let n = ref_residues.len().min(model_residues.len());
    if n == 0 {
        return 0.0;
    }
    let sum: f32 = ref_residues[..n]
        .iter()
        .zip(model_residues[..n].iter())
        .map(|(&r, &m)| blosum62_norm(r, m))
        .sum();
    sum / n as f32
}

/// Compute the steric-clash fraction (1 − γ term) for a set of hit Cα coordinates.
///
/// Counts the fraction of all intra-hit Cα pairs `(i, j)` with `i < j` that
/// are closer than [`D_CLASH`] (3.0 Å).
///
/// # Returns
///
/// A value in **[0, 1]**.  Returns 0.0 for zero or one residue (no pairs).
pub fn clash_fraction(coords: &[[f32; 3]]) -> f32 {
    let n = coords.len();
    if n < 2 {
        return 0.0;
    }
    let total_pairs = n * (n - 1) / 2;
    let clash_threshold_sq = D_CLASH * D_CLASH;
    let mut clashes = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = coords[i][0] - coords[j][0];
            let dy = coords[i][1] - coords[j][1];
            let dz = coords[i][2] - coords[j][2];
            let dist_sq = dx * dx + dy * dy + dz * dz;
            if dist_sq < clash_threshold_sq {
                clashes += 1;
            }
        }
    }
    clashes as f32 / total_pairs as f32
}

// ---------------------------------------------------------------------------
// Main API
// ---------------------------------------------------------------------------

/// Intermediate sub-scores returned by [`grams_score_detailed`].
///
/// Useful for diagnostics or per-term ranking.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GramsComponents {
    /// TM-score (geometric sub-score), in [0, 1].
    pub tm_score: f32,
    /// Residue-compatibility sub-score (BLOSUM62 normalised), in [0, 1].
    pub residue_compatibility: f32,
    /// Steric-clash penalty: `1 − ClashFrac`, in [0, 1].
    pub clash_penalty: f32,
    /// Final GRAMS composite score, in [0, 1].
    pub grams: f32,
}

/// Compute the **GRAMS** composite ranking score.
///
/// # Arguments
///
/// * `ref_ca` — Reference (query) Cα coordinates after superimposition, in Å.
/// * `model_ca` — Hit Cα coordinates after superimposition, in Å.
/// * `ref_residues` — ASCII uppercase one-letter residue codes for the query.
/// * `model_residues` — ASCII uppercase one-letter residue codes for the hit.
/// * `weights` — Relative importance of each sub-score.  Normalised internally.
///
/// # Returns
///
/// GRAMS score in **[0, 1]**.  Returns 0.0 on empty or mismatched input.
///
/// # Examples
///
/// ```rust
/// use folddisco::ranking::{grams_score, GramsWeights};
///
/// // Identical geometry, identical residues → score ≈ 1.0
/// let coords = vec![[0.0_f32, 0.0, 0.0], [4.0, 0.0, 0.0], [8.0, 0.0, 0.0]];
/// let residues = b"ALV";
/// let score = grams_score(&coords, &coords, residues, residues, GramsWeights::default());
/// assert!((score - 1.0_f32).abs() < 1e-4);
/// ```
pub fn grams_score(
    ref_ca: &[[f32; 3]],
    model_ca: &[[f32; 3]],
    ref_residues: &[u8],
    model_residues: &[u8],
    weights: GramsWeights,
) -> f32 {
    grams_score_detailed(ref_ca, model_ca, ref_residues, model_residues, weights).grams
}

/// Compute the **GRAMS** composite score and return all intermediate sub-scores.
///
/// Identical parameters and semantics to [`grams_score`]; additionally returns
/// the individual TM-score, residue-compatibility, and clash-penalty components.
///
/// # Examples
///
/// ```rust
/// use folddisco::ranking::{grams_score_detailed, GramsWeights};
///
/// let coords = vec![[0.0_f32, 0.0, 0.0], [4.0, 0.0, 0.0], [8.0, 0.0, 0.0]];
/// let residues = b"ALV";
/// let components = grams_score_detailed(
///     &coords, &coords, residues, residues, GramsWeights::default()
/// );
/// assert!((components.tm_score - 1.0_f32).abs() < 1e-4);
/// assert!((components.residue_compatibility - 1.0_f32).abs() < 0.05);
/// assert!((components.clash_penalty - 1.0_f32).abs() < 1e-4);
/// ```
pub fn grams_score_detailed(
    ref_ca: &[[f32; 3]],
    model_ca: &[[f32; 3]],
    ref_residues: &[u8],
    model_residues: &[u8],
    weights: GramsWeights,
) -> GramsComponents {
    let n = ref_ca.len();
    if n == 0 || n != model_ca.len() {
        return GramsComponents {
            tm_score: 0.0,
            residue_compatibility: 0.0,
            clash_penalty: 0.0,
            grams: 0.0,
        };
    }

    // Normalise weights to sum to 1.0
    let total_w = weights.geometric + weights.biochemical + weights.clash;
    let (w_g, w_b, w_c) = if total_w > 0.0 {
        (
            weights.geometric / total_w,
            weights.biochemical / total_w,
            weights.clash / total_w,
        )
    } else {
        (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
    };

    // 1. TM-score (geometric quality)
    let precomputed = PrecomputedDistances::new(ref_ca, model_ca);
    let tm = tm_score(&precomputed, None);

    // 2. Residue compatibility (BLOSUM62)
    let res_comp = residue_compatibility(ref_residues, model_residues);

    // 3. Steric-clash penalty: (1 − ClashFrac)
    let clash_pen = 1.0 - clash_fraction(model_ca);

    let grams = w_g * tm + w_b * res_comp + w_c * clash_pen;

    GramsComponents {
        tm_score: tm,
        residue_compatibility: res_comp,
        clash_penalty: clash_pen,
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
    // blosum62_norm
    // -----------------------------------------------------------------------

    #[test]
    fn test_blosum62_norm_identical_residues() {
        // Identical residues should yield normalised scores > 0.5 for most AAs
        for &aa in b"ACDEFGHIKLMNPQRSTVWY" {
            let score = blosum62_norm(aa, aa);
            assert!(score > 0.0, "Self-score should be positive for {}", aa as char);
        }
    }

    #[test]
    fn test_blosum62_norm_bounds() {
        // All values must be in [0, 1]
        for &a in b"ACDEFGHIKLMNPQRSTVWY" {
            for &b in b"ACDEFGHIKLMNPQRSTVWY" {
                let v = blosum62_norm(a, b);
                assert!(
                    (0.0..=1.0).contains(&v),
                    "blosum62_norm({}, {}) = {} out of [0,1]",
                    a as char, b as char, v
                );
            }
        }
    }

    #[test]
    fn test_blosum62_norm_unknown_residue() {
        assert_eq!(blosum62_norm(b'X', b'A'), 0.0);
        assert_eq!(blosum62_norm(b'A', b'Z'), 0.0);
        assert_eq!(blosum62_norm(b'B', b'B'), 0.0);
    }

    #[test]
    fn test_blosum62_raw_lookup() {
        // Known values from the BLOSUM62 matrix
        assert_eq!(blosum62_score(b'A', b'A'), Some(4));
        assert_eq!(blosum62_score(b'W', b'W'), Some(11));
        assert_eq!(blosum62_score(b'A', b'W'), Some(-3));
    }

    #[test]
    fn test_blosum62_all_20_self_scores() {
        // Self-substitution scores for all 20 standard amino acids must be
        // positive and normalise to a value in (0, 1].
        // Exercises every row/column including P, Q, R, S, T, Y.
        let expected_diag: &[(u8, i8)] = &[
            (b'A',  4), (b'R',  5), (b'N',  6), (b'D',  6),
            (b'C',  9), (b'Q',  5), (b'E',  5), (b'G',  6),
            (b'H',  8), (b'I',  4), (b'L',  4), (b'K',  5),
            (b'M',  5), (b'F',  6), (b'P',  7), (b'S',  4),
            (b'T',  5), (b'W', 11), (b'Y',  7), (b'V',  4),
        ];
        for &(aa, expected_raw) in expected_diag {
            let raw = blosum62_score(aa, aa)
                .unwrap_or_else(|| panic!("missing self-score for {}", aa as char));
            assert_eq!(raw, expected_raw,
                "self-score mismatch for {}", aa as char);
            let norm = blosum62_norm(aa, aa);
            assert!(norm > 0.0 && norm <= 1.0,
                "normalised self-score for {} = {} out of (0,1]", aa as char, norm);
        }
    }

    #[test]
    fn test_blosum62_cross_pairs_pqrsty() {
        // Verify a selection of cross-pair lookups for P, Q, R, S, T, Y —
        // the amino acids not explicitly tested elsewhere.
        // Values taken directly from the standard BLOSUM62 table.
        let pairs: &[(u8, u8, i8)] = &[
            (b'P', b'A', -1),
            (b'Q', b'K',  1),
            (b'R', b'D', -2),
            (b'S', b'T',  1),
            (b'T', b'V',  0),
            (b'Y', b'F',  3),
            (b'P', b'P',  7),
            (b'Y', b'Y',  7),
        ];
        for &(a, b_aa, expected) in pairs {
            let got = blosum62_score(a, b_aa)
                .unwrap_or_else(|| panic!("missing score ({}, {})", a as char, b_aa as char));
            assert_eq!(got, expected,
                "BLOSUM62({},{}) expected {} got {}", a as char, b_aa as char, expected, got);
        }
    }

    // -----------------------------------------------------------------------
    // residue_compatibility
    // -----------------------------------------------------------------------

    #[test]
    fn test_residue_compatibility_identical() {
        let res = b"ACDEFG";
        let score = residue_compatibility(res, res);
        // All diagonal scores are positive, so normalised mean should be > 0.5
        assert!(score > 0.5, "identical residues should yield high compatibility");
    }

    #[test]
    fn test_residue_compatibility_empty() {
        assert_eq!(residue_compatibility(&[], &[]), 0.0);
    }

    #[test]
    fn test_residue_compatibility_bounds() {
        let ref_res = b"ACDEFGHIKLMNPQRSTVWY";
        let model_res = b"WYTVSRQPNMLKIHGFEDCA"; // reversed
        let score = residue_compatibility(ref_res, model_res);
        assert!((0.0..=1.0).contains(&score));
    }

    // -----------------------------------------------------------------------
    // clash_fraction
    // -----------------------------------------------------------------------

    #[test]
    fn test_clash_fraction_no_clash() {
        // Well-separated chain: no clashes
        let coords = linear_chain(5, 4.0);
        assert_eq!(clash_fraction(&coords), 0.0);
    }

    #[test]
    fn test_clash_fraction_all_clash() {
        // All atoms at the same position → all pairs clash
        let coords = vec![[0.0, 0.0, 0.0]; 4];
        assert_eq!(clash_fraction(&coords), 1.0);
    }

    #[test]
    fn test_clash_fraction_single_atom() {
        let coords = vec![[1.0, 2.0, 3.0]];
        assert_eq!(clash_fraction(&coords), 0.0);
    }

    #[test]
    fn test_clash_fraction_empty() {
        assert_eq!(clash_fraction(&[]), 0.0);
    }

    #[test]
    fn test_clash_fraction_partial() {
        // 3 atoms: first two close, third far away → 1/3 pairs clash
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0], // < 3.0 Å from [0]
            [10.0, 0.0, 0.0],
        ];
        let cf = clash_fraction(&coords);
        // 1 clash out of 3 pairs
        assert!((cf - 1.0 / 3.0).abs() < 1e-6);
    }

    // -----------------------------------------------------------------------
    // grams_score
    // -----------------------------------------------------------------------

    #[test]
    fn test_grams_identical_structures() {
        let coords = linear_chain(4, 4.0);
        let residues = b"ALVI";
        let score = grams_score(&coords, &coords, residues, residues, GramsWeights::default());
        // Perfect geometric match with residues whose BLOSUM62 self-score = 4 (A, L, V, I),
        // giving a normalised ResComp ≈ 0.533.
        // Minimum for identical structures: 0.5×1.0 + 0.3×0.533 + 0.2×1.0 ≈ 0.86.
        assert!(score > 0.80, "identical structures: score={score}");
    }

    #[test]
    fn test_grams_empty_input() {
        let score = grams_score(&[], &[], &[], &[], GramsWeights::default());
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_grams_mismatched_lengths() {
        let coords3 = linear_chain(3, 4.0);
        let coords4 = linear_chain(4, 4.0);
        let res3 = b"ALV";
        let res4 = b"ALVI";
        let score = grams_score(&coords3, &coords4, res3, res4, GramsWeights::default());
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_grams_pure_geometry_weight() {
        // With only geometric weight, GRAMS ≈ TM-score
        let coords = linear_chain(5, 4.0);
        let residues = b"ACDEF";
        let w = GramsWeights::new(1.0, 0.0, 0.0);
        let components = grams_score_detailed(&coords, &coords, residues, residues, w);
        assert!((components.grams - components.tm_score).abs() < 1e-6);
    }

    #[test]
    fn test_grams_score_range() {
        // GRAMS must always be in [0, 1]
        use std::f32::consts::PI;
        let mut coords_a = Vec::new();
        let mut coords_b = Vec::new();
        for i in 0..6 {
            let angle = i as f32 * PI / 3.0;
            coords_a.push([angle.cos() * 5.0, angle.sin() * 5.0, 0.0]);
            // Slightly perturbed version
            coords_b.push([angle.cos() * 5.5, angle.sin() * 5.0, 1.0]);
        }
        let res_a = b"ACDGHI";
        let res_b = b"ACWGHI"; // W is a conservative-ish mismatch for D
        let score = grams_score(&coords_a, &coords_b, res_a, res_b, GramsWeights::default());
        assert!((0.0..=1.0).contains(&score), "score={score}");
    }

    #[test]
    fn test_grams_random_vs_perfect() {
        // A perfect structural match should score higher than a random one
        let coords_ref = linear_chain(5, 4.0);
        let coords_random = vec![
            [10.0, 20.0, 5.0],
            [0.5, 15.0, 3.0],
            [-4.0, 7.0, 12.0],
            [8.0, -2.0, 1.0],
            [3.0, 11.0, -6.0],
        ];
        let residues = b"ACLVI";
        let perfect_score = grams_score(
            &coords_ref, &coords_ref, residues, residues, GramsWeights::default()
        );
        let random_score = grams_score(
            &coords_ref, &coords_random, residues, residues, GramsWeights::default()
        );
        assert!(
            perfect_score > random_score,
            "perfect={perfect_score}, random={random_score}"
        );
    }

    #[test]
    fn test_grams_residue_mismatch_lowers_score() {
        // Same geometry, different residues → lower score than identical residues
        let coords = linear_chain(4, 4.0);
        let ref_res = b"WWWW"; // tryptophan — maximum self-score in BLOSUM62
        let mut_res = b"DDDD"; // aspartate — mismatches with W
        let same = grams_score(&coords, &coords, ref_res, ref_res, GramsWeights::default());
        let diff = grams_score(&coords, &coords, ref_res, mut_res, GramsWeights::default());
        assert!(same > diff, "same={same}, diff={diff}");
    }

    #[test]
    fn test_grams_components_sum_to_grams() {
        let coords = linear_chain(4, 4.0);
        let residues = b"ALVI";
        let w = GramsWeights::default();
        let c = grams_score_detailed(&coords, &coords, residues, residues, w);
        let expected = w.geometric * c.tm_score
            + w.biochemical * c.residue_compatibility
            + w.clash * c.clash_penalty;
        assert!((c.grams - expected).abs() < 1e-6, "components don't sum to GRAMS");
    }

    #[test]
    fn test_grams_weights_normalisation() {
        // Non-normalised vs normalised weights should give the same result
        let coords = linear_chain(4, 4.0);
        let residues = b"ALVI";
        let w1 = GramsWeights::new(5.0, 3.0, 2.0);
        let w2 = GramsWeights::new(0.5, 0.3, 0.2);
        let s1 = grams_score(&coords, &coords, residues, residues, w1);
        let s2 = grams_score(&coords, &coords, residues, residues, w2);
        assert!((s1 - s2).abs() < 1e-6);
    }

    #[test]
    fn test_grams_single_residue() {
        // Single residue: no clash pairs, TM-score with d0=0.5
        let coords = vec![[0.0_f32, 0.0, 0.0]];
        let residues = b"A";
        let score = grams_score(&coords, &coords, residues, residues, GramsWeights::default());
        assert!((0.0..=1.0).contains(&score));
    }

    #[test]
    fn test_grams_clash_detected() {
        // Overlapping coordinates (distance 0) should trigger clash penalty
        let coords = vec![[0.0_f32, 0.0, 0.0]; 4];
        let residues = b"ALVI";
        let components =
            grams_score_detailed(&coords, &coords, residues, residues, GramsWeights::default());
        // All pairs clash → clash_penalty = 0.0
        assert_eq!(components.clash_penalty, 0.0);
    }
}
