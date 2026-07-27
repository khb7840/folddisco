// File: expand.rs
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2026 Hyunbin Kim, All rights reserved
//
// Tolerance driven expansion of a query feature vector.
//
// Folddisco discretizes the geometry of a residue pair into bins and looks the
// resulting hash up in the inverted index. A target pair that is geometrically
// close to the query still gets a different hash whenever one of its features sits
// on the other side of a bin boundary, so the query has to be expanded into the
// neighbourhood of the observed geometry.
//
// The expansion here differs from a plain per-dimension offset in four ways:
//
// 1. Dimensions are combined. Two features straddling a boundary at the same time
//    is common, and a query that only ever moves one dimension misses it. The
//    number of dimensions allowed to deviate at once is the `radius`.
// 2. A tolerance wider than one bin is sub-stepped, so it reaches every bin in the
//    interval instead of only the two at its ends.
// 3. Distance tolerance can grow with the distance itself (`dist_ratio`), which is
//    what a non-rigid deformation does to a motif: a 16 A pair drifts much further
//    than a 5 A one.
// 4. Perturbed angles are pulled back into their domain (torsions wrap, bounded
//    angles reflect) instead of running off the end of their bit field or asking for
//    a `sin` sign no real structure can produce. Distances stay untouched so the
//    query keeps mirroring the encoding the index was built with.
//
// Variants are produced nearest-first: every single-dimension neighbour before any
// two-dimension neighbour, and closer offsets before farther ones. Callers can
// therefore stop early and keep the most useful part of the neighbourhood.

use crate::geometry::core::HashType;

/// Cap on the sub-steps generated for one dimension, so a very large threshold
/// cannot blow up the query on its own.
const MAX_SUBSTEPS: usize = 8;

/// How far a query hash is expanded around the observed geometry.
#[derive(Debug, Clone, PartialEq)]
pub struct ToleranceConfig {
    /// Distance offsets in Angstroms.
    pub dist_thresholds: Vec<f32>,
    /// Angle offsets in degrees, as given on the command line.
    pub angle_thresholds: Vec<f32>,
    /// Extra distance tolerance as a fraction of the pair distance. Models elastic
    /// (non-rigid) deformation, where longer pairs move further.
    pub dist_ratio: f32,
    /// Number of feature dimensions allowed to deviate from the observed bin at
    /// the same time. 1 reproduces the classic one-dimension-at-a-time expansion.
    pub radius: usize,
}

impl ToleranceConfig {
    pub fn new(
        dist_thresholds: Vec<f32>, angle_thresholds: Vec<f32>,
        dist_ratio: f32, radius: usize,
    ) -> Self {
        Self { dist_thresholds, angle_thresholds, dist_ratio, radius }
    }

    /// Tolerances used by `folddisco query` when nothing is given.
    pub fn default_query() -> Self {
        Self::new(vec![0.5], vec![5.0], 0.0, 1)
    }

    /// No expansion at all: only the observed geometry is queried.
    pub fn none() -> Self {
        Self::new(Vec::new(), Vec::new(), 0.0, 0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum DimKind {
    Distance,
    Angle,
}

#[derive(Debug, Clone)]
struct TolerantDim {
    /// Position of this dimension inside the feature vector
    index: usize,
    kind: DimKind,
    /// Widest tolerance for this dimension, in the unit the feature uses. Only the
    /// widest one is needed: sub-stepping it reaches a contiguous run of bins that
    /// already contains every bin a smaller tolerance could reach.
    tolerance: f32,
    /// Offsets for a fixed tolerance, nearest first. Distance dimensions rebuild
    /// theirs per residue pair when `dist_ratio` is in play.
    offsets: Vec<f32>,
}

/// Enumerates the tolerance neighbourhood of a feature vector for one hash type.
///
/// Holds its own scratch buffers, so expanding thousands of residue pairs does not
/// allocate after construction.
pub struct FeatureExpander {
    hash_type: HashType,
    dims: Vec<TolerantDim>,
    radius: usize,
    dist_ratio: f32,
    dist_bin_width: f32,
    /// Feature vector with offsets applied, before domain fix-up
    raw: Vec<f32>,
    /// `raw` pulled back into the valid domain; this is what callers see
    fixed: Vec<f32>,
    /// Per-dimension offsets for the pair currently being expanded
    active: Vec<Vec<f32>>,
}

impl FeatureExpander {
    pub fn new(
        hash_type: HashType, nbin_dist: usize, nbin_angle: usize, tolerance: &ToleranceConfig,
    ) -> Self {
        let dist_bin_width = hash_type.dist_bin_width(nbin_dist);
        let angle_bin_width = hash_type.angle_bin_width(nbin_angle);
        let dist_ratio = tolerance.dist_ratio.max(0.0);
        let mut dims: Vec<TolerantDim> = Vec::new();

        if let Some(dist_indices) = hash_type.dist_index() {
            let widest = widest_tolerance(&tolerance.dist_thresholds);
            if widest > 0.0 || dist_ratio > 0.0 {
                let offsets = substep_offsets(widest, dist_bin_width);
                for index in dist_indices {
                    dims.push(TolerantDim {
                        index, kind: DimKind::Distance, tolerance: widest, offsets: offsets.clone(),
                    });
                }
            }
        }

        if let Some(angle_indices) = hash_type.angle_index() {
            // Command line thresholds are degrees; every hash type except PDBMotif
            // stores angles in radians.
            let widest = widest_tolerance(&tolerance.angle_thresholds);
            let widest = if hash_type.angle_in_degrees() { widest } else { widest.to_radians() };
            if widest > 0.0 {
                let offsets = substep_offsets(widest, angle_bin_width);
                for index in angle_indices {
                    dims.push(TolerantDim {
                        index, kind: DimKind::Angle, tolerance: widest, offsets: offsets.clone(),
                    });
                }
            }
        }

        let active = vec![Vec::new(); dims.len()];
        Self {
            hash_type,
            dims,
            radius: tolerance.radius,
            dist_ratio,
            dist_bin_width,
            raw: vec![0.0; 9],
            fixed: vec![0.0; 9],
            active,
        }
    }

    /// True when this expander would never produce a neighbour.
    pub fn is_noop(&self) -> bool {
        self.radius == 0
            || self.dims.is_empty()
            || (self.dist_ratio == 0.0 && self.dims.iter().all(|d| d.offsets.is_empty()))
    }

    /// Number of feature dimensions that carry a tolerance.
    pub fn dim_count(&self) -> usize {
        self.dims.len()
    }

    /// Visit every neighbour of `feature`, nearest first. The observed feature
    /// itself is *not* visited: the caller owns that hash because it also carries
    /// the IDF of the edge.
    ///
    /// `visit` returns `false` to stop the expansion early, which keeps the closest
    /// neighbours and drops the rest.
    pub fn for_each_neighbor<F>(&mut self, feature: &[f32], mut visit: F)
    where
        F: FnMut(&[f32]) -> bool,
    {
        if self.is_noop() {
            return;
        }
        // Offsets for this pair. Distance dimensions pick up the proportional part,
        // which depends on the observed distance; angles reuse the fixed offsets
        // computed at construction.
        for i in 0..self.dims.len() {
            let elastic = match self.dims[i].kind {
                DimKind::Distance => self.dist_ratio * feature[self.dims[i].index].abs(),
                DimKind::Angle => 0.0,
            };
            self.active[i].clear();
            if elastic > 0.0 {
                let widened = substep_offsets(
                    self.dims[i].tolerance + elastic, self.dist_bin_width
                );
                self.active[i].extend_from_slice(&widened);
            } else {
                self.active[i].extend_from_slice(&self.dims[i].offsets);
            }
        }

        let radius = self.radius.min(self.dims.len());
        self.raw[..feature.len()].copy_from_slice(feature);
        for level in 1..=radius {
            if !self.visit_level(level, 0, feature, &mut visit) {
                return;
            }
        }
    }

    /// Choose `level` distinct dimensions from `dims[start..]` and offset each one.
    /// Returns false once `visit` asked to stop.
    fn visit_level<F>(&mut self, level: usize, start: usize, feature: &[f32], visit: &mut F) -> bool
    where
        F: FnMut(&[f32]) -> bool,
    {
        if level == 0 {
            let len = feature.len();
            self.fixed[..len].copy_from_slice(&self.raw[..len]);
            self.hash_type.sanitize_perturbed_feature(&mut self.fixed[..len]);
            return visit(&self.fixed[..len]);
        }
        // Leave room for the remaining levels
        let last = self.dims.len() - level;
        for dim in start..=last {
            let index = self.dims[dim].index;
            let original = self.raw[index];
            for offset_pos in 0..self.active[dim].len() {
                self.raw[index] = original + self.active[dim][offset_pos];
                if !self.visit_level(level - 1, dim + 1, feature, visit) {
                    self.raw[index] = original;
                    return false;
                }
            }
            self.raw[index] = original;
        }
        true
    }
}

/// Widest magnitude in a threshold list. Smaller thresholds are redundant: the
/// sub-stepped offsets of the widest one cover a contiguous run of bins that
/// already includes every bin a smaller threshold could reach.
fn widest_tolerance(thresholds: &[f32]) -> f32 {
    thresholds.iter().fold(0.0f32, |widest, t| widest.max(t.abs()))
}

/// Signed offsets covering `[-tolerance, +tolerance]`, nearest first.
///
/// Consecutive offsets are kept within one bin of each other so the covered bins
/// form an unbroken run; a single jump of 2.5 bins would silently skip the bins in
/// between and lose the targets that sit there.
fn substep_offsets(tolerance: f32, bin_width: f32) -> Vec<f32> {
    if !(tolerance > 0.0) {
        return Vec::new();
    }
    let steps = if bin_width.is_finite() && bin_width > 0.0 {
        ((tolerance / bin_width).ceil() as usize).clamp(1, MAX_SUBSTEPS)
    } else {
        1
    };
    let mut offsets = Vec::with_capacity(steps * 2);
    for step in 1..=steps {
        let magnitude = tolerance * step as f32 / steps as f32;
        offsets.push(-magnitude);
        offsets.push(magnitude);
    }
    offsets
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::convert::{reflect_into_range, wrap_to_pi};
    use crate::geometry::core::GeometricHash;
    use rustc_hash::FxHashSet as HashSet;

    const PI: f32 = std::f32::consts::PI;

    /// PDBTrRosetta feature: aa1, aa2, ca_dist, cb_dist, ca_cb_angle, theta1, theta2
    fn pdbtr_feature() -> Vec<f32> {
        let mut feature = vec![0.0; 9];
        feature[0] = 8.0;   // HIS
        feature[1] = 3.0;   // ASP
        feature[2] = 9.4;
        feature[3] = 8.1;
        feature[4] = 1.2;
        feature[5] = -0.7;
        feature[6] = 2.3;
        feature
    }

    fn hash_of(variant: &[f32], hash_type: HashType, nbin_dist: usize, nbin_angle: usize) -> u32 {
        GeometricHash::perfect_hash_as_u32(&variant.to_vec(), hash_type, nbin_dist, nbin_angle)
    }

    #[test]
    fn radius_one_moves_one_dimension_at_a_time() {
        let tolerance = ToleranceConfig::new(vec![0.5], vec![5.0], 0.0, 1);
        let mut expander = FeatureExpander::new(HashType::PDBTrRosetta, 16, 4, &tolerance);
        // 2 distance + 3 angle dimensions, two signs each
        assert_eq!(expander.dim_count(), 5);
        let mut count = 0;
        expander.for_each_neighbor(&pdbtr_feature(), |_| {
            count += 1;
            true
        });
        assert_eq!(count, 10);
    }

    #[test]
    fn radius_two_adds_joint_neighbors() {
        let feature = pdbtr_feature();
        let mut hashes = Vec::new();
        for radius in [1usize, 2] {
            let tolerance = ToleranceConfig::new(vec![0.5], vec![5.0], 0.0, radius);
            let mut found = HashSet::default();
            FeatureExpander::new(HashType::PDBTrRosetta, 16, 4, &tolerance)
                .for_each_neighbor(&feature, |variant| {
                    found.insert(hash_of(variant, HashType::PDBTrRosetta, 16, 4));
                    true
                });
            hashes.push(found);
        }
        // Radius 2 is a strict superset: it still walks every single-dimension
        // neighbour, then adds the joint ones.
        assert!(hashes[0].is_subset(&hashes[1]));
        assert!(hashes[1].len() > hashes[0].len());
    }

    #[test]
    fn wide_tolerance_is_substepped_and_leaves_no_gap() {
        // 16 distance bins over [2, 20] give a 1.2 A bin, so a 3.0 A tolerance has
        // to reach the bins in between and not only the ones at +-3.0 A.
        let tolerance = ToleranceConfig::new(vec![3.0], vec![], 0.0, 1);
        let mut expander = FeatureExpander::new(HashType::PDBTrRosetta, 16, 4, &tolerance);
        let feature = pdbtr_feature();
        let mut values: Vec<f32> = vec![feature[2]];
        expander.for_each_neighbor(&feature, |variant| {
            values.push(variant[2]);
            true
        });
        values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        values.dedup_by(|a, b| (*a - *b).abs() < 1e-4);
        assert!(values.len() >= 7, "expected sub-steps, got {:?}", values);
        assert!((values[0] - (feature[2] - 3.0)).abs() < 1e-3);
        assert!((values[values.len() - 1] - (feature[2] + 3.0)).abs() < 1e-3);
        // No gap wider than one bin anywhere in the covered interval
        for pair in values.windows(2) {
            assert!(pair[1] - pair[0] <= 1.2001, "gap {:?} skips a bin", pair);
        }
    }

    #[test]
    fn smaller_thresholds_reach_no_extra_bin() {
        // `-d 0.5,1.0` and `-d 1.0` have to reach the same set of bins: only the
        // widest tolerance is expanded, sub-stepped to stay gapless.
        let feature = pdbtr_feature();
        let mut sets = Vec::new();
        for thresholds in [vec![0.5, 1.0], vec![1.0]] {
            let tolerance = ToleranceConfig::new(thresholds, vec![], 0.0, 2);
            let mut found = HashSet::default();
            FeatureExpander::new(HashType::PDBTrRosetta, 16, 4, &tolerance)
                .for_each_neighbor(&feature, |variant| {
                    found.insert(hash_of(variant, HashType::PDBTrRosetta, 16, 4));
                    true
                });
            sets.push(found);
        }
        assert_eq!(sets[0], sets[1]);
    }

    #[test]
    fn elastic_tolerance_scales_with_distance() {
        let tolerance = ToleranceConfig::new(vec![0.5], vec![], 0.1, 1);
        let mut expander = FeatureExpander::new(HashType::PDBTrRosetta, 16, 4, &tolerance);

        let mut widest = |ca_dist: f32| {
            let mut feature = pdbtr_feature();
            feature[2] = ca_dist;
            feature[3] = ca_dist;
            let mut reach = 0.0f32;
            expander.for_each_neighbor(&feature, |variant| {
                reach = reach.max((variant[2] - ca_dist).abs());
                true
            });
            reach
        };
        // 0.5 + 0.1 * 5 = 1.0 against 0.5 + 0.1 * 15 = 2.0
        let short = widest(5.0);
        let long = widest(15.0);
        assert!((short - 1.0).abs() < 1e-3, "short: {}", short);
        assert!((long - 2.0).abs() < 1e-3, "long: {}", long);
    }

    #[test]
    fn torsion_offsets_wrap_around_pi() {
        // A query torsion at 179 degrees has to be able to reach -179 degrees,
        // which is 2 degrees away and not 358.
        let mut feature = pdbtr_feature();
        feature[5] = 179.0_f32.to_radians();
        let tolerance = ToleranceConfig::new(vec![], vec![5.0], 0.0, 1);
        let mut expander = FeatureExpander::new(HashType::FolddiscoDist, 32, 16, &tolerance);
        let mut wrapped = false;
        expander.for_each_neighbor(&feature, |variant| {
            // Every visited value stays inside the encodable domain
            assert!(variant[5] >= -PI - 1e-4 && variant[5] <= PI + 1e-4, "{}", variant[5]);
            if variant[5] < -3.0 {
                wrapped = true;
            }
            true
        });
        assert!(wrapped, "torsion offset past PI did not wrap");
    }

    #[test]
    fn perturbed_angles_stay_inside_their_bit_field() {
        // FolddiscoDist packs theta1/theta2 into 4 bits each and the Ca-Cb angle
        // into 3. A wide angle tolerance used to push a bin index past its field and
        // corrupt the neighbouring one; wrapping and reflecting keep every value in
        // the domain the encoding was sized for.
        let tolerance = ToleranceConfig::new(vec![], vec![60.0], 0.0, 3);
        let mut expander = FeatureExpander::new(HashType::FolddiscoDist, 32, 16, &tolerance);
        for &(ca_cb, theta1, theta2) in &[
            (0.0f32, PI, -PI), (PI, -PI, PI), (0.05, 3.10, -3.10),
        ] {
            let mut feature = pdbtr_feature();
            feature[4] = ca_cb;
            feature[5] = theta1;
            feature[6] = theta2;
            expander.for_each_neighbor(&feature, |variant| {
                assert!((0.0..=PI).contains(&variant[4]), "Ca-Cb angle {}", variant[4]);
                assert!((-PI..=PI).contains(&variant[5]), "theta1 {}", variant[5]);
                assert!((-PI..=PI).contains(&variant[6]), "theta2 {}", variant[6]);
                let hash = hash_of(variant, HashType::FolddiscoDist, 32, 16);
                // The residue pair occupies bits 21..29 and the distances 11..20;
                // an angle overflow would leak upwards into them.
                assert_eq!((hash >> 21) & 0x1FF, 8 * 20 + 3,
                    "angle overflowed into the residue field");
                assert_eq!((hash >> 11) & 0x3FF, expected_distance_bits(&feature),
                    "angle overflowed into the distance fields");
                true
            });
        }
    }

    /// Distance bits of an unperturbed FolddiscoDist hash, for overflow checks.
    fn expected_distance_bits(feature: &[f32]) -> u32 {
        (hash_of(feature, HashType::FolddiscoDist, 32, 16) >> 11) & 0x3FF
    }

    /// Hashes a query reaches for one residue pair, observed geometry included.
    fn reachable_hashes(
        feature: &[f32], hash_type: HashType, nbin_dist: usize, nbin_angle: usize,
        tolerance: &ToleranceConfig,
    ) -> HashSet<u32> {
        let mut found = HashSet::default();
        found.insert(hash_of(feature, hash_type, nbin_dist, nbin_angle));
        FeatureExpander::new(hash_type, nbin_dist, nbin_angle, tolerance)
            .for_each_neighbor(feature, |variant| {
                found.insert(hash_of(variant, hash_type, nbin_dist, nbin_angle));
                true
            });
        found
    }

    #[test]
    fn radius_two_finds_a_pair_that_crosses_two_boundaries_at_once() {
        // 16 distance bins over [2, 20] put a boundary at 2 + 1.2 * 4.5 = 7.4 A.
        // Query sits just below it in both Ca and Cb distance, target just above:
        // 0.1 A of drift, but in two dimensions at the same time.
        let mut query = pdbtr_feature();
        query[2] = 7.35;
        query[3] = 7.35;
        let mut target = query.clone();
        target[2] = 7.45;
        target[3] = 7.45;
        let target_hash = hash_of(&target, HashType::PDBTrRosetta, 16, 4);
        assert_ne!(target_hash, hash_of(&query, HashType::PDBTrRosetta, 16, 4));

        let radius_one = reachable_hashes(
            &query, HashType::PDBTrRosetta, 16, 4,
            &ToleranceConfig::new(vec![0.5], vec![5.0], 0.0, 1),
        );
        let radius_two = reachable_hashes(
            &query, HashType::PDBTrRosetta, 16, 4,
            &ToleranceConfig::new(vec![0.5], vec![5.0], 0.0, 2),
        );
        assert!(!radius_one.contains(&target_hash), "radius 1 should not reach a joint shift");
        assert!(radius_two.contains(&target_hash), "radius 2 missed the joint shift");
    }

    #[test]
    fn elastic_tolerance_finds_a_stretched_long_range_pair() {
        // A 16.9 A pair stretched to 18.6 A: 1.7 A of drift is far outside a flat
        // 0.5 A tolerance, but well inside 0.5 + 10% of the pair distance.
        let mut query = pdbtr_feature();
        query[2] = 16.9;
        query[3] = 16.9;
        let mut target = query.clone();
        target[2] = 18.6;
        target[3] = 18.6;
        let target_hash = hash_of(&target, HashType::PDBTrRosetta, 16, 4);

        let flat = reachable_hashes(
            &query, HashType::PDBTrRosetta, 16, 4,
            &ToleranceConfig::new(vec![0.5], vec![5.0], 0.0, 2),
        );
        let elastic = reachable_hashes(
            &query, HashType::PDBTrRosetta, 16, 4,
            &ToleranceConfig::new(vec![0.5], vec![5.0], 0.1, 2),
        );
        assert!(!flat.contains(&target_hash), "a flat 0.5 A tolerance should not reach 1.7 A");
        assert!(elastic.contains(&target_hash), "elastic tolerance missed the stretched pair");
        // Short pairs must stay tight: the same 10% is only 0.5 A at 5 A
        let mut short_query = pdbtr_feature();
        short_query[2] = 5.0;
        short_query[3] = 5.0;
        let mut short_target = short_query.clone();
        short_target[2] = 6.7;
        short_target[3] = 6.7;
        let short_elastic = reachable_hashes(
            &short_query, HashType::PDBTrRosetta, 16, 4,
            &ToleranceConfig::new(vec![0.5], vec![5.0], 0.1, 2),
        );
        assert!(
            !short_elastic.contains(&hash_of(&short_target, HashType::PDBTrRosetta, 16, 4)),
            "elastic tolerance must not loosen short pairs by the same absolute amount"
        );
    }

    #[test]
    fn wrapping_finds_a_torsion_across_the_pi_boundary() {
        // 179 degrees against -179 degrees: 2 degrees apart, opposite ends of the
        // encoded range. Without wrapping no angle tolerance can bridge them.
        let mut query = pdbtr_feature();
        query[5] = 179.0_f32.to_radians();
        let mut target = query.clone();
        target[5] = -179.0_f32.to_radians();
        let target_hash = hash_of(&target, HashType::FolddiscoDist, 32, 16);
        assert_ne!(target_hash, hash_of(&query, HashType::FolddiscoDist, 32, 16));

        let reachable = reachable_hashes(
            &query, HashType::FolddiscoDist, 32, 16,
            &ToleranceConfig::new(vec![0.5], vec![5.0], 0.0, 1),
        );
        assert!(reachable.contains(&target_hash), "wrapped torsion neighbour not reached");
    }

    #[test]
    fn visit_can_stop_early() {
        let tolerance = ToleranceConfig::new(vec![0.5], vec![5.0], 0.0, 3);
        let mut expander = FeatureExpander::new(HashType::PDBTrRosetta, 16, 4, &tolerance);
        let mut count = 0;
        expander.for_each_neighbor(&pdbtr_feature(), |_| {
            count += 1;
            count < 4
        });
        assert_eq!(count, 4);
    }

    #[test]
    fn wrap_and_reflect_helpers() {
        assert!((wrap_to_pi(PI + 0.2) - (-PI + 0.2)).abs() < 1e-5);
        assert!((wrap_to_pi(-PI - 0.2) - (PI - 0.2)).abs() < 1e-5);
        assert!((wrap_to_pi(0.3) - 0.3).abs() < 1e-6);
        assert!((reflect_into_range(-0.2, 0.0, PI) - 0.2).abs() < 1e-5);
        assert!((reflect_into_range(PI + 0.2, 0.0, PI) - (PI - 0.2)).abs() < 1e-5);
        assert!((reflect_into_range(1.0, 0.0, PI) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn empty_tolerance_produces_nothing() {
        let mut expander = FeatureExpander::new(
            HashType::PDBTrRosetta, 16, 4, &ToleranceConfig::none()
        );
        assert!(expander.is_noop());
        let mut count = 0;
        expander.for_each_neighbor(&pdbtr_feature(), |_| { count += 1; true });
        assert_eq!(count, 0);
    }

    #[test]
    fn every_hash_type_keeps_wide_angle_tolerance_inside_its_encoding() {
        // Sanity net for the per-type angle metadata: whatever a wide angle
        // tolerance does to a feature, the packed hash stays inside the declared
        // width. Distances are kept mid-range on purpose — those mirror the index
        // encoding out of the window as well, so they are not bounded here.
        let tolerance = ToleranceConfig::new(vec![1.0], vec![60.0], 0.05, 2);
        for hash_type in [
            HashType::PDBMotif, HashType::PDBMotifSinCos, HashType::TrRosetta,
            HashType::PDBTrRosetta, HashType::PointPairFeature,
            HashType::TertiaryInteraction, HashType::Hybrid,
            HashType::FolddiscoAngle, HashType::FolddiscoDist,
        ] {
            let nbin_dist = hash_type.default_dist_bin();
            let nbin_angle = hash_type.default_angle_bin();
            let mut expander = FeatureExpander::new(hash_type, nbin_dist, nbin_angle, &tolerance);
            let mut feature = pdbtr_feature();
            if hash_type == HashType::PDBMotif {
                feature[4] = 120.0; // this one keeps angles in degrees
            }
            if hash_type == HashType::TertiaryInteraction {
                feature[7] = 12.0;  // Ca distance
                feature[8] = 3.0;   // sequence separation
            }
            let bits = hash_type.encoding_bits() as u32;
            expander.for_each_neighbor(&feature, |variant| {
                let hash = hash_of(variant, hash_type, nbin_dist, nbin_angle);
                if bits < 32 {
                    assert_eq!(hash >> bits, 0, "{:?} exceeded {} bits", hash_type, bits);
                }
                true
            });
        }
    }
}
