// File: core.rs
// Author: Hyunbin Kim (khb7840@gmail.com)
// Description: Core geometric hash enum and types

use std::{fmt, hash::Hash, io::{BufRead, Write}};

use crate::utils::traits::HashableSync;

#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub enum HashType {
    PDBMotif,
    PDBMotifSinCos,
    TrRosetta,
    PDBTrRosetta,
    PointPairFeature,
    TertiaryInteraction,
    Hybrid,
    FolddiscoAngle,
    FolddiscoDist,
    // append new hash type here
    Other,
}

impl HashType {
    
    pub fn get_with_index(index: usize) -> Self {
        match index {
            0 => HashType::PDBMotif,
            1 => HashType::PDBMotifSinCos,
            2 => HashType::TrRosetta,
            3 => HashType::PDBTrRosetta,
            4 => HashType::PointPairFeature,
            5 => HashType::TertiaryInteraction,
            6 => HashType::Hybrid,
            7 => HashType::FolddiscoAngle,
            8 => HashType::FolddiscoDist,
            // append new hash type here
            _ => HashType::Other,
        }
    }
    
    pub fn get_with_str(hash_type: &str) -> Self {
        match hash_type {
            "0" | "PDBMotif" | "pyscomotif" | "orig_pdb" => HashType::PDBMotif,
            "1" | "PDBMotifSinCos" | "pdb" => HashType::PDBMotifSinCos,
            "2" | "TrRosetta" | "trrosetta" | "tr" => HashType::TrRosetta,
            "3" | "PDBTrRosetta" | "pdbtr" | "default" | "folddisco" => HashType::PDBTrRosetta,
            "4" | "PointPairFeature" | "ppf" => HashType::PointPairFeature,
            "5" | "TertiaryInteraction" | "tertiary" | "3di" => HashType::TertiaryInteraction,
            "6" | "Hybrid" | "hybrid" => HashType::Hybrid,
            "7" | "FolddiscoAngle"| "angle" | "folddisco_angle" => HashType::FolddiscoAngle,
            "8" | "FolddiscoDist" | "distance" | "dist" | "folddisco_dist" => HashType::FolddiscoDist,
            // append new hash type here
            _ => HashType::Other,
        }
    }
    
    pub fn to_string(&self) -> String {
        match self {
            HashType::PDBMotif => "PDBMotif".to_string(),
            HashType::PDBMotifSinCos => "PDBMotifSinCos".to_string(),
            HashType::TrRosetta => "TrRosetta".to_string(),
            HashType::PDBTrRosetta => "PDBTrRosetta".to_string(),
            HashType::PointPairFeature => "PointPairFeature".to_string(),
            HashType::TertiaryInteraction => "TertiaryInteraction".to_string(),
            HashType::Hybrid => "Hybrid".to_string(),
            HashType::FolddiscoAngle => "FolddiscoAngle".to_string(),
            HashType::FolddiscoDist => "FolddiscoDist".to_string(),
            // append new hash type here
            HashType::Other => "Other".to_string(),
        }
    }

    pub fn encoding_type(&self) -> usize {
        // Unified to u32 encoding
        32usize
    }

    pub fn encoding_bits(&self) -> usize {
        match self {
            HashType::PDBMotif => 25usize,
            HashType::PDBMotifSinCos => 26usize,
            HashType::TrRosetta => 32usize,
            HashType::PDBTrRosetta => 30usize,
            HashType::PointPairFeature => 32usize,
            HashType::TertiaryInteraction => 29usize,
            HashType::Hybrid => 32usize,
            HashType::FolddiscoAngle => 32usize,
            HashType::FolddiscoDist => 32usize,
            // append new hash type here
            HashType::Other => 32usize,
        }
    }
    
    pub fn save_to_file(&self, path: &str) {
        let mut file = std::fs::File::create(path).unwrap();
        file.write_all(format!("{:?}", self).as_bytes()).unwrap();
    }

    pub fn load_from_file(path: &str) -> Self {
        let file = std::fs::File::open(path).unwrap();
        let reader = std::io::BufReader::new(file);
        let mut hash_type = HashType::PDBTrRosetta;
        for line in reader.lines() {
            let line = line.unwrap();
            hash_type = match line.as_str() {
                "PDBMotif" => HashType::PDBMotif,
                "PDBMotifSinCos" => HashType::PDBMotifSinCos,
                "TrRosetta" => HashType::TrRosetta,
                "PDBTrRosetta" => HashType::PDBTrRosetta,
                "PointPairFeature" => HashType::PointPairFeature,
                "TertiaryInteraction" => HashType::TertiaryInteraction,
                "Hybrid" => HashType::Hybrid,
                "FolddiscoAngle" => HashType::FolddiscoAngle,
                "FolddiscoDist" => HashType::FolddiscoDist,
                // append new hash type here
                _ => HashType::Other,
            };
        }
        hash_type
    }
    
    pub fn default_dist_bin(&self) -> usize {
        match self {
            HashType::PDBMotif => super::pdb_motif::NBIN_DIST as usize,
            HashType::PDBMotifSinCos | HashType::TrRosetta | HashType::PointPairFeature |
            HashType::TertiaryInteraction => crate::utils::convert::NBIN_DIST as usize,
            HashType::PDBTrRosetta => super::pdb_tr::PDBTR_NBIN_DIST as usize,
            HashType::Hybrid => super::hybrid::HYBRID_NBIN_DIST as usize,
            HashType::FolddiscoAngle => super::folddisco_angle::NBIN_DIST as usize,
            HashType::FolddiscoDist => super::folddisco_dist::NBIN_DIST as usize,
            // append new hash type here
            HashType::Other => 0,
        }
    }
    
    pub fn default_angle_bin(&self) -> usize {
        match self {
            HashType::PDBMotif => super::pdb_motif::NBIN_ANGLE as usize,
            HashType::PDBMotifSinCos | HashType::TrRosetta | HashType::PointPairFeature |
            HashType::TertiaryInteraction => crate::utils::convert::NBIN_SIN_COS as usize,
            HashType::PDBTrRosetta => super::pdb_tr::PDBTR_NBIN_SIN_COS as usize,
            HashType::Hybrid => super::hybrid::HYBRID_NBIN_SIN_COS as usize,
            HashType::FolddiscoAngle => super::folddisco_angle::NBIN_ANGLE_360 as usize,
            HashType::FolddiscoDist => super::folddisco_dist::NBIN_ANGLE_360 as usize,
            // append new hash type here
            HashType::Other => 0,
        }
    }

    /// Widest bin count the bit layout of this hash type can hold. `perfect_hash`
    /// clamps its arguments to these, so anything reasoning about bin boundaries has
    /// to use the same numbers.
    pub fn max_dist_bin(&self) -> usize {
        match self {
            HashType::PDBMotif => super::pdb_motif::MAX_NBIN_DIST as usize,
            HashType::PDBMotifSinCos => super::pdb_motif_sincos::MAX_NBIN_DIST as usize,
            HashType::TrRosetta => super::trrosetta::MAX_NBIN_DIST as usize,
            HashType::PDBTrRosetta => super::pdb_tr::PDBTR_MAX_NBIN_DIST as usize,
            HashType::PointPairFeature => super::ppf::MAX_NBIN_DIST as usize,
            HashType::TertiaryInteraction => super::tertiary_interaction::MAX_NBIN_DIST as usize,
            HashType::Hybrid => super::hybrid::HYBRID_MAX_NBIN_DIST as usize,
            // These two clamp against their own default bin count
            HashType::FolddiscoAngle => super::folddisco_angle::NBIN_DIST as usize,
            HashType::FolddiscoDist => super::folddisco_dist::NBIN_DIST as usize,
            // append new hash type here
            HashType::Other => 0,
        }
    }

    pub fn max_angle_bin(&self) -> usize {
        match self {
            HashType::PDBMotif => super::pdb_motif::MAX_NBIN_ANGLE as usize,
            HashType::PDBMotifSinCos => super::pdb_motif_sincos::MAX_NBIN_SIN_COS as usize,
            HashType::TrRosetta => super::trrosetta::MAX_NBIN_SIN_COS as usize,
            HashType::PDBTrRosetta => super::pdb_tr::PDBTR_MAX_NBIN_SIN_COS as usize,
            HashType::PointPairFeature => super::ppf::MAX_NBIN_SIN_COS as usize,
            HashType::TertiaryInteraction => super::tertiary_interaction::MAX_NBIN_SIN_COS as usize,
            HashType::Hybrid => super::hybrid::HYBRID_MAX_NBIN_SIN_COS as usize,
            // These two clamp against their own default bin count
            HashType::FolddiscoAngle => super::folddisco_angle::NBIN_ANGLE_360 as usize,
            HashType::FolddiscoDist => super::folddisco_dist::NBIN_ANGLE_360 as usize,
            // append new hash type here
            HashType::Other => 0,
        }
    }

    /// Bin count actually used for a requested one, mirroring `perfect_hash`
    /// (0 means "use the default for this hash type").
    pub fn effective_dist_bin(&self, nbin_dist: usize) -> usize {
        if nbin_dist == 0 { self.default_dist_bin() } else { nbin_dist.min(self.max_dist_bin()) }
    }

    pub fn effective_angle_bin(&self, nbin_angle: usize) -> usize {
        if nbin_angle == 0 { self.default_angle_bin() } else { nbin_angle.min(self.max_angle_bin()) }
    }

    /// Distance window this hash type discretizes over. `PDBMotif` carries its own
    /// copy of the constants, the rest share the ones in `utils::convert`.
    pub fn dist_range(&self) -> (f32, f32) {
        match self {
            HashType::PDBMotif => (super::pdb_motif::MIN_DIST, super::pdb_motif::MAX_DIST),
            // append new hash type here if it uses its own window
            _ => (crate::utils::convert::MIN_DIST, crate::utils::convert::MAX_DIST),
        }
    }

    /// Width of one distance bin, in Angstroms.
    pub fn dist_bin_width(&self, nbin_dist: usize) -> f32 {
        let nbin = self.effective_dist_bin(nbin_dist);
        if nbin < 2 {
            return f32::INFINITY;
        }
        let (min_dist, max_dist) = self.dist_range();
        (max_dist - min_dist) / (nbin as f32 - 1.0)
    }

    /// Angular step guaranteed not to skip over a bin, in the unit this hash type
    /// stores angles in.
    ///
    /// The types that discretize the angle value directly get the plain bin width.
    /// Sin-cos encoded types bin `sin` and `cos` over `[MIN_SIN_COS, MAX_SIN_COS]`
    /// instead; since `|d sin/d theta| <= 1`, an angular step no wider than one
    /// sin-cos bin cannot jump past a bin there either.
    pub fn angle_bin_width(&self, nbin_angle: usize) -> f32 {
        let nbin = self.effective_angle_bin(nbin_angle);
        if nbin < 2 {
            return f32::INFINITY;
        }
        let span = match self {
            HashType::PDBMotif => super::pdb_motif::MAX_ANGLE - super::pdb_motif::MIN_ANGLE,
            HashType::FolddiscoAngle => {
                super::folddisco_angle::MAX_ANGLE_RAD - super::folddisco_angle::MIN_ANGLE_RAD
            }
            HashType::FolddiscoDist => {
                super::folddisco_dist::MAX_ANGLE_RAD - super::folddisco_dist::MIN_ANGLE_RAD
            }
            // sin-cos encoded types
            _ => crate::utils::convert::MAX_SIN_COS - crate::utils::convert::MIN_SIN_COS,
        };
        span / (nbin as f32 - 1.0)
    }

    /// True when the angle features of this hash type are stored in degrees.
    /// Everything except the original `PDBMotif` encoding works in radians.
    pub fn angle_in_degrees(&self) -> bool {
        matches!(self, HashType::PDBMotif)
    }

    /// Feature dimensions holding a periodic torsion, i.e. an `atan2` value on
    /// `[-PI, PI]`. A tolerance offset that runs past a bound belongs at the other
    /// end of the range rather than clamped against it.
    ///
    /// This only changes anything for the types that bin the angle value directly
    /// (`FolddiscoAngle`, `FolddiscoDist`), where 179 degrees is the top bin and
    /// -179 the bottom one, so an unwrapped +5 degree offset runs off the end of the
    /// field instead of reaching the neighbour 2 degrees away. Sin-cos encoded types
    /// are periodic already -- `sin(184 deg) == sin(-176 deg)` -- so wrapping their
    /// value first is a harmless no-op that keeps the handling uniform.
    pub fn periodic_angle_index(&self) -> Option<Vec<usize>> {
        match self {
            // omega, theta1, theta2
            HashType::TrRosetta => Some(vec![3, 4, 5]),
            // theta1, theta2
            HashType::PDBTrRosetta | HashType::FolddiscoAngle | HashType::FolddiscoDist => Some(vec![5, 6]),
            // theta1, theta2 and the two backbone torsions
            HashType::Hybrid => Some(vec![5, 6, 7, 8]),
            // append new hash type here
            _ => None,
        }
    }

    /// Angle dimensions confined to a range, together with that range. These come
    /// out of `acos`, so a tolerance offset crossing a bound is reflected back into
    /// the range: one radian below zero is one radian above it.
    pub fn bounded_angle_index(&self) -> Option<Vec<(usize, f32, f32)>> {
        const PI: f32 = std::f32::consts::PI;
        match self {
            HashType::PDBMotif => Some(vec![
                (4, super::pdb_motif::MIN_ANGLE, super::pdb_motif::MAX_ANGLE)
            ]),
            // Ca-Cb angle
            HashType::PDBMotifSinCos | HashType::PDBTrRosetta | HashType::Hybrid |
            HashType::FolddiscoAngle | HashType::FolddiscoDist => Some(vec![(4, 0.0, PI)]),
            // phi1, phi2
            HashType::TrRosetta => Some(vec![(6, 0.0, PI), (7, 0.0, PI)]),
            HashType::PointPairFeature => Some(vec![(3, 0.0, PI), (4, 0.0, PI), (5, 0.0, PI)]),
            HashType::TertiaryInteraction => Some((0..=6).map(|i| (i, 0.0, PI)).collect()),
            // append new hash type here
            _ => None,
        }
    }

    /// Pull the angle dimensions of a perturbed feature vector back into the domain
    /// they are defined on: torsions wrap at +-PI, bounded angles reflect at their
    /// ends. Without this a wide angle tolerance pushes the bin index past the end of
    /// its bit field on the types that discretize the angle directly, which corrupts
    /// the neighbouring field, and on sin-cos types it asks for a `sin` sign no real
    /// structure can produce.
    ///
    /// Distances are deliberately left alone. A structure can hold a Cb-Cb distance
    /// beyond `MAX_DIST` -- only the Ca distance is checked against the cutoff -- so
    /// the index itself contains those out-of-window encodings, aliasing included.
    /// Indexing and querying run through the same discretizer, so mirroring it
    /// exactly is what makes a perturbed query hash mean the same thing as an index
    /// hash; clamping here would be a *different* encoding and would stop the
    /// tolerance from reaching long pairs at all.
    pub fn sanitize_perturbed_feature(&self, feature: &mut [f32]) {
        if let Some(periodic) = self.periodic_angle_index() {
            for idx in periodic {
                feature[idx] = crate::utils::convert::wrap_to_pi(feature[idx]);
            }
        }
        if let Some(bounded) = self.bounded_angle_index() {
            for (idx, lo, hi) in bounded {
                feature[idx] = crate::utils::convert::reflect_into_range(feature[idx], lo, hi);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_hash_type() {
        let path = "data/test.type";
        let hash_type_vec = vec![
            HashType::PDBMotif,
            HashType::PDBMotifSinCos,
            HashType::TrRosetta,
            HashType::PDBTrRosetta,
            HashType::PointPairFeature,
            HashType::TertiaryInteraction,
            HashType::Hybrid,
            HashType::FolddiscoAngle,
            HashType::FolddiscoDist,
            // append new hash type here
        ];
        for hash_type in hash_type_vec {
            hash_type.save_to_file(path);
            let loaded_hash_type = HashType::load_from_file(path);
            assert_eq!(hash_type, loaded_hash_type);
        }
    }
}

#[derive(Clone, Copy, Hash, Eq, PartialEq, Ord, PartialOrd)]
pub enum GeometricHash {
    PDBMotif(super::pdb_motif::HashValue),
    PDBMotifSinCos(super::pdb_motif_sincos::HashValue),
    TrRosetta(super::trrosetta::HashValue),
    PDBTrRosetta(super::pdb_tr::HashValue),
    PointPairFeature(super::ppf::HashValue),
    TertiaryInteraction(super::tertiary_interaction::HashValue),
    Hybrid(super::hybrid::HashValue),
    FolddiscoAngle(super::folddisco_angle::HashValue),
    FolddiscoDist(super::folddisco_dist::HashValue),
    // append new hash type here
}

impl HashableSync for GeometricHash {}

impl GeometricHash {
    pub fn perfect_hash_default_as_u32(feature: &Vec<f32>, hash_type: HashType) -> u32 {
        match hash_type {
            HashType::PDBMotif => super::pdb_motif::HashValue::perfect_hash_default(feature),
            HashType::PDBMotifSinCos => super::pdb_motif_sincos::HashValue::perfect_hash_default(feature),
            HashType::TrRosetta => super::trrosetta::HashValue::perfect_hash_default(feature),
            HashType::PDBTrRosetta => super::pdb_tr::HashValue::perfect_hash_default(feature),
            HashType::PointPairFeature => super::ppf::HashValue::perfect_hash_default(feature),
            HashType::TertiaryInteraction => super::tertiary_interaction::HashValue::perfect_hash_default(feature),
            HashType::Hybrid => super::hybrid::HashValue::perfect_hash_default(feature),
            HashType::FolddiscoAngle => super::folddisco_angle::HashValue::perfect_hash_default(feature),
            HashType::FolddiscoDist => super::folddisco_dist::HashValue::perfect_hash_default(feature),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }

    pub fn perfect_hash_as_u32(
        feature: &Vec<f32>, hash_type: HashType, nbin_dist: usize, nbin_angle: usize
    ) -> u32 {
        match hash_type {
            HashType::PDBMotif => super::pdb_motif::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            HashType::PDBMotifSinCos => super::pdb_motif_sincos::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            HashType::TrRosetta => super::trrosetta::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            HashType::PDBTrRosetta => super::pdb_tr::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            HashType::PointPairFeature => super::ppf::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            HashType::TertiaryInteraction => super::tertiary_interaction::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            HashType::Hybrid => super::hybrid::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            HashType::FolddiscoAngle => super::folddisco_angle::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            HashType::FolddiscoDist => super::folddisco_dist::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }

    pub fn perfect_hash_with_shifts_dedup_inline(feature: &Vec<f32>, hash_type: HashType) -> (u8, [u32; 8]) {
        match hash_type {
            HashType::PDBTrRosetta => super::pdb_tr::HashValue::perfect_hash_with_shifts_dedup_inline(feature),
            // Use exhaustive deduplication for now
            // HashType::PDBTrRosetta => super::pdb_tr::HashValue::perfect_hash_with_all_shifts_exhaustive(feature),
            _ => panic!("Hash type does not support shift deduplication"),
        }
    }
    
    pub fn perfect_hash_default(feature: &Vec<f32>, hash_type: HashType) -> Self {
        match hash_type {
            HashType::PDBMotif => GeometricHash::PDBMotif(
                super::pdb_motif::HashValue(
                    super::pdb_motif::HashValue::perfect_hash_default(feature)
                )
            ),
            HashType::PDBMotifSinCos => GeometricHash::PDBMotifSinCos(
                super::pdb_motif_sincos::HashValue(
                    super::pdb_motif_sincos::HashValue::perfect_hash_default(feature)
                )
            ),
            HashType::TrRosetta => GeometricHash::TrRosetta(
                super::trrosetta::HashValue(
                    super::trrosetta::HashValue::perfect_hash_default(feature)
                )
            ),
            HashType::PDBTrRosetta => GeometricHash::PDBTrRosetta(
                super::pdb_tr::HashValue(
                    super::pdb_tr::HashValue::perfect_hash_default(feature)
                )
            ),
            HashType::PointPairFeature => GeometricHash::PointPairFeature(
                super::ppf::HashValue(
                    super::ppf::HashValue::perfect_hash_default(feature)
                )
            ),
            HashType::TertiaryInteraction => GeometricHash::TertiaryInteraction(
                super::tertiary_interaction::HashValue(
                    super::tertiary_interaction::HashValue::perfect_hash_default(feature)
                )
            ),
            HashType::Hybrid => GeometricHash::Hybrid(
                super::hybrid::HashValue(
                    super::hybrid::HashValue::perfect_hash_default(feature)
                )
            ),
            HashType::FolddiscoAngle => GeometricHash::FolddiscoAngle(
                super::folddisco_angle::HashValue(
                    super::folddisco_angle::HashValue::perfect_hash_default(feature)
                )
            ),
            HashType::FolddiscoDist => GeometricHash::FolddiscoDist(
                super::folddisco_dist::HashValue(
                    super::folddisco_dist::HashValue::perfect_hash_default(feature)
                )
            ),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }

    pub fn perfect_hash(
        feature: &Vec<f32>, hash_type: HashType, nbin_dist: usize, nbin_angle: usize
    ) -> Self {
        match hash_type {
            HashType::PDBMotif => GeometricHash::PDBMotif(
                super::pdb_motif::HashValue(
                    super::pdb_motif::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            HashType::PDBMotifSinCos => GeometricHash::PDBMotifSinCos(
                super::pdb_motif_sincos::HashValue(
                    super::pdb_motif_sincos::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            HashType::TrRosetta => GeometricHash::TrRosetta(
                super::trrosetta::HashValue(
                    super::trrosetta::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            HashType::PDBTrRosetta => GeometricHash::PDBTrRosetta(
                super::pdb_tr::HashValue(
                    super::pdb_tr::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            HashType::PointPairFeature => GeometricHash::PointPairFeature(
                super::ppf::HashValue(
                    super::ppf::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            HashType::TertiaryInteraction => GeometricHash::TertiaryInteraction(
                super::tertiary_interaction::HashValue(
                    super::tertiary_interaction::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            HashType::Hybrid => GeometricHash::Hybrid(
                super::hybrid::HashValue(
                    super::hybrid::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            HashType::FolddiscoAngle => GeometricHash::FolddiscoAngle(
                super::folddisco_angle::HashValue(
                    super::folddisco_angle::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            HashType::FolddiscoDist => GeometricHash::FolddiscoDist(
                super::folddisco_dist::HashValue(
                    super::folddisco_dist::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }

    pub fn reverse_hash_default(&self, output: &mut Vec<f32>) {
        match self {
            GeometricHash::PDBMotif(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::PDBMotifSinCos(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::TrRosetta(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::PDBTrRosetta(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            }
            GeometricHash::PointPairFeature(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::TertiaryInteraction(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::Hybrid(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::FolddiscoAngle(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::FolddiscoDist(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            // append new hash type here
            // _ => panic!("Invalid hash type"),
        }
    }


    pub fn reverse_hash(&self, nbin_dist: usize, nbin_angle: usize, output: &mut Vec<f32>) {
        match self {
            GeometricHash::PDBMotif(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::PDBMotifSinCos(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::TrRosetta(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::PDBTrRosetta(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::PointPairFeature(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::TertiaryInteraction(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::Hybrid(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::FolddiscoAngle(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::FolddiscoDist(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            // append new hash type here
            // _ => panic!("Invalid hash type"),
        }
    }


    pub fn hash_type(&self) -> HashType {
        match self {
            GeometricHash::PDBMotif(hash) => hash.hash_type(),
            GeometricHash::PDBMotifSinCos(hash) => hash.hash_type(),
            GeometricHash::TrRosetta(hash) => hash.hash_type(),
            GeometricHash::PointPairFeature(hash) => hash.hash_type(),
            GeometricHash::PDBTrRosetta(hash) => hash.hash_type(),
            GeometricHash::TertiaryInteraction(hash) => hash.hash_type(),
            GeometricHash::Hybrid(hash) => hash.hash_type(),
            GeometricHash::FolddiscoAngle(hash) => hash.hash_type(),
            GeometricHash::FolddiscoDist(hash) => hash.hash_type(),
            // append new hash type here
            // _ => panic!("Invalid hash type"),
        }
    }


    pub fn from_u32(hashvalue: u32, hash_type: HashType) -> Self {
        match hash_type {
            HashType::PDBMotif => GeometricHash::PDBMotif(
                super::pdb_motif::HashValue::from_u32(hashvalue)
            ),
            HashType::PDBMotifSinCos => GeometricHash::PDBMotifSinCos(
                super::pdb_motif_sincos::HashValue::from_u32(hashvalue)
            ),
            HashType::TrRosetta => GeometricHash::TrRosetta(
                super::trrosetta::HashValue::from_u32(hashvalue)
            ),
            HashType::PDBTrRosetta => GeometricHash::PDBTrRosetta(
                super::pdb_tr::HashValue::from_u32(hashvalue)
            ),
            HashType::PointPairFeature => GeometricHash::PointPairFeature(
                super::ppf::HashValue::from_u32(hashvalue)
            ),
            HashType::TertiaryInteraction => GeometricHash::TertiaryInteraction(
                super::tertiary_interaction::HashValue::from_u32(hashvalue)
            ),
            HashType::Hybrid => GeometricHash::Hybrid(
                super::hybrid::HashValue::from_u32(hashvalue)
            ),
            HashType::FolddiscoAngle => GeometricHash::FolddiscoAngle(
                super::folddisco_angle::HashValue::from_u32(hashvalue)
            ),
            HashType::FolddiscoDist => GeometricHash::FolddiscoDist(
                super::folddisco_dist::HashValue::from_u32(hashvalue)
            ),
            // append new hash type here if it is encoded as u32
            _ => panic!("Invalid hash type"),
        }
    }
    
    pub fn from_u64(hashvalue: u64, hash_type: HashType) -> Self {
        match hash_type {
            HashType::PDBMotif => GeometricHash::PDBMotif(
                super::pdb_motif::HashValue::from_u64(hashvalue)
            ),
            HashType::PDBMotifSinCos => GeometricHash::PDBMotifSinCos(
                super::pdb_motif_sincos::HashValue::from_u64(hashvalue)
            ),
            HashType::TrRosetta => GeometricHash::TrRosetta(
                super::trrosetta::HashValue::from_u64(hashvalue)
            ),
            HashType::PDBTrRosetta => GeometricHash::PDBTrRosetta(
                super::pdb_tr::HashValue::from_u64(hashvalue)
            ),
            HashType::PointPairFeature => GeometricHash::PointPairFeature(
                super::ppf::HashValue::from_u64(hashvalue)
            ),
            HashType::TertiaryInteraction => GeometricHash::TertiaryInteraction(
                super::tertiary_interaction::HashValue::from_u64(hashvalue)
            ),
            HashType::Hybrid => GeometricHash::Hybrid(
                super::hybrid::HashValue::from_u64(hashvalue)
            ),
            HashType::FolddiscoAngle => GeometricHash::FolddiscoAngle(
                super::folddisco_angle::HashValue::from_u64(hashvalue)
            ),
            HashType::FolddiscoDist => GeometricHash::FolddiscoDist(
                super::folddisco_dist::HashValue::from_u64(hashvalue)
            ),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }
    
    pub fn as_u32(&self) -> u32 {
        match self {
            GeometricHash::PDBMotif(hash) => hash.as_u32(),
            GeometricHash::PDBMotifSinCos(hash) => hash.as_u32(),
            GeometricHash::TrRosetta(hash) => hash.as_u32(),
            GeometricHash::PDBTrRosetta(hash) => hash.as_u32(),
            GeometricHash::PointPairFeature(hash) => hash.as_u32(),
            GeometricHash::TertiaryInteraction(hash) => hash.as_u32(),
            GeometricHash::Hybrid(hash) => hash.as_u32(),
            GeometricHash::FolddiscoAngle(hash) => hash.as_u32(),
            GeometricHash::FolddiscoDist(hash) => hash.as_u32(),
            // append new hash type here
        }
    }
    pub fn as_u64(&self) -> u64 {
        match self {
            GeometricHash::PDBMotif(hash) => hash.as_u64(),
            GeometricHash::PDBMotifSinCos(hash) => hash.as_u64(),
            GeometricHash::TrRosetta(hash) => hash.as_u64(),
            GeometricHash::PDBTrRosetta(hash) => hash.as_u64(),
            GeometricHash::PointPairFeature(hash) => hash.as_u64(),
            GeometricHash::TertiaryInteraction(hash) => hash.as_u64(),
            GeometricHash::Hybrid(hash) => hash.as_u64(),
            GeometricHash::FolddiscoAngle(hash) => hash.as_u64(),
            GeometricHash::FolddiscoDist(hash) => hash.as_u64(),
            // append new hash type here
        }
    }
    
    pub fn is_symmetric(&self) -> bool {
        match self {
            GeometricHash::PDBMotif(hash) => hash.is_symmetric(),
            GeometricHash::PDBMotifSinCos(hash) => hash.is_symmetric(),
            GeometricHash::TrRosetta(hash) => hash.is_symmetric(),
            GeometricHash::PDBTrRosetta(hash) => hash.is_symmetric(),
            GeometricHash::PointPairFeature(hash) => hash.is_symmetric(),
            GeometricHash::TertiaryInteraction(hash) => hash.is_symmetric(),
            GeometricHash::Hybrid(hash) => hash.is_symmetric(),
            GeometricHash::FolddiscoAngle(hash) => hash.is_symmetric(),
            GeometricHash::FolddiscoDist(hash) => hash.is_symmetric(),
            // append new hash type here
        }
    }
    
    pub fn downcast_pdb_motif(&self) -> super::pdb_motif::HashValue {
        match self {
            GeometricHash::PDBMotif(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_pdb_motif_sincos(&self) -> super::pdb_motif_sincos::HashValue {
        match self {
            GeometricHash::PDBMotifSinCos(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_default_32bit(&self) -> super::trrosetta::HashValue {
        match self {
            GeometricHash::TrRosetta(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_point_pair_feature(&self) -> super::ppf::HashValue {
        match self {
            GeometricHash::PointPairFeature(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_pdb_tr(&self) -> super::pdb_tr::HashValue {
        match self {
            GeometricHash::PDBTrRosetta(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_tertiary_interaction(&self) -> super::tertiary_interaction::HashValue {
        match self {
            GeometricHash::TertiaryInteraction(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_hybrid(&self) -> super::hybrid::HashValue {
        match self {
            GeometricHash::Hybrid(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_folddisco_angle(&self) -> super::folddisco_angle::HashValue {
        match self {
            GeometricHash::FolddiscoAngle(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_folddisco_dist(&self) -> super::folddisco_dist::HashValue {
        match self {
            GeometricHash::FolddiscoDist(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    // append the downcast method for new hash type here

}

impl fmt::Debug for GeometricHash {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GeometricHash::PDBMotif(hash) => {
                write!(f, "PDBMotif({:?})", hash)
            },
            GeometricHash::PDBMotifSinCos(hash) => {
                write!(f, "PDBMotifSinCos({:?})", hash)
            },
            GeometricHash::TrRosetta(hash) => {
                write!(f, "TrRosetta({:?})", hash)
            },
            GeometricHash::PDBTrRosetta(hash) => {
                write!(f, "PDBTrRosetta({:?})", hash)
            },
            GeometricHash::PointPairFeature(hash) => {
                write!(f, "PointPairFeature({:?})", hash)
            },
            GeometricHash::TertiaryInteraction(hash) => {
                write!(f, "TertiaryInteraction({:?})", hash)
            },
            GeometricHash::Hybrid(hash) => {
                write!(f, "Hybrid({:?})", hash)
            },
            GeometricHash::FolddiscoAngle(hash) => {
                write!(f, "FolddiscoAngle({:?})", hash)
            },
            GeometricHash::FolddiscoDist(hash) => {
                write!(f, "FolddiscoDist({:?})", hash)
            },  
            // append new hash type here
            // _ => panic!("Invalid hash type"),
        }
    }
}

impl fmt::Display for GeometricHash {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GeometricHash::PDBMotif(hash) => {
                write!(f, "PDBMotif\t{:?}", hash)
            },
            GeometricHash::PDBMotifSinCos(hash) => {
                write!(f, "PDBMotifSinCos\t{:?}", hash)
            },
            GeometricHash::TrRosetta(hash) => {
                write!(f, "TrRosetta\t{:?}", hash)
            },
            GeometricHash::PDBTrRosetta(hash) => {
                write!(f, "PDBTrRosetta\t{:?}", hash)
            },
            GeometricHash::PointPairFeature(hash) => {
                write!(f, "PointPairFeature\t{:?}", hash)
            },
            GeometricHash::TertiaryInteraction(hash) => {
                write!(f, "TertiaryInteraction\t{:?}", hash)
            },
            GeometricHash::Hybrid(hash) => {
                write!(f, "Hybrid\t{:?}", hash)
            },
            GeometricHash::FolddiscoAngle(hash) => {
                write!(f, "FolddiscoAngle\t{:?}", hash)
            },
            GeometricHash::FolddiscoDist(hash) => {
                write!(f, "FolddiscoDist\t{:?}", hash)
            },
            // append new hash type here
            // _ => panic!("Invalid hash type"),
        }
    }
}