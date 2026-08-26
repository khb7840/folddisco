//! Torsion-angle ENM/NMA conformational sampling for non-rigid query search.
//!
//! This module builds a residue-level elastic network in torsion space (phi/psi),
//! samples low-frequency torsional normal modes, and applies the perturbations by
//! rotating downstream backbone atoms around peptide bond axes.
//!
//! Rotating about the N-CA and CA-C axes leaves both pivot atoms on the axis, so
//! every bonded distance is preserved exactly - the conformers are valid backbones
//! by construction and there is nothing to validate afterwards. What does need care
//! is where a rotation stops: chains are not bonded to each other, so a torsion in
//! one chain must not displace another, and the displacement of every rotation is
//! rescaled to the backbone RMSD the caller asked for.

use nalgebra::{DMatrix, DVector, SymmetricEigen};
use rayon::prelude::*;

use crate::structure::coordinate::{approx_cb, Coordinate};
use crate::structure::core::CompactStructure;

// Provenance of the tuned values below. Only two of them have a source outside this
// module; the rest were carried over from the branch that introduced the torsion ENM
// and were never swept, so treat them as arbitrary-but-working starting points. Saying
// so is the point: a constant nobody can re-derive is a constant nobody can change.

/// Contact cutoff of the elastic network, in Angstroms. Matches the conventional
/// anisotropic-network-model cutoff (Atilgan et al. 2001 use 12-15 A between C-alphas).
/// Not swept here.
const TORSION_ENM_CUTOFF_ANGSTROM: f32 = 12.0;
/// Eigenvalues at or below this count as the trivial (zero-frequency) modes. Chosen as
/// a round float-noise floor for a matrix of this size; not swept.
const EIGENVALUE_EPS: f64 = 1e-8;
/// Smoothing of the sampled torsion field: passes and the weight each neighbour gets.
/// Keeps a single mode from putting all of its displacement on one torsion. Arbitrary,
/// never swept.
const TORSION_SMOOTHING_PASSES: usize = 2;
const TORSION_SMOOTHING_NEIGHBOR_WEIGHT: f32 = 0.2;
/// Spring constants of the two coupling terms in the torsion Hessian: sequence
/// neighbours (bonded, stiffer) against spatial contacts. The 1.5:1 ratio is arbitrary
/// and was never swept; only their ratio matters, since the mode shapes are scale free.
const SEQUENCE_COUPLING_WEIGHT: f64 = 1.5;
const SPATIAL_COUPLING_WEIGHT: f64 = 1.0;
/// Added to the diagonal so the Hessian of a structure with an isolated residue is
/// still positive definite. Small enough not to move the low-frequency modes; arbitrary.
const DIAGONAL_REGULARIZATION: f64 = 1e-5;
/// Local phi/psi coupling of the same residue: a self term and the off-diagonal that
/// makes the pair prefer to rotate in opposite directions (a crankshaft, which moves
/// the downstream chain least). Both arbitrary, never swept.
const LOCAL_TORSION_SELF_WEIGHT: f64 = 0.25;
const LOCAL_PHI_PSI_COUPLING: f64 = 0.10;
/// A mode counts as a uniform (trivial) torsion mode when the phi and psi fields
/// deviate from their own means by less than this. Eigenvectors are unit norm, so it
/// is a fraction of the mode's length; anything this flat carries no deformation.
const UNIFORM_MODE_RESIDUAL: f64 = 1e-3;
/// Torsion RMS, in radians, at which the response of the backbone to the sampled mode
/// is probed before scaling it to the requested displacement. A quarter of a degree is
/// small enough for the response to be linear in the amplitude (which is what makes a
/// single rescale land on the target) and large enough to stay clear of float noise.
const TORSION_PROBE_RAD: f32 = 0.005;

/// Ensemble shape used by `--enm-sample`: how many conformers to sample and how many
/// low-frequency modes to draw them from. The benchmark on this branch varied the
/// conformer count (5 against 10) and never varied the mode count.
pub const ENSEMBLE_CONFORMERS: usize = 5;
pub const ENSEMBLE_TORSION_MODES: usize = 3;

type TorsionField = Vec<[f32; 2]>;

/// SplitMix64. Sampling has to be reproducible from the seed alone, and this is short
/// enough to audit here rather than resting on a dependency's stability across
/// versions. Constants from Steele et al., "Fast splittable pseudorandom number
/// generators" (2014).
struct SplitMix64(u64);

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }

    /// Uniform in [-1, 1), from the top 53 bits.
    fn next_signed_unit(&mut self) -> f64 {
        let unit = (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
        unit * 2.0 - 1.0
    }
}

/// Half-open residue ranges of the chains, in the order the residues are stored.
///
/// Everything downstream of a torsion is rotated with it, and "downstream" stops at
/// the end of the chain: two chains share no covalent bond, so a phi in one cannot
/// move the other. Without this the dominant term of a multi-chain query's
/// displacement is one chain sliding rigidly away from the rest.
fn chain_segments(chain_per_residue: &[u8]) -> Vec<(usize, usize)> {
    let mut segments = Vec::new();
    if chain_per_residue.is_empty() {
        return segments;
    }
    let mut start = 0usize;
    for i in 1..chain_per_residue.len() {
        if chain_per_residue[i] != chain_per_residue[i - 1] {
            segments.push((start, i));
            start = i;
        }
    }
    segments.push((start, chain_per_residue.len()));
    segments
}

fn build_backbone_arrays(
    structure: &CompactStructure,
) -> Result<(Vec<Coordinate>, Vec<Coordinate>, Vec<Coordinate>), String> {
    let mut n_vec = Vec::with_capacity(structure.num_residues);
    let mut ca_vec = Vec::with_capacity(structure.num_residues);
    let mut c_vec = Vec::with_capacity(structure.num_residues);

    for i in 0..structure.num_residues {
        let ca = structure
            .ca_vector
            .get_coord(i)
            .ok_or_else(|| format!("Missing CA coordinate at residue index {}", i))?;
        let n = structure.n_vector.get_coord(i).unwrap_or(ca);
        let c = structure.c_vector.get_coord(i).unwrap_or(ca);
        n_vec.push(n);
        ca_vec.push(ca);
        c_vec.push(c);
    }

    Ok((n_vec, ca_vec, c_vec))
}

fn add_laplacian_coupling(h: &mut DMatrix<f64>, i: usize, j: usize, weight: f64) {
    for torsion_idx in 0..2 {
        let a = 2 * i + torsion_idx;
        let b = 2 * j + torsion_idx;
        h[(a, a)] += weight;
        h[(b, b)] += weight;
        h[(a, b)] -= weight;
        h[(b, a)] -= weight;
    }
}

fn build_torsion_hessian(ca_vec: &[Coordinate], segments: &[(usize, usize)]) -> DMatrix<f64> {
    let n = ca_vec.len();
    let mut h = DMatrix::<f64>::zeros(2 * n, 2 * n);

    // Spatial ENM couplings in torsion space. These are non-bonded contacts, so they
    // are wanted between chains as much as within one.
    let cutoff2 = TORSION_ENM_CUTOFF_ANGSTROM * TORSION_ENM_CUTOFF_ANGSTROM;
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = ca_vec[j].x - ca_vec[i].x;
            let dy = ca_vec[j].y - ca_vec[i].y;
            let dz = ca_vec[j].z - ca_vec[i].z;
            let dist2 = dx * dx + dy * dy + dz * dz;
            if dist2 <= cutoff2 {
                add_laplacian_coupling(&mut h, i, j, SPATIAL_COUPLING_WEIGHT);
            }
        }
    }

    // Sequence-neighbor couplings (chain smoothness). Bonded, so they stop at the end
    // of each chain.
    for (start, end) in segments {
        for i in *start..end.saturating_sub(1) {
            add_laplacian_coupling(&mut h, i, i + 1, SEQUENCE_COUPLING_WEIGHT);
        }
    }

    // Mild phi/psi local coupling + regularization for numerical stability
    for i in 0..n {
        let phi = 2 * i;
        let psi = phi + 1;
        h[(phi, phi)] += DIAGONAL_REGULARIZATION + LOCAL_TORSION_SELF_WEIGHT;
        h[(psi, psi)] += DIAGONAL_REGULARIZATION + LOCAL_TORSION_SELF_WEIGHT;
        h[(phi, psi)] -= LOCAL_PHI_PSI_COUPLING;
        h[(psi, phi)] -= LOCAL_PHI_PSI_COUPLING;
    }

    h
}

/// True for a mode that only shifts every phi (or every psi) by the same amount.
///
/// Such a mode deforms nothing: it is the torsion-space equivalent of a rigid motion,
/// and the centering step in the sampler removes it exactly. It has to be recognised by
/// its shape rather than by a near-zero eigenvalue, because the local phi/psi
/// self-coupling in the Hessian lifts these modes to an eigenvalue of order
/// `LOCAL_TORSION_SELF_WEIGHT`. Left in, they crowd out the genuine low-frequency modes
/// and every conformer of the ensemble comes out as the same displacement.
fn is_uniform_torsion_mode(mode: &DVector<f64>) -> bool {
    let n_res = mode.len() / 2;
    if n_res == 0 {
        return true;
    }
    let mut residual = 0.0f64;
    for torsion_idx in 0..2 {
        let mean = (0..n_res).map(|i| mode[2 * i + torsion_idx]).sum::<f64>() / n_res as f64;
        residual += (0..n_res)
            .map(|i| (mode[2 * i + torsion_idx] - mean).powi(2))
            .sum::<f64>();
    }
    residual.sqrt() < UNIFORM_MODE_RESIDUAL
}

fn select_nontrivial_modes(hessian: DMatrix<f64>, mode_count: usize) -> Vec<DVector<f64>> {
    if mode_count == 0 {
        return Vec::new();
    }

    let eig = SymmetricEigen::new(hessian);
    let mut ranked: Vec<(f64, usize)> = eig
        .eigenvalues
        .iter()
        .copied()
        .enumerate()
        .map(|(idx, value)| (value, idx))
        .collect();
    ranked.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    ranked
        .into_iter()
        .filter(|(val, _)| *val > EIGENVALUE_EPS)
        .map(|(_, idx)| eig.eigenvectors.column(idx).into_owned())
        .filter(|mode| !is_uniform_torsion_mode(mode))
        .take(mode_count)
        .collect()
}

/// Sample one conformer's torsion displacement from the normal modes.
///
/// `seed` makes the ensemble reproducible: the same query, conformer count and mode
/// count must give the same hashes on every run, otherwise a search is not repeatable
/// and neither are the benchmarks measuring it. Everything from here to the finished
/// conformer is therefore serial - a rayon `sum` over floats reorders the additions
/// with the thread count, which was enough to move a conformer by 2e-5 A.
fn sample_torsion_displacement(
    modes: &[DVector<f64>], n_residues: usize, seed: u64,
) -> TorsionField {
    let mut disp = vec![[0.0f32; 2]; n_residues];
    if modes.is_empty() {
        return disp;
    }

    let mut rng = SplitMix64::new(seed);
    let coeffs: Vec<f64> = (0..modes.len()).map(|_| rng.next_signed_unit()).collect();

    for (i, d) in disp.iter_mut().enumerate() {
        for torsion_idx in 0..2 {
            let idx = 2 * i + torsion_idx;
            let mut value = 0.0f64;
            for (mode, c) in modes.iter().zip(coeffs.iter()) {
                value += *c * mode[idx];
            }
            d[torsion_idx] = value as f32;
        }
    }

    disp
}

fn torsion_rmsd(disp: &TorsionField) -> f32 {
    if disp.is_empty() {
        return 0.0;
    }

    let sum: f32 = disp.iter().map(|d| d[0] * d[0] + d[1] * d[1]).sum();
    (sum / ((disp.len() * 2) as f32)).sqrt()
}

fn scale_torsion_displacement(disp: &mut TorsionField, scale: f32) {
    for d in disp.iter_mut() {
        d[0] *= scale;
        d[1] *= scale;
    }
}

/// Rescale the field to a given torsion RMS, in radians. Returns false when the field
/// is all but zero and cannot be scaled to anything.
fn set_torsion_rmsd(disp: &mut TorsionField, target_rad: f32) -> bool {
    let current = torsion_rmsd(disp);
    if current <= f32::EPSILON {
        return false;
    }
    scale_torsion_displacement(disp, target_rad / current);
    true
}

fn center_and_smooth_torsion_displacement(disp: &mut TorsionField, segments: &[(usize, usize)]) {
    if disp.is_empty() {
        return;
    }

    let n = disp.len() as f32;
    let mean_phi = disp.iter().map(|d| d[0]).sum::<f32>() / n;
    let mean_psi = disp.iter().map(|d| d[1]).sum::<f32>() / n;

    for d in disp.iter_mut() {
        d[0] -= mean_phi;
        d[1] -= mean_psi;
    }

    if TORSION_SMOOTHING_PASSES == 0 {
        return;
    }

    // Smoothing runs along the chain, so it stops where the chain does: averaging the
    // last torsion of one chain with the first of the next has no physical meaning.
    let mut scratch = disp.clone();
    let self_weight = 1.0 - (2.0 * TORSION_SMOOTHING_NEIGHBOR_WEIGHT);
    for _ in 0..TORSION_SMOOTHING_PASSES {
        for (start, end) in segments {
            if end - start < 3 {
                continue;
            }
            for i in (start + 1)..(end - 1) {
                for torsion_idx in 0..2 {
                    scratch[i][torsion_idx] = disp[i][torsion_idx] * self_weight
                        + disp[i - 1][torsion_idx] * TORSION_SMOOTHING_NEIGHBOR_WEIGHT
                        + disp[i + 1][torsion_idx] * TORSION_SMOOTHING_NEIGHBOR_WEIGHT;
                }
            }
            scratch[*start] = disp[*start];
            scratch[end - 1] = disp[end - 1];
        }
        std::mem::swap(disp, &mut scratch);
    }
}

fn rotate_point_around_axis(
    point: &Coordinate,
    pivot: &Coordinate,
    axis: &Coordinate,
    angle: f32,
) -> Coordinate {
    let axis_norm = axis.norm();
    if axis_norm <= f32::EPSILON || angle.abs() <= f32::EPSILON {
        return *point;
    }

    let k = axis.scale(1.0 / axis_norm);
    let v = point.sub(pivot);
    let cos_t = angle.cos();
    let sin_t = angle.sin();

    // Rodrigues' rotation formula
    let term1 = v.scale(cos_t);
    let term2 = k.cross(&v).scale(sin_t);
    let term3 = k.scale(k.dot(&v) * (1.0 - cos_t));
    pivot.add(&term1.add(&term2).add(&term3))
}

/// Rotate `coords[start..end]` about the axis from `pivot_a` to `pivot_b`.
fn apply_axis_rotation(
    coords: &mut [Coordinate],
    start: usize,
    end: usize,
    pivot_a: &Coordinate,
    pivot_b: &Coordinate,
    angle: f32,
) {
    if start >= end || end > coords.len() || angle.abs() <= f32::EPSILON {
        return;
    }

    let axis = pivot_b.sub(pivot_a);
    for coord in coords[start..end].iter_mut() {
        *coord = rotate_point_around_axis(coord, pivot_a, &axis, angle);
    }
}

fn build_backbone_atom_coords(
    n_vec: &[Coordinate],
    ca_vec: &[Coordinate],
    c_vec: &[Coordinate],
) -> Vec<Coordinate> {
    let mut atoms = Vec::with_capacity(n_vec.len() * 3);
    for i in 0..n_vec.len() {
        atoms.push(n_vec[i]);
        atoms.push(ca_vec[i]);
        atoms.push(c_vec[i]);
    }
    atoms
}

/// Apply a torsion displacement to a backbone atom array laid out as
/// `[N, CA, C]` per residue, one chain at a time.
fn displace_backbone(
    base_atoms: &[Coordinate], segments: &[(usize, usize)], disp: &TorsionField,
) -> Vec<Coordinate> {
    let mut atoms = base_atoms.to_vec();
    for (start, end) in segments {
        let (start, end) = (*start, *end);
        for i in start..end {
            // phi rotation around N-CA; rotate C(i) and everything after it in this chain
            if i > start {
                let n_i = atoms[3 * i];
                let ca_i = atoms[3 * i + 1];
                apply_axis_rotation(&mut atoms, 3 * i + 2, 3 * end, &n_i, &ca_i, disp[i][0]);
            }

            // psi rotation around CA-C; rotate the next residue and everything after it
            if i + 1 < end {
                let ca_i = atoms[3 * i + 1];
                let c_i = atoms[3 * i + 2];
                apply_axis_rotation(&mut atoms, 3 * (i + 1), 3 * end, &ca_i, &c_i, disp[i][1]);
            }
        }
    }
    atoms
}

/// Backbone displacement RMSD between two atom arrays in the same frame, in Angstroms.
///
/// No superposition: the rotations keep the first residue of every chain fixed, so the
/// two conformations already share a frame and this is the distance the atoms actually
/// moved. It is the quantity `--nma-rmsd` names.
fn backbone_displacement_rmsd(base_atoms: &[Coordinate], atoms: &[Coordinate]) -> f32 {
    if base_atoms.is_empty() || base_atoms.len() != atoms.len() {
        return 0.0;
    }
    let mut sum = 0.0f64;
    for (a, b) in base_atoms.iter().zip(atoms.iter()) {
        sum += a.calc_distance(b).powi(2) as f64;
    }
    (sum / base_atoms.len() as f64).sqrt() as f32
}

/// Write a displaced backbone back into a copy of the structure, C-beta included.
fn structure_with_backbone(structure: &CompactStructure, atoms: &[Coordinate]) -> CompactStructure {
    let mut sampled = structure.clone();
    for i in 0..structure.num_residues {
        let n = atoms[3 * i];
        let ca = atoms[3 * i + 1];
        let c = atoms[3 * i + 2];

        sampled.n_vector.x[i] = Some(n.x);
        sampled.n_vector.y[i] = Some(n.y);
        sampled.n_vector.z[i] = Some(n.z);

        sampled.ca_vector.x[i] = Some(ca.x);
        sampled.ca_vector.y[i] = Some(ca.y);
        sampled.ca_vector.z[i] = Some(ca.z);

        sampled.c_vector.x[i] = Some(c.x);
        sampled.c_vector.y[i] = Some(c.y);
        sampled.c_vector.z[i] = Some(c.z);

        let cb = approx_cb(&ca, &n, &c);
        sampled.cb_vector.x[i] = Some(cb.x);
        sampled.cb_vector.y[i] = Some(cb.y);
        sampled.cb_vector.z[i] = Some(cb.z);
    }
    sampled
}

/// Generate an ensemble consisting of the original structure plus torsion-ENM sampled
/// conformations, each displaced by `target_backbone_rmsd` Angstroms of backbone RMSD.
///
/// The requested displacement is delivered rather than approximated: the response of a
/// backbone to a torsion perturbation depends on the length of the chain downstream of
/// it, so there is no fixed radians-per-Angstrom. Each conformer is therefore probed at
/// a small amplitude and scaled once by the ratio to the target, which lands inside a
/// few percent while the response is linear.
pub fn generate_ensemble(
    query_structure: &CompactStructure,
    num_confs: usize,
    target_backbone_rmsd: f32,
    nma_modes: usize,
) -> Result<Vec<CompactStructure>, String> {
    let mut ensemble = Vec::with_capacity(num_confs.saturating_add(1));
    ensemble.push(query_structure.clone());

    let target_backbone_rmsd = target_backbone_rmsd.abs();
    if num_confs == 0 || nma_modes == 0 || query_structure.num_residues < 3
        || target_backbone_rmsd <= 0.0 {
        return Ok(ensemble);
    }

    let (n_vec, ca_vec, c_vec) = build_backbone_arrays(query_structure)?;
    let segments = chain_segments(&query_structure.chain_per_residue[..query_structure.num_residues]);
    let hessian = build_torsion_hessian(&ca_vec, &segments);
    let modes = select_nontrivial_modes(hessian, nma_modes);

    if modes.is_empty() {
        return Ok(ensemble);
    }
    let base_atoms = build_backbone_atom_coords(&n_vec, &ca_vec, &c_vec);

    let sampled: Vec<CompactStructure> = (0..num_confs)
        .into_par_iter()
        .map(|conf_index| {
            // Deterministic per conformer, so the ensemble does not depend on how
            // rayon schedules the work.
            let mut disp =
                sample_torsion_displacement(&modes, query_structure.num_residues, conf_index as u64);
            center_and_smooth_torsion_displacement(&mut disp, &segments);
            if !set_torsion_rmsd(&mut disp, TORSION_PROBE_RAD) {
                return query_structure.clone();
            }
            let probe = displace_backbone(&base_atoms, &segments, &disp);
            let probe_rmsd = backbone_displacement_rmsd(&base_atoms, &probe);
            if probe_rmsd <= f32::EPSILON {
                return query_structure.clone();
            }
            scale_torsion_displacement(&mut disp, target_backbone_rmsd / probe_rmsd);
            let atoms = displace_backbone(&base_atoms, &segments, &disp);
            structure_with_backbone(query_structure, &atoms)
        })
        .collect();

    ensemble.extend(sampled);
    Ok(ensemble)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::controller::io::read_compact_structure;

    /// Both queries shipped with the tool. Multi-chain on purpose: 1G2F has 2 chains
    /// and 4CHA has 6, and treating them as one continuous chain is exactly the bug
    /// these tests exist to keep out.
    const EXAMPLE_QUERIES: [&str; 2] = ["query/1G2F.pdb", "query/4CHA.pdb"];

    fn load(path: &str) -> CompactStructure {
        read_compact_structure(path).expect("failed to read structure").0
    }

    /// Backbone (N, CA, C) displacement RMSD between two conformations of the same
    /// structure, in the frame they are stored in.
    fn backbone_rmsd(a: &CompactStructure, b: &CompactStructure) -> f32 {
        let mut sum = 0.0f64;
        let mut count = 0usize;
        for i in 0..a.num_residues {
            for (va, vb) in [
                (&a.n_vector, &b.n_vector), (&a.ca_vector, &b.ca_vector), (&a.c_vector, &b.c_vector),
            ] {
                if let (Some(pa), Some(pb)) = (va.get_coord(i), vb.get_coord(i)) {
                    sum += pa.calc_distance(&pb).powi(2) as f64;
                    count += 1;
                }
            }
        }
        if count == 0 { 0.0 } else { (sum / count as f64).sqrt() as f32 }
    }

    #[test]
    fn rotating_one_chain_leaves_the_others_where_they_were() {
        let structure = load("query/1G2F.pdb");
        let (n_vec, ca_vec, c_vec) = build_backbone_arrays(&structure).unwrap();
        let base = build_backbone_atom_coords(&n_vec, &ca_vec, &c_vec);
        let segments = chain_segments(&structure.chain_per_residue[..structure.num_residues]);
        assert!(segments.len() > 1, "the test structure must have more than one chain");
        let boundary = segments[1].0;

        // A psi at the very start of the first chain: everything after it inside that
        // chain swings, and nothing in any later chain may move at all.
        let mut disp = vec![[0.0f32; 2]; structure.num_residues];
        disp[0][1] = 0.3;
        let moved = displace_backbone(&base, &segments, &disp);

        for atom in (3 * boundary)..base.len() {
            let drift = base[atom].calc_distance(&moved[atom]);
            assert!(drift < 1e-4,
                "atom {} of a later chain moved {} A for a torsion in the first chain",
                atom, drift);
        }
        let last_of_first_chain = base[3 * boundary - 1]
            .calc_distance(&moved[3 * boundary - 1]);
        assert!(last_of_first_chain > 1.0,
            "the rotation did not propagate along its own chain ({} A)", last_of_first_chain);
    }

    #[test]
    fn torsion_rotation_preserves_bonded_geometry() {
        // Rotating about the N-CA and CA-C axes leaves both pivot atoms on the axis,
        // so every bonded distance is invariant by construction, at any amplitude.
        // This is the invariant that makes a plausibility check on bond lengths
        // pointless, and it has to hold for a displacement far larger than any the
        // sampler asks for.
        for path in EXAMPLE_QUERIES {
            let structure = load(path);
            let ensemble = generate_ensemble(&structure, 2, 20.0, 3).expect("sampling failed");
            for conformer in ensemble.iter().skip(1) {
                for i in 0..structure.num_residues {
                    let before_n_ca = structure.n_vector.get_coord(i).unwrap()
                        .calc_distance(&structure.ca_vector.get_coord(i).unwrap());
                    let after_n_ca = conformer.n_vector.get_coord(i).unwrap()
                        .calc_distance(&conformer.ca_vector.get_coord(i).unwrap());
                    assert!((before_n_ca - after_n_ca).abs() < 1e-2,
                        "{} residue {}: N-CA {} -> {}", path, i, before_n_ca, after_n_ca);
                    let before_ca_c = structure.ca_vector.get_coord(i).unwrap()
                        .calc_distance(&structure.c_vector.get_coord(i).unwrap());
                    let after_ca_c = conformer.ca_vector.get_coord(i).unwrap()
                        .calc_distance(&conformer.c_vector.get_coord(i).unwrap());
                    assert!((before_ca_c - after_ca_c).abs() < 1e-2,
                        "{} residue {}: CA-C {} -> {}", path, i, before_ca_c, after_ca_c);
                    // Peptide bond and adjacent CA-CA, inside a chain only
                    if i + 1 < structure.num_residues
                        && structure.chain_per_residue[i] == structure.chain_per_residue[i + 1] {
                        let before_c_n = structure.c_vector.get_coord(i).unwrap()
                            .calc_distance(&structure.n_vector.get_coord(i + 1).unwrap());
                        let after_c_n = conformer.c_vector.get_coord(i).unwrap()
                            .calc_distance(&conformer.n_vector.get_coord(i + 1).unwrap());
                        assert!((before_c_n - after_c_n).abs() < 1e-2,
                            "{} residue {}: C-N {} -> {}", path, i, before_c_n, after_c_n);
                        let before_ca_ca = structure.ca_vector.get_coord(i).unwrap()
                            .calc_distance(&structure.ca_vector.get_coord(i + 1).unwrap());
                        let after_ca_ca = conformer.ca_vector.get_coord(i).unwrap()
                            .calc_distance(&conformer.ca_vector.get_coord(i + 1).unwrap());
                        assert!((before_ca_ca - after_ca_ca).abs() < 1e-2,
                            "{} residue {}: CA-CA {} -> {}", path, i, before_ca_ca, after_ca_ca);
                    }
                }
            }
        }
    }

    #[test]
    fn achieved_backbone_rmsd_matches_the_request() {
        // --nma-rmsd is documented in Angstroms, so the conformers have to land on
        // the value asked for rather than on some length-dependent multiple of it.
        for path in EXAMPLE_QUERIES {
            let structure = load(path);
            for requested in [0.25f32, 0.5, 1.0] {
                let ensemble = generate_ensemble(&structure, 3, requested, 3).expect("sampling failed");
                assert_eq!(ensemble.len(), 4);
                for conformer in ensemble.iter().skip(1) {
                    let achieved = backbone_rmsd(&structure, conformer);
                    let error = (achieved - requested).abs() / requested;
                    assert!(error < 0.15,
                        "{}: requested {} A, achieved {:.4} A ({:.0}% off)",
                        path, requested, achieved, error * 100.0);
                }
            }
        }
    }

    #[test]
    fn conformers_deform_the_query_without_being_the_query() {
        let structure = load("query/1G2F.pdb");
        let ensemble = generate_ensemble(&structure, 3, 0.5, 3).expect("sampling failed");
        for conformer in ensemble.iter().skip(1) {
            assert!(backbone_rmsd(&structure, conformer) > 1e-3, "conformer is a copy of the query");
        }
        // Distinct conformers, not the same displacement three times
        for i in 1..ensemble.len() {
            for j in (i + 1)..ensemble.len() {
                assert!(backbone_rmsd(&ensemble[i], &ensemble[j]) > 1e-3,
                    "conformers {} and {} are identical", i, j);
            }
        }
    }

    #[test]
    fn ensemble_is_reproducible_and_thread_count_independent() {
        let structure = load("query/1G2F.pdb");
        let first = generate_ensemble(&structure, 4, 0.5, 3).expect("sampling failed");
        let second = generate_ensemble(&structure, 4, 0.5, 3).expect("sampling failed");
        let single_threaded = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap()
            .install(|| generate_ensemble(&structure, 4, 0.5, 3).expect("sampling failed"));
        for (i, conformer) in first.iter().enumerate() {
            assert_eq!(backbone_rmsd(conformer, &second[i]), 0.0, "conformer {} not reproducible", i);
            assert_eq!(backbone_rmsd(conformer, &single_threaded[i]), 0.0,
                "conformer {} depends on the thread count", i);
        }
    }

    #[test]
    fn degenerate_requests_return_the_query_alone() {
        let structure = load("query/1G2F.pdb");
        assert_eq!(generate_ensemble(&structure, 0, 0.5, 3).unwrap().len(), 1);
        assert_eq!(generate_ensemble(&structure, 5, 0.5, 0).unwrap().len(), 1);
        // A structure too short to have an interior torsion
        let mut tiny = structure.clone();
        tiny.num_residues = 2;
        assert_eq!(generate_ensemble(&tiny, 5, 0.5, 3).unwrap().len(), 1);
    }
}
