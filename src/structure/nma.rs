//! Torsion-angle ENM/NMA conformational sampling for non-rigid query search.
//!
//! This module builds a residue-level elastic network in torsion space (phi/psi),
//! samples low-frequency torsional normal modes, and applies the perturbations by
//! rotating downstream backbone atoms around peptide bond axes.

use nalgebra::{DMatrix, DVector, SymmetricEigen};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
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

/// Ensemble shape used by `--enm-sample`: how many conformers to sample and how many
/// low-frequency modes to draw them from. The benchmark on this branch varied the
/// conformer count (5 against 10) and never varied the mode count.
pub const ENSEMBLE_CONFORMERS: usize = 5;
pub const ENSEMBLE_TORSION_MODES: usize = 3;

const TORSION_RADIANS_PER_ANGSTROM: f32 = 0.06;
const MAX_TORSION_STEP_RAD: f32 = 0.12;
const MAX_TARGET_RMSD_FOR_WIGGLE: f32 = 0.75;
const MIN_TARGET_RMSD_FOR_WIGGLE: f32 = 0.05;

// Plausibility limits of a backbone, in Angstroms: bonded distances and the widest
// adjacent C-alpha separation a peptide can show. Physical ranges, deliberately loose.
const MIN_N_CA_BOND: f32 = 0.7;
const MAX_N_CA_BOND: f32 = 2.6;
const MIN_CA_C_BOND: f32 = 0.7;
const MAX_CA_C_BOND: f32 = 2.6;
const MIN_C_N_PEPTIDE_BOND: f32 = 0.7;
const MAX_C_N_PEPTIDE_BOND: f32 = 2.6;
const MAX_ADJACENT_CA_DISTANCE: f32 = 6.5;
const MAX_CONFORMER_SAMPLING_ATTEMPTS: usize = 8;
const CONFORMER_RETRY_SCALE_FACTOR: f32 = 0.85;

type TorsionField = Vec<[f32; 2]>;

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

fn build_torsion_hessian(ca_vec: &[Coordinate]) -> DMatrix<f64> {
    let n = ca_vec.len();
    let mut h = DMatrix::<f64>::zeros(2 * n, 2 * n);

    // Spatial ENM couplings in torsion space
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

    // Sequence-neighbor couplings (chain smoothness)
    for i in 0..n.saturating_sub(1) {
        add_laplacian_coupling(&mut h, i, i + 1, SEQUENCE_COUPLING_WEIGHT);
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
        .take(mode_count)
        .map(|(_, idx)| eig.eigenvectors.column(idx).into_owned())
        .collect()
}

/// Sample one conformer's torsion displacement from the normal modes.
///
/// `seed` makes the ensemble reproducible: the same query, conformer count and mode
/// count must give the same hashes on every run, otherwise a search is not repeatable
/// and neither are the benchmarks measuring it.
fn sample_torsion_displacement(
    modes: &[DVector<f64>], n_residues: usize, seed: u64,
) -> TorsionField {
    let mut disp = vec![[0.0f32; 2]; n_residues];
    if modes.is_empty() {
        return disp;
    }

    let mut rng = StdRng::seed_from_u64(seed);
    let coeffs: Vec<f64> = (0..modes.len())
        .map(|_| rng.gen_range(-1.0f64..1.0f64))
        .collect();

    disp.par_iter_mut().enumerate().for_each(|(i, d)| {
        for torsion_idx in 0..2 {
            let idx = 2 * i + torsion_idx;
            let mut value = 0.0f64;
            for (mode, c) in modes.iter().zip(coeffs.iter()) {
                value += *c * mode[idx];
            }
            d[torsion_idx] = value as f32;
        }
    });

    disp
}

fn torsion_rmsd(disp: &TorsionField) -> f32 {
    if disp.is_empty() {
        return 0.0;
    }

    let sum: f32 = disp.par_iter().map(|d| d[0] * d[0] + d[1] * d[1]).sum();
    (sum / ((disp.len() * 2) as f32)).sqrt()
}

fn scale_torsion_displacement(disp: &mut TorsionField, target_rmsd_angstrom: f32) {
    let target_torsion_rmsd = (target_rmsd_angstrom.abs() * TORSION_RADIANS_PER_ANGSTROM).max(0.0);
    let current = torsion_rmsd(disp);
    if current <= f32::EPSILON || target_torsion_rmsd <= 0.0 {
        return;
    }

    let scale = target_torsion_rmsd / current;
    disp.par_iter_mut().for_each(|d| {
        d[0] *= scale;
        d[1] *= scale;
    });
}

fn cap_torsion_step(disp: &mut TorsionField, cap: f32) {
    if cap <= 0.0 {
        return;
    }

    disp.par_iter_mut().for_each(|d| {
        d[0] = d[0].clamp(-cap, cap);
        d[1] = d[1].clamp(-cap, cap);
    });
}

fn center_and_smooth_torsion_displacement(disp: &mut TorsionField) {
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

    if disp.len() < 3 || TORSION_SMOOTHING_PASSES == 0 {
        return;
    }

    let mut scratch = disp.clone();
    let self_weight = 1.0 - (2.0 * TORSION_SMOOTHING_NEIGHBOR_WEIGHT);
    for _ in 0..TORSION_SMOOTHING_PASSES {
        for i in 1..(disp.len() - 1) {
            for torsion_idx in 0..2 {
                scratch[i][torsion_idx] = disp[i][torsion_idx] * self_weight
                    + disp[i - 1][torsion_idx] * TORSION_SMOOTHING_NEIGHBOR_WEIGHT
                    + disp[i + 1][torsion_idx] * TORSION_SMOOTHING_NEIGHBOR_WEIGHT;
            }
        }
        scratch[0] = disp[0];
        scratch[disp.len() - 1] = disp[disp.len() - 1];
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

fn apply_axis_rotation(
    coords: &mut [Coordinate],
    start_idx: usize,
    pivot_a: &Coordinate,
    pivot_b: &Coordinate,
    angle: f32,
) {
    if start_idx >= coords.len() || angle.abs() <= f32::EPSILON {
        return;
    }

    let axis = pivot_b.sub(pivot_a);
    for coord in coords.iter_mut().skip(start_idx) {
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

fn apply_torsion_displacement(
    structure: &CompactStructure,
    n_vec: &[Coordinate],
    ca_vec: &[Coordinate],
    c_vec: &[Coordinate],
    disp: &TorsionField,
) -> CompactStructure {
    let mut sampled = structure.clone();
    let mut atoms = build_backbone_atom_coords(n_vec, ca_vec, c_vec);
    let n_res = structure.num_residues;

    for i in 0..n_res {
        // phi rotation around N-CA; rotate C(i) and all downstream atoms
        if i > 0 {
            let n_i = atoms[3 * i];
            let ca_i = atoms[3 * i + 1];
            let phi_delta = disp[i][0];
            apply_axis_rotation(&mut atoms, 3 * i + 2, &n_i, &ca_i, phi_delta);
        }

        // psi rotation around CA-C; rotate next residue and downstream atoms
        if i + 1 < n_res {
            let ca_i = atoms[3 * i + 1];
            let c_i = atoms[3 * i + 2];
            let psi_delta = disp[i][1];
            apply_axis_rotation(&mut atoms, 3 * (i + 1), &ca_i, &c_i, psi_delta);
        }
    }

    for i in 0..n_res {
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

fn is_backbone_plausible(structure: &CompactStructure) -> bool {
    if structure.num_residues < 2 {
        return true;
    }

    for i in 0..structure.num_residues {
        let (n, ca, c) = (
            structure.n_vector.get_coord(i),
            structure.ca_vector.get_coord(i),
            structure.c_vector.get_coord(i),
        );
        if let (Some(n), Some(ca), Some(c)) = (n, ca, c) {
            let n_ca = n.calc_distance(&ca);
            let ca_c = ca.calc_distance(&c);
            if !(MIN_N_CA_BOND..=MAX_N_CA_BOND).contains(&n_ca) {
                return false;
            }
            if !(MIN_CA_C_BOND..=MAX_CA_C_BOND).contains(&ca_c) {
                return false;
            }
            if i + 1 < structure.num_residues {
                if let Some(next_n) = structure.n_vector.get_coord(i + 1) {
                    let c_n = c.calc_distance(&next_n);
                    if !(MIN_C_N_PEPTIDE_BOND..=MAX_C_N_PEPTIDE_BOND).contains(&c_n) {
                        return false;
                    }
                }
                if let Some(next_ca) = structure.ca_vector.get_coord(i + 1) {
                    let ca_ca = ca.calc_distance(&next_ca);
                    if ca_ca > MAX_ADJACENT_CA_DISTANCE {
                        return false;
                    }
                }
            }
        }
    }
    true
}

/// Generate an ensemble consisting of the original structure plus torsion-ENM sampled conformations.
pub fn generate_ensemble(
    query_structure: &CompactStructure,
    num_confs: usize,
    target_rmsd: f32,
    nma_modes: usize,
) -> Result<Vec<CompactStructure>, String> {
    let mut ensemble = Vec::with_capacity(num_confs.saturating_add(1));
    ensemble.push(query_structure.clone());

    if num_confs == 0 || nma_modes == 0 || query_structure.num_residues < 3 {
        return Ok(ensemble);
    }
    let effective_target_rmsd = target_rmsd
        .abs()
        .clamp(MIN_TARGET_RMSD_FOR_WIGGLE, MAX_TARGET_RMSD_FOR_WIGGLE);

    let (n_vec, ca_vec, c_vec) = build_backbone_arrays(query_structure)?;
    let hessian = build_torsion_hessian(&ca_vec);
    let modes = select_nontrivial_modes(hessian, nma_modes);

    if modes.is_empty() {
        return Ok(ensemble);
    }

    let sampled: Vec<CompactStructure> = (0..num_confs)
        .into_par_iter()
        .map(|conf_index| {
            let mut best: Option<(f32, CompactStructure)> = None;
            for attempt in 0..MAX_CONFORMER_SAMPLING_ATTEMPTS {
                // Deterministic per (conformer, attempt) so the ensemble is reproducible
                // regardless of how rayon schedules the work.
                let seed = (conf_index as u64) << 8 | attempt as u64;
                let mut disp =
                    sample_torsion_displacement(&modes, query_structure.num_residues, seed);
                let attempt_scale = CONFORMER_RETRY_SCALE_FACTOR.powi(attempt as i32);
                scale_torsion_displacement(&mut disp, effective_target_rmsd * attempt_scale);
                center_and_smooth_torsion_displacement(&mut disp);
                cap_torsion_step(&mut disp, MAX_TORSION_STEP_RAD);
                let conformer =
                    apply_torsion_displacement(query_structure, &n_vec, &ca_vec, &c_vec, &disp);
                if is_backbone_plausible(&conformer) {
                    return conformer;
                }

                let eff = torsion_rmsd(&disp);
                if best.as_ref().map_or(true, |(curr, _)| eff < *curr) {
                    best = Some((eff, conformer));
                }
            }
            best.map(|(_, conf)| conf)
                .unwrap_or_else(|| query_structure.clone())
        })
        .collect();

    ensemble.extend(sampled);
    Ok(ensemble)
}
