//! Torsion-angle ENM/NMA conformational sampling for non-rigid query search.
//!
//! This module builds a residue-level elastic network in torsion space (phi/psi),
//! samples low-frequency torsional normal modes, and applies the perturbations by
//! rotating downstream backbone atoms around peptide bond axes.

use nalgebra::{DMatrix, DVector, SymmetricEigen};
use rand::Rng;
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;

use crate::structure::coordinate::{approx_cb, calc_torsion_radian, Coordinate};
use crate::structure::core::CompactStructure;

const TORSION_ENM_CUTOFF_ANGSTROM: f32 = 12.0;
const EIGENVALUE_EPS: f64 = 1e-8;
const TORSION_RADIANS_PER_ANGSTROM: f32 = 0.35;
const MAX_TORSION_STEP_RAD: f32 = 0.65;
const SEQUENCE_COUPLING_WEIGHT: f64 = 1.5;
const SPATIAL_COUPLING_WEIGHT: f64 = 1.0;
const DIAGONAL_REGULARIZATION: f64 = 1e-5;

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

fn extract_backbone_torsions(
    n_vec: &[Coordinate],
    ca_vec: &[Coordinate],
    c_vec: &[Coordinate],
) -> TorsionField {
    let n_res = ca_vec.len();
    let mut torsions = vec![[0.0f32, 0.0f32]; n_res];

    for i in 0..n_res {
        // phi(i) = C(i-1)-N(i)-CA(i)-C(i)
        if i > 0 {
            torsions[i][0] = calc_torsion_radian(&c_vec[i - 1], &n_vec[i], &ca_vec[i], &c_vec[i]);
        }
        // psi(i) = N(i)-CA(i)-C(i)-N(i+1)
        if i + 1 < n_res {
            torsions[i][1] = calc_torsion_radian(&n_vec[i], &ca_vec[i], &c_vec[i], &n_vec[i + 1]);
        }
    }

    torsions
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
        h[(phi, phi)] += DIAGONAL_REGULARIZATION + 0.25;
        h[(psi, psi)] += DIAGONAL_REGULARIZATION + 0.25;
        h[(phi, psi)] -= 0.10;
        h[(psi, phi)] -= 0.10;
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

fn sample_torsion_displacement(modes: &[DVector<f64>], n_residues: usize) -> TorsionField {
    let mut disp = vec![[0.0f32; 2]; n_residues];
    if modes.is_empty() {
        return disp;
    }

    let mut rng = rand::thread_rng();
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

    let (n_vec, ca_vec, c_vec) = build_backbone_arrays(query_structure)?;
    let _base_torsions = extract_backbone_torsions(&n_vec, &ca_vec, &c_vec);
    let hessian = build_torsion_hessian(&ca_vec);
    let modes = select_nontrivial_modes(hessian, nma_modes);

    if modes.is_empty() {
        return Ok(ensemble);
    }

    let sampled: Vec<CompactStructure> = (0..num_confs)
        .into_par_iter()
        .map(|_| {
            let mut best: Option<(f32, CompactStructure)> = None;
            for attempt in 0..MAX_CONFORMER_SAMPLING_ATTEMPTS {
                let mut disp = sample_torsion_displacement(&modes, query_structure.num_residues);
                let attempt_scale = CONFORMER_RETRY_SCALE_FACTOR.powi(attempt as i32);
                scale_torsion_displacement(&mut disp, target_rmsd * attempt_scale);
                cap_torsion_step(&mut disp, MAX_TORSION_STEP_RAD);
                let conformer =
                    apply_torsion_displacement(query_structure, &n_vec, &ca_vec, &c_vec, &disp);
                if is_backbone_plausible(&conformer) {
                    return conformer;
                }

                let eff = torsion_rmsd(&disp);
                if best.as_ref().map_or(true, |(curr, _)| eff > *curr) {
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

fn atom_line(
    serial: usize,
    atom_name: &str,
    res_name: &[u8; 3],
    chain: u8,
    res_seq: u64,
    coord: &Coordinate,
) -> String {
    let residue = std::str::from_utf8(res_name).unwrap_or("UNK");
    let element = atom_name.trim().chars().next().unwrap_or('C');
    format!(
        "ATOM  {:>5} {:<4} {:>3} {:>1}{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00           {:>2}",
        serial, atom_name, residue, chain as char, res_seq, coord.x, coord.y, coord.z, element
    )
}

pub fn write_conformer_as_pdb(structure: &CompactStructure, path: &str) -> Result<(), String> {
    let mut file = File::create(path).map_err(|e| format!("Failed to create {}: {}", path, e))?;
    let mut serial = 1usize;
    for i in 0..structure.num_residues {
        let chain = structure.chain_per_residue[i];
        let res_seq = structure.residue_serial[i];
        let res_name = &structure.residue_name[i];
        if let Some(n) = structure.n_vector.get_coord(i) {
            writeln!(
                file,
                "{}",
                atom_line(serial, "N", res_name, chain, res_seq, &n)
            )
            .map_err(|e| format!("Failed to write {}: {}", path, e))?;
            serial += 1;
        }
        if let Some(ca) = structure.ca_vector.get_coord(i) {
            writeln!(
                file,
                "{}",
                atom_line(serial, "CA", res_name, chain, res_seq, &ca)
            )
            .map_err(|e| format!("Failed to write {}: {}", path, e))?;
            serial += 1;
        }
        if let Some(c) = structure.c_vector.get_coord(i) {
            writeln!(
                file,
                "{}",
                atom_line(serial, "C", res_name, chain, res_seq, &c)
            )
            .map_err(|e| format!("Failed to write {}: {}", path, e))?;
            serial += 1;
        }
        if let Some(cb) = structure.cb_vector.get_coord(i) {
            writeln!(
                file,
                "{}",
                atom_line(serial, "CB", res_name, chain, res_seq, &cb)
            )
            .map_err(|e| format!("Failed to write {}: {}", path, e))?;
            serial += 1;
        }
    }
    writeln!(file, "END").map_err(|e| format!("Failed to finalize {}: {}", path, e))?;
    Ok(())
}
