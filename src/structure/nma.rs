//! Native ANM/NMA conformational sampling for non-rigid query search.
//!
//! This module implements a lightweight ProDy-like ANM pipeline:
//! 1) build a 3N×3N Hessian from Cα coordinates (cutoff ENM),
//! 2) eigendecompose and discard rigid-body modes,
//! 3) sample random linear combinations of low-frequency modes,
//! 4) scale displacements to a target RMSD and apply them to N/CA/CB atoms.

use nalgebra::{DMatrix, DVector, SymmetricEigen};
use rand::Rng;

use crate::structure::coordinate::Coordinate;
use crate::structure::core::CompactStructure;

const ANM_CUTOFF_ANGSTROM: f64 = 15.0;
const EIGENVALUE_EPS: f64 = 1e-8;

type DisplacementField = Vec<[f32; 3]>;

fn extract_ca_coordinates(structure: &CompactStructure) -> Result<Vec<Coordinate>, String> {
    let mut coords = Vec::with_capacity(structure.num_residues);
    for i in 0..structure.num_residues {
        let ca = structure
            .ca_vector
            .get_coord(i)
            .ok_or_else(|| format!("Missing CA coordinate at residue index {}", i))?;
        coords.push(ca);
    }
    Ok(coords)
}

fn build_anm_hessian(ca_coords: &[Coordinate], cutoff: f64) -> DMatrix<f64> {
    let n = ca_coords.len();
    let size = 3 * n;
    let mut h = DMatrix::<f64>::zeros(size, size);
    let cutoff2 = cutoff * cutoff;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = (ca_coords[j].x - ca_coords[i].x) as f64;
            let dy = (ca_coords[j].y - ca_coords[i].y) as f64;
            let dz = (ca_coords[j].z - ca_coords[i].z) as f64;
            let dist2 = dx * dx + dy * dy + dz * dz;
            if dist2 <= f64::EPSILON || dist2 > cutoff2 {
                continue;
            }

            let inv_dist2 = 1.0 / dist2;
            let block = [
                [
                    dx * dx * inv_dist2,
                    dx * dy * inv_dist2,
                    dx * dz * inv_dist2,
                ],
                [
                    dy * dx * inv_dist2,
                    dy * dy * inv_dist2,
                    dy * dz * inv_dist2,
                ],
                [
                    dz * dx * inv_dist2,
                    dz * dy * inv_dist2,
                    dz * dz * inv_dist2,
                ],
            ];

            for a in 0..3 {
                for b in 0..3 {
                    let ii = 3 * i + a;
                    let jj = 3 * j + b;
                    h[(ii, jj)] -= block[a][b];
                    h[(jj, ii)] -= block[a][b];
                    h[(ii, ii)] += block[a][b];
                    h[(jj, jj)] += block[a][b];
                }
            }
        }
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
        .skip(6)
        .filter(|(val, _)| *val > EIGENVALUE_EPS)
        .take(mode_count)
        .map(|(_, idx)| eig.eigenvectors.column(idx).into_owned())
        .collect()
}

fn sample_displacement(modes: &[DVector<f64>], n_res: usize) -> DisplacementField {
    let mut rng = rand::thread_rng();
    let mut disp = vec![[0.0f32; 3]; n_res];
    if modes.is_empty() {
        return disp;
    }

    let coeffs: Vec<f64> = (0..modes.len())
        .map(|_| rng.gen_range(-1.0f64..1.0f64))
        .collect();

    for i in 0..n_res {
        for axis in 0..3 {
            let idx = 3 * i + axis;
            let mut value = 0.0f64;
            for (mode, c) in modes.iter().zip(coeffs.iter()) {
                value += *c * mode[idx];
            }
            disp[i][axis] = value as f32;
        }
    }
    disp
}

fn displacement_rmsd(disp: &DisplacementField) -> f32 {
    if disp.is_empty() {
        return 0.0;
    }
    let mut sum = 0.0f32;
    for d in disp {
        sum += d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    }
    (sum / (disp.len() as f32)).sqrt()
}

fn scale_displacement_to_rmsd(disp: &mut DisplacementField, target_rmsd: f32) {
    let current = displacement_rmsd(disp);
    if current <= f32::EPSILON || target_rmsd <= 0.0 {
        return;
    }
    let scale = target_rmsd / current;
    for d in disp.iter_mut() {
        d[0] *= scale;
        d[1] *= scale;
        d[2] *= scale;
    }
}

fn apply_displacement_to_vector(
    x: &mut [Option<f32>],
    y: &mut [Option<f32>],
    z: &mut [Option<f32>],
    disp: &DisplacementField,
) {
    for i in 0..disp.len() {
        if let (Some(xi), Some(yi), Some(zi)) = (x[i], y[i], z[i]) {
            x[i] = Some(xi + disp[i][0]);
            y[i] = Some(yi + disp[i][1]);
            z[i] = Some(zi + disp[i][2]);
        }
    }
}

fn apply_residue_displacement(
    structure: &CompactStructure,
    disp: &DisplacementField,
) -> CompactStructure {
    let mut sampled = structure.clone();
    apply_displacement_to_vector(
        &mut sampled.ca_vector.x,
        &mut sampled.ca_vector.y,
        &mut sampled.ca_vector.z,
        disp,
    );
    apply_displacement_to_vector(
        &mut sampled.n_vector.x,
        &mut sampled.n_vector.y,
        &mut sampled.n_vector.z,
        disp,
    );
    apply_displacement_to_vector(
        &mut sampled.cb_vector.x,
        &mut sampled.cb_vector.y,
        &mut sampled.cb_vector.z,
        disp,
    );
    sampled
}

/// Generate an ensemble consisting of the original structure plus ANM-sampled conformations.
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

    let ca_coords = extract_ca_coordinates(query_structure)?;
    let hessian = build_anm_hessian(&ca_coords, ANM_CUTOFF_ANGSTROM);
    let modes = select_nontrivial_modes(hessian, nma_modes);

    if modes.is_empty() {
        return Ok(ensemble);
    }

    for _ in 0..num_confs {
        let mut disp = sample_displacement(&modes, query_structure.num_residues);
        scale_displacement_to_rmsd(&mut disp, target_rmsd);
        ensemble.push(apply_residue_displacement(query_structure, &disp));
    }

    Ok(ensemble)
}
