//! Native ANM/NMA conformational sampling for non-rigid query search.
//!
//! This module implements a lightweight ProDy-like ANM pipeline:
//! 1) build a 3N×3N Hessian from backbone node coordinates (N/CA/C, cutoff ENM),
//! 2) eigendecompose and discard rigid-body modes,
//! 3) sample random linear combinations of low-frequency modes,
//! 4) scale displacements to a target RMSD and apply them to N/CA/C atoms,
//! 5) reconstruct CB from perturbed backbone orientation.

use nalgebra::{DMatrix, DVector, SymmetricEigen};
use rand::Rng;
use std::fs::File;
use std::io::Write;

use crate::structure::coordinate::{approx_cb, Coordinate};
use crate::structure::core::CompactStructure;

const ANM_CUTOFF_ANGSTROM: f64 = 15.0;
const EIGENVALUE_EPS: f64 = 1e-8;

type DisplacementField = Vec<[f32; 3]>;

fn extract_backbone_coordinates(structure: &CompactStructure) -> Result<Vec<Coordinate>, String> {
    let mut coords = Vec::with_capacity(structure.num_residues * 3);
    for i in 0..structure.num_residues {
        let ca = structure
            .ca_vector
            .get_coord(i)
            .ok_or_else(|| format!("Missing CA coordinate at residue index {}", i))?;
        let n = structure.n_vector.get_coord(i).unwrap_or(ca);
        let c = structure.c_vector.get_coord(i).unwrap_or(ca);
        coords.push(n);
        coords.push(ca);
        coords.push(c);
    }
    Ok(coords)
}

fn build_anm_hessian(coords: &[Coordinate], cutoff: f64) -> DMatrix<f64> {
    let n = coords.len();
    let size = 3 * n;
    let mut h = DMatrix::<f64>::zeros(size, size);
    let cutoff2 = cutoff * cutoff;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = (coords[j].x - coords[i].x) as f64;
            let dy = (coords[j].y - coords[i].y) as f64;
            let dz = (coords[j].z - coords[i].z) as f64;
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
        // The first six ANM eigenmodes are rigid-body motions (3 translations + 3 rotations).
        // They do not represent internal conformational changes and are excluded from sampling.
        .skip(6)
        .filter(|(val, _)| *val > EIGENVALUE_EPS)
        .take(mode_count)
        .map(|(_, idx)| eig.eigenvectors.column(idx).into_owned())
        .collect()
}

fn sample_displacement(modes: &[DVector<f64>], n_nodes: usize) -> DisplacementField {
    let mut rng = rand::thread_rng();
    let mut disp = vec![[0.0f32; 3]; n_nodes];
    if modes.is_empty() {
        return disp;
    }

    let coeffs: Vec<f64> = (0..modes.len())
        .map(|_| rng.gen_range(-1.0f64..1.0f64))
        .collect();

    for i in 0..n_nodes {
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
    debug_assert_eq!(disp.len(), structure.num_residues * 3);

    let mut n_disp = vec![[0.0f32; 3]; structure.num_residues];
    let mut ca_disp = vec![[0.0f32; 3]; structure.num_residues];
    let mut c_disp = vec![[0.0f32; 3]; structure.num_residues];
    for i in 0..structure.num_residues {
        n_disp[i] = disp[3 * i];
        ca_disp[i] = disp[3 * i + 1];
        c_disp[i] = disp[3 * i + 2];
    }

    apply_displacement_to_vector(
        &mut sampled.ca_vector.x,
        &mut sampled.ca_vector.y,
        &mut sampled.ca_vector.z,
        &ca_disp,
    );
    apply_displacement_to_vector(
        &mut sampled.n_vector.x,
        &mut sampled.n_vector.y,
        &mut sampled.n_vector.z,
        &n_disp,
    );
    apply_displacement_to_vector(
        &mut sampled.c_vector.x,
        &mut sampled.c_vector.y,
        &mut sampled.c_vector.z,
        &c_disp,
    );

    // Reconstruct CB from backbone orientation after perturbation.
    for i in 0..sampled.num_residues {
        if let (Some(n), Some(ca), Some(c)) = (
            sampled.n_vector.get_coord(i),
            sampled.ca_vector.get_coord(i),
            sampled.c_vector.get_coord(i),
        ) {
            let cb = approx_cb(&ca, &n, &c);
            sampled.cb_vector.x[i] = Some(cb.x);
            sampled.cb_vector.y[i] = Some(cb.y);
            sampled.cb_vector.z[i] = Some(cb.z);
        }
    }
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

    // Fewer than 3 residues cannot provide meaningful internal normal modes for non-rigid sampling.
    if num_confs == 0 || nma_modes == 0 || query_structure.num_residues < 3 {
        return Ok(ensemble);
    }

    let backbone_coords = extract_backbone_coordinates(query_structure)?;
    let hessian = build_anm_hessian(&backbone_coords, ANM_CUTOFF_ANGSTROM);
    let modes = select_nontrivial_modes(hessian, nma_modes);

    if modes.is_empty() {
        return Ok(ensemble);
    }

    for _ in 0..num_confs {
        let mut disp = sample_displacement(&modes, query_structure.num_residues * 3);
        scale_displacement_to_rmsd(&mut disp, target_rmsd);
        ensemble.push(apply_residue_displacement(query_structure, &disp));
    }

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
