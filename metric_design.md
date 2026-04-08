# Structural Metric Design: DMS, PAS, SOS

## Overview

The current motif quality design is based on three explicit geometry metrics:

- **DMS** (Distance Matrix Score): pairwise Cα distance consistency
- **PAS** (Pseudo-Bond Angle Score): local Cα backbone curvature consistency
- **SOS** (Side-chain Orientation Score): Cα→Cβ orientation consistency

These metrics are computed in existing structure modules (primarily
`src/structure/metrics.rs`) and propagated through retrieval/sorting/output.

## Definitions

### 1) DMS — Distance Matrix Score

For two aligned motifs with Cα coordinates \(Q\) and \(H\):

\[
s(i,j)=\max\left(0,1-\frac{|d_Q(i,j)-d_H(i,j)|}{D_{TOL}}\right),\quad D_{TOL}=2.0\ \text{Å}
\]

\[
DMS(Q,H)=\frac{1}{\binom{N}{2}}\sum_{i<j}s(i,j)
\]

- Range: \([0,1]\)
- Complexity: \(O(N^2)\)

### 2) PAS — Pseudo-Bond Angle Score

\[
PAS(Q,H)=\frac{1}{N-2}\sum_{k=0}^{N-3}\frac{1+\cos(\theta^Q_k-\theta^H_k)}{2}
\]

where \(\theta_k\) is the angle formed by \((C\alpha_k, C\alpha_{k+1}, C\alpha_{k+2})\).

- Range: \([0,1]\)
- Complexity: \(O(N)\)

### 3) SOS — Side-chain Orientation Score

Using unit vectors \(u_i=\widehat{C\beta_i-C\alpha_i}\):

\[
SOS(Q,H)=\frac{1}{M}\sum_{i\in valid}\frac{1+\left(u^Q_i\cdot u^H_i\right)}{2}
\]

where `valid` residues have non-degenerate Cα→Cβ vectors in both motifs.

- Range: \([0,1]\)
- Complexity: \(O(N)\)

## Input conventions

- DMS and PAS operate on aligned Cα coordinate arrays.
- SOS operates on aligned Cα and Cβ arrays.
- In retrieval, DMS/PAS/SOS are computed from interleaved `[CA, CB, ...]` pairs via
  `StructureSimilarityMetrics::calculate_dms_pas_sos_from_interleaved_ca_cb`.

## Edge-case policy

- Empty/mismatched coordinates in retrieval-level interleaved computation:
  DMS = PAS = SOS = 0.0.
- Direct standalone functions:
  - DMS returns 1.0 for \(N < 2\)
  - PAS returns 1.0 for \(N < 3\)
  - SOS returns 1.0 when no valid orientation vectors exist

## Integration status

- **Computation**: `src/structure/metrics.rs`
- **Retrieval path wiring**: `src/controller/retrieve.rs`
- **Sorting keys**: `src/controller/sort.rs` (`dms`, `pas`, `sos`)
- **Result columns**: `src/controller/result.rs` (`dms`, `pas`, `sos`)
- **CLI help text**: `src/cli/workflows/query_pdb.rs`
