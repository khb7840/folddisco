# GRAMS: Geometric Alignment Motif Score

## Overview

GRAMS is a fully geometry-focused composite ranking metric for structural motif
search results in folddisco.  It scores hits using only Cα coordinates — no
sequence identity or residue type information is needed.

---

## Three Candidate Metric Concepts (Geometry-Focused)

### Concept 1 — GRAMS with TM + DMS + PAS (selected)

**Theory**: Combines three complementary, purely geometric quality measures:
global backbone quality (TM-score), internal distance-matrix consistency (DMS),
and local backbone curvature agreement (PAS).

**Mathematical formulation**:

```
GRAMS(Q, H) = α · TM(Q, H) + β · DMS(Q, H) + γ · PAS(Q, H)
```

where

| Symbol | Description |
|--------|-------------|
| `TM(Q, H)` | TM-score computed from post-superimposition Cα distances |
| `DMS(Q, H)` | Distance Matrix Score: pairwise Cα distance agreement |
| `PAS(Q, H)` | Pseudo-Bond Angle Score: Cα-Cα-Cα bond-angle agreement |
| α, β, γ | Non-negative weights with α+β+γ = 1 (defaults: 0.5, 0.3, 0.2) |

**DMS**:

```
s(i,j) = max(0, 1 − |d_Q(i,j) − d_H(i,j)| / D_TOL)     D_TOL = 2.0 Å
DMS(Q, H) = (1 / C(N,2)) · Σ_{i<j} s(i,j)
```

**PAS**:

```
PAS(Q, H) = (1/(N−2)) · Σ_{k=0}^{N−3} (1 + cos(θ_Q_k − θ_H_k)) / 2
```

where θ_k is the angle at Cα_{k+1} formed by the triplet (Cα_k, Cα_{k+1}, Cα_{k+2}).

**Time complexity**: TM and PAS are O(N); DMS is O(N²).
For small motifs (N ≤ 10) common in folddisco, N² ≤ 45 distance evaluations.

---

### Concept 2 — Torsion-Angle Profile Similarity

**Theory**: For each consecutive pair of backbone Cα positions, derive a
pseudo-torsion angle from four successive Cα atoms and compare the resulting
torsion-angle profiles between query and hit.

**Mathematical formulation**:

```
TAPS(Q, H) = (1/(N−3)) · Σ_{k=0}^{N−4} (1 + cos(τ_Q_k − τ_H_k)) / 2
```

where τ_k is the torsion angle of the quadruple (Cα_k, Cα_{k+1}, Cα_{k+2}, Cα_{k+3}).

**Limitation**: Requires N ≥ 4; carries less information for very short motifs
(3–4 residues) where only 0–1 torsion angles are available.

---

### Concept 3 — Voronoi-Weighted Shape Descriptor

**Theory**: Tessellate the motif Cα point cloud and compare normalised Voronoi
cell volumes between query and hit as a rotation/translation-invariant shape
descriptor.

**Limitation**: Computing Voronoi tessellations in 3D requires additional
algorithmic complexity and a third-party crate (e.g. `voro-rs`), conflicting
with the constraint against heavy external dependencies.

---

## Rationale for Selecting GRAMS (Concept 1)

1. **Purely geometric**: All three terms depend only on Cα coordinates; no
   sequence information or residue type tables are consulted.
2. **Complementary coverage**: TM captures global alignment quality; DMS detects
   distortion of the internal distance geometry; PAS flags incorrect local
   backbone curvature.
3. **Normalisation**: Every term is independently bounded in [0, 1], making scores
   directly comparable across motifs of different sizes.
4. **Computational efficiency**: O(N²) at most — negligible for motifs of size ≤ 10.
5. **Pure Rust**: No C bindings, no external ML runtimes.

---

## Inputs

| Parameter | Type | Description |
|-----------|------|-------------|
| `ref_ca` | `&[[f32; 3]]` | Query Cα coordinates (post-superimposition, Å) |
| `model_ca` | `&[[f32; 3]]` | Hit Cα coordinates (post-superimposition, Å) |
| `weights` | `GramsWeights` | tm, distance, angle weights (auto-normalised) |

## Edge Cases

| Scenario | Handling |
|----------|----------|
| N = 0 | Returns 0.0 |
| Length mismatch | Returns 0.0 |
| N = 1 | DMS = 1.0, PAS = 1.0 (no pairs/triplets); TM uses d₀ = 0.5 |
| N = 2 | PAS = 1.0 (no triplets); DMS and TM computed normally |
| Coincident atoms | `pseudo_bond_angle` returns 0.0 (degenerate check) |
| NaN coordinates | Propagates; guarded by `clamp(-1,1)` in acos call |
