# GRAMS: Geometric-Residue Alignment Motif Score

## Overview

GRAMS is a composite ranking metric for structural motif search results in folddisco.
It addresses a fundamental limitation of pure backbone RMSD metrics: two motifs can
have identical backbone geometry yet differ significantly in biochemical function if
their side-chain residue identities are dissimilar.

---

## Three Candidate Metric Concepts

### Concept 1 — GRAMS (selected)

**Theory**: Combines backbone geometric quality (TM-score), residue-level biochemical
compatibility (normalised BLOSUM62), and a steric-clash penalty.

**Mathematical formulation**:

```
GRAMS(Q, H) = α · TM(Q, H) + β · ResComp(Q, H) + γ · (1 − ClashFrac(H))
```

where

| Symbol | Description |
|--------|-------------|
| `TM(Q, H)` | TM-score computed from post-superimposition Cα distances |
| `ResComp(Q, H)` | Mean normalised BLOSUM62 score over matched residue pairs |
| `ClashFrac(H)` | Fraction of intra-hit Cα pairs closer than 3.0 Å |
| α, β, γ | Non-negative weights with α + β + γ = 1 (defaults: 0.5, 0.3, 0.2) |

**ResComp**:

```
ResComp(Q, H) = (1/N) · Σ_i  b62_norm(q_i, h_i)

b62_norm(a, b) = (BLOSUM62[a][b] − B62_MIN) / (B62_MAX − B62_MIN)   ∈ [0, 1]
```

where `B62_MIN = −4` (minimum value across the entire BLOSUM62 matrix, used as the
normalisation lower bound) and `B62_MAX = 11` (maximum value — the Trp self-score,
used as the normalisation upper bound).

**ClashFrac**:

```
ClashFrac(H) = |{(i, j) : i < j, d_Cα(h_i, h_j) < D_CLASH}| / C(N, 2)
D_CLASH = 3.0 Å
```

Because all three components are in [0, 1], GRAMS ∈ [0, 1].
A higher GRAMS indicates a better match.

**Time complexity**: O(N) for TM and ResComp; O(N²) for ClashFrac.
For small motifs (N ≤ 10) typical in folddisco, the N² term is negligible.

---

### Concept 2 — SP-Extended with Side-Chain Dihedral Deviation

**Theory**: Extends the SP (Sum-of-Pairs) contact-map score with a χ₁ side-chain
dihedral-angle agreement term to penalise rotameric mismatch.

**Mathematical formulation**:

```
SPX(Q, H) = SP(Q, H) − δ · (1/N) · Σ_i min(|Δχ₁_i|, π) / π
```

where `SP(Q, H)` is the standard SP-score, `Δχ₁_i` is the difference in χ₁ dihedral
angles for residue i, and δ is a weighting factor.

**Limitation**: Requires explicit side-chain coordinates and a reliable rotamer
library, increasing data requirements and complexity.

---

### Concept 3 — Weighted Distance-Physicochemical Score (WDPS)

**Theory**: Reweights RMSD contributions residue-by-residue using a physicochemical
similarity factor derived from hydrophobicity, charge, and size differences.

**Mathematical formulation**:

```
WDPS(Q, H) = 1 − sqrt( (1/N) · Σ_i  w_i · d_i² ) / d_max

w_i = Physico_sim(q_i, h_i)   (normalised to Σ w_i = N)
```

**Limitation**: The physicochemical weight scheme requires careful calibration and
does not capture co-evolutionary substitution tolerance as well as BLOSUM62 does.

---

## Rationale for Selecting GRAMS

1. **Computational efficiency**: All operations are O(N) or O(N²) with small constants;
   no rotamer libraries or external data files are needed.
2. **Biological soundness**: BLOSUM62 encodes evolutionary substitution frequencies,
   providing a principled biochemical similarity measure.  The clash penalty naturally
   filters steric impossibilities.
3. **Normalisation**: Every term is independently bounded in [0, 1], making scores
   directly comparable across motifs of different sizes.
4. **Pure Rust**: No C bindings, no external ML runtimes — the entire implementation
   is a single self-contained Rust module.

---

## Inputs

| Parameter | Type | Description |
|-----------|------|-------------|
| `ref_ca` | `&[[f32; 3]]` | Query Cα coordinates (post-superimposition, Å) |
| `model_ca` | `&[[f32; 3]]` | Hit Cα coordinates (post-superimposition, Å) |
| `ref_residues` | `&[u8]` | Query residue one-letter codes (ASCII uppercase) |
| `model_residues` | `&[u8]` | Hit residue one-letter codes (ASCII uppercase) |
| `weights` | `GramsWeights` | α, β, γ weights (must sum to 1.0) |

## Edge Cases

| Scenario | Handling |
|----------|----------|
| N = 0 | Returns 0.0 |
| Length mismatch | Returns 0.0 |
| Unknown residue code | BLOSUM62 lookup returns 0.0 (neutral) |
| Single residue | ClashFrac = 0.0 (no pairs); TM-score uses d₀ = 0.5 |
| NaN in coordinates | Propagates to 0.0 via `unwrap_or(0.0)` |
| All clashing | ClashFrac = 1.0, γ term contributes 0.0 |
