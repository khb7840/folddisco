# Folddisco 3.0 — update report

Web version: see the link in the release notes. Figures live in
`fd-branchbench/figure/`, `fd-branchbench/defaults/figure/` and
`fd-branchbench/index_expansion/figure/`; measurement detail is in
[feature_evaluation.md](feature_evaluation.md), branch state in [UPDATE_REPORT.md](UPDATE_REPORT.md).

## What 3.0 adds

| feature | flag | what it does |
| --- | --- | --- |
| Sensitive search | `query --sensitive` (`--expand-radius`) | a residue pair may fall in neighbouring distance/angle bins, several features at once |
| Amino acid substitution | `query --aa-subst blosum62\|group\|size`, `:*` per residue | matches chemically similar residues, scored below exact ones |
| Confident hit list | `query --confident` | keeps only full, low-RMSD matches |
| Index-time expansion | `index --expand-radius/--expand-distance/--expand-angle/--aa-subst` | stores the neighbourhood in the index instead of expanding the query |
| Novelty evidence | `query --novelty-mode` | one evidence row per query: coverage, best hit, RMSD |
| Multi-character chain IDs | `-q AA_250`, `10_250`, `--chain-sep` | mmCIF chains that are not a single letter |
| Deformation metrics | `drmsd`, `max_dist_deviation` | superposition-free deviation, as column, sort key and filter |
| Mapped lookup caches | automatic | `*.lookup.cache`, `*.fdcache` built on first use |

Indices built with 2.x are read unchanged.

## Benchmarks

| set | index | queries | answers | metric |
| --- | --- | --- | --- | --- |
| Motif | human AFDB, 23,391 structures | 8 commands over 4 motifs (3-23 residues) | zinc finger 761, MEROPS S01 124 accessions | precision/recall/F1, Sens@1FP |
| M-CSA | PDB-derived, 62,122 entries | 250 catalytic sites (3-21 residues), 492 held out | M-CSA homologues per site | Sens@1FP, average precision |
| Mutant M-CSA | same | the same 250 sites with one residue renamed to its closest BLOSUM62 alternative | unchanged | Sens@1FP |
| Motif-only index | 24,762 M-CSA motifs | 250 sites, published and mutant | same-entry motifs | Sens@1FP, top-1 |

Machine: 20 cores, local NVMe, warm cache. M-CSA timings run 5 queries in parallel with 4
threads each; motif timings are serial medians of 5 repeats.

## Compatibility

- Index files are unchanged; 2.x indices need no rebuild.
- `--nonrigid` is now `--sensitive`.
- Default ordering changed (see below), so hit lists come back in a different order than 2.x.
  On the eight motif commands the branch returns exactly the same structures as
  `origin/master`; the four `--per-structure` commands return them in a new order.
- The residue that ends a chain is no longer labelled with the next chain's ID, which changes
  matched-residue labels on multi-chain entries.
- `-d`/`-a` with several values now use only the widest one.

## Limits

- Sens@1FP is a k=1 metric: single accessions move it. Paired win/loss counts and average
  precision are reported alongside.
- M-CSA answer sets are homologous entries, so a partial match can be a true answer; the zinc
  and serine sets are defined by exact residue identity.
- `--confident` past 12 residues rests on four queries only.
- Index-time expansion was measured on one motif-only index.
- Single machine, one hash type (`PDBTrRosetta`, default binning).

## Reproduce

```bash
cd fd-branchbench
REPEATS=5 WARMUP=1 ./scripts/run_all.sh 250     # motif, M-CSA, mutant, --confident, figures
defaults/scripts/collect.sh <binary>            # raw output for sort/filter selection
defaults/scripts/analyze.py && defaults/scripts/plot.py
index_expansion/scripts/run_queries.sh          # motif-only index study
```
