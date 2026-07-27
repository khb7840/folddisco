# Non-rigid motif search: design and measurements

Motifs are rarely rigid. The same catalytic site in two homologs, or the same site
before and after a conformational change, keeps its residues in the same arrangement
while the distances and angles between them drift. Folddisco discretizes that geometry
into bins, so a pair whose features land on the far side of a bin boundary produces a
different hash and is missed.

This document records what was changed to address that, and what each change measured.

## Benchmark

Everything below uses the same setup:

- **Index**: human proteome, 23,391 AlphaFold structures (`index/h_sapiens_folddisco`,
  hash type `PDBTrRosetta`, default bins). Download:
  `aria2c https://opendata.mmseqs.org/folddisco/h_sapiens_folddisco.tar.lz4`
- **Query**: zinc-finger motif from `query/1G2F.pdb`, three residue selections —
  `F207,F212,F225,F229` (4 residues), `F207,F212,F225` (3), `F205-212,F223-230` (16)
- **Answers**: `data/zinc_answer.tsv`, 1,816 human proteins annotated as zinc fingers
- **Metric**: `TP@k FP` — true positives found while walking the ranked result list
  until *k* false positives. Rank-sensitive and independent of any output cutoff.
- **Machine**: 8 threads, warm page cache, medians over 9–15 repeats

Reproduce with:

```bash
folddisco query -i index/h_sapiens_folddisco -p query/1G2F.pdb -q F207,F212,F225,F229 \
  -t 8 --skip-match --per-structure --format-output tid > result.tsv
folddisco benchmark -r result.tsv -a data/zinc_answer.tsv -i index/h_sapiens_folddisco \
  --afdb-to-uniprot --fp 5 -f default
```

## Results

| motif | configuration | hits | TP@5FP | TP@10FP | TP@100FP | recall |
| --- | --- | --- | --- | --- | --- | --- |
| 4 residues | before these changes | 1645 | 464 | 649 | 782 | 0.496 |
| | default | 1607 | **578** | **691** | 789 | 0.495 |
| | `--nonrigid` | 1850 | **684** | **731** | 795 | 0.518 |
| | `--nonrigid --enm-sample --num-confs 10` | 2021 | **709** | **746** | **823** | **0.542** |
| 3 residues | before these changes | 595 | 94 | 113 | 268 | 0.197 |
| | default | 590 | 94 | 113 | 269 | 0.197 |
| | `--nonrigid` | 1248 | **108** | **145** | 301 | 0.484 |
| | `--nonrigid --enm-sample --num-confs 10` | 1417 | 98 | 118 | **304** | **0.509** |
| 16 residues | before these changes | 22491 | 559 | 585 | 663 | 0.9923 |
| | `--nonrigid` | 22739 | 514 | 597 | 663 | 0.9934 |

Long motifs are already saturated (recall 0.99) and gain nothing; the wins are on
short motifs, which is the common case.

### Runtime

| | prefilter (`--skip-match`) | with residue matching |
| --- | --- | --- |
| before these changes | 9.1 ms | 405 ms |
| default | 9.5 ms (1.04x) | 403 ms (0.99x) |
| `--nonrigid` | 9.6 ms (1.05x) | 470 ms (1.16x) |
| `--nonrigid --enm-sample --num-confs 10` | 32.8 ms (3.59x) | 531 ms (1.31x) |

Indexing costs 1.4% more (19.31 s -> 19.58 s over 27 structures) for the backbone
carbonyl carbon the ENM needs.

## What changed

### 1. Joint bin expansion — `--expand-radius` (`src/controller/expand.rs`)

The tolerance expansion used to perturb one feature dimension at a time, so a target
pair whose distance *and* angle both drifted across a boundary was unreachable.
`--expand-radius 2` searches pairs of dimensions together. `--nonrigid` is a preset for
it. Radius 3 measured no better than 2 and slightly worse on some motifs.

Three further corrections in the same module, all active at every setting:

- **Wide tolerances are sub-stepped.** A single jump of 2.5 bins skipped the bins in
  between; offsets are now spaced within one bin so the covered bins are contiguous.
- **Angles are pulled back into their domain.** Torsions wrap at ±π, `acos`-derived
  angles reflect at their ends. Previously a wide `-a` pushed the bin index past the end
  of its bit field on the hash types that discretize radians directly, and on sin-cos
  types it asked for a `sin` sign no real structure can produce. Sensitivity used to
  *fall* past `-a 30` (12 → 10 → 4 → 2 matches at `-a 30/45/60/90`); it is now monotone.
- **Substitutions compose with tolerance.** `164:H` used to be hashed only at the
  observed geometry, so the alternative residue got no distance or angle slack.

Distances are deliberately *not* clamped. A structure can hold a Cβ–Cβ distance beyond
`MAX_DIST` — only the Cα distance is checked against the cutoff — so the index itself
contains those out-of-window encodings. Indexing and querying share one discretizer, so
mirroring it exactly is what makes a perturbed query hash mean the same thing as an
index hash.

### 2. Rare-hash filter (`src/controller/count_query.rs`)

An expanded hash far rarer in the database than the observed hash it came from is
dropped. A rare hash carries a large IDF, so one spurious hit on it can outrank several
real ones. The threshold is a ratio between two IDFs from the same index, so it needs no
tuning per database size.

Sweeping the margin on the 4-residue motif with `--expand-radius 2`:

| margin (2^m × rarer) | TP@5FP | TP@10FP | recall |
| --- | --- | --- | --- |
| 0 | 509 | 636 | 0.503 |
| 1 | 681 | 728 | 0.515 |
| 2 | 681 | 729 | 0.515 |
| **3 (chosen)** | **684** | **731** | **0.518** |
| 5 | 655 | 692 | 0.519 |
| 8+ | 535 | 682 | 0.520 *(no filtering)* |

Margin 0 discards the most informative expansions and is worse than no filter. Margin 2
is worse than 3 on the 16-residue motif, so 3 is the floor.

An absolute IDF ceiling (`--max-idf 8`, i.e. "appears in ≥0.39% of the database")
measured better still at the top of the ranking — 704 vs 684 TP@5FP here, and
124/150/351 vs 108/145/301 on the 3-residue motif. It is not implemented: the relative
rule needs no flag. Worth revisiting if short-motif precision matters more than
simplicity.

### 3. Torsion-angle ENM query sampling — `--enm-sample` (`src/structure/nma.rs`)

The query is wiggled along its low-frequency torsional normal modes and the union of the
ensemble's hashes is searched in one pass. Tunable with `--num-confs`, `--nma-rmsd`,
`--nma-modes`.

Best deep recall measured (0.542 vs 0.518 on the 4-residue motif) at ~4x the prefilter
runtime, but it *loses* ground at the top of the ranking for the 3-residue motif
(98 vs 108 TP@5FP), which is why it is opt-in rather than part of `--nonrigid`.

Sampling is seeded per (conformer, attempt), so an ensemble is reproducible regardless
of thread scheduling. A search tool whose results change between runs cannot be
benchmarked.

This needed a backbone carbonyl carbon on `CompactStructure`, which the parser already
read for the glycine Cβ approximation but did not store. See the note below.

### 4. Deformation metrics (`src/structure/metrics.rs`)

`drmsd` (distance-matrix RMSD) and `max_dist_deviation` compare the internal distances
of a match instead of superposing it, so a motif whose halves swung apart on a hinge
keeps a low dRMSD where its RMSD is large. Available to `--sort-by`, `--format-output`
and `--drmsd`, plus `min_drmsd` per structure.

On a synthetic hinge (chain C of `4CHA` rotated, motif `B55-58,C193-196` spanning it),
`--nonrigid` recovers all 8 residues at 0.97 Å RMSD where the default search collapses
to a 6-residue mis-assignment at 7.0 Å.

Ranking honestly did not improve: sorting by dRMSD does not beat RMSD after
`node_count` on the benchmark. These are reporting and filtering metrics.

## Alternatives measured and rejected

Evaluated from `khb7840/folddisco` branches on the same benchmark.

| feature | verdict |
| --- | --- |
| ANM/NMA query sampling | Rejected. Never beat torsion-ENM on quality and ran at 1.4 s vs 9 ms baseline (~110x), most likely a full 3N×3N eigendecomposition where torsion space is far smaller. |
| Binary lookup cache | Conditional. 2.4x faster at `-t 1`, 1.4x *slower* at `-t 8`; crossover at 3–4 threads. The text path is `par_lines()` over an mmap and scales with threads, while the cache read is serial. Worth taking if the read is parallelised — the format is fixed-layout. |
| DMS / PAS / SOS ranking metrics | Rejected. With an identical 1635-structure candidate set none beats plain RMSD as the tiebreaker after `node_count` (RMSD 730/742/843/898 vs DMS 729/741/834/896, PAS 730/740/837/880, SOS 668/741/844/895). DMS is a normalised dRMSD, and that does not improve ranking either. |
| `--dist-ratio` (elastic, distance-proportional tolerance) | Kept but off by default. Raises total recall while losing true positives at every early-precision point: widening the tolerance most for long pairs adds the least specific matches. |
| Perturbation analytics (`analyze-perturb`) | Not a search feature — byte-identical search output. Judge on the diagnostics alone. |

## Known issue found but not fixed

`CompactStructure::build` never resets its `c` variable between residues, so the Cβ
approximation for a residue lacking Cβ can fall back on an **earlier residue's**
carbonyl carbon. Fixing it changes every hash and would invalidate published indices, so
the existing behaviour is preserved verbatim and the ENM reads a separate per-residue
value instead. Worth addressing on a future index-format change.

## Caveat on all of the above

One motif family. The ordering of configurations held across three residue selections of
the zinc finger, but that is three views of one site. A second answer set — serine
peptidases can be built from `data/serine_peptidases/info_serhisasp.tsv` — would be
needed before treating these numbers as general.
