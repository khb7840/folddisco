# Feature evaluation

Measurement record for the `feature-integration` branch. It ships:

1. **Sensitive search**: joint bin expansion, `--sensitive` / `--expand-radius` (§2–§4)
2. **Binary lookup cache**: automatic, no flag (§5)
3. **Novelty evidence mode**: `--novelty-mode`, one evidence row per query, no verdict (§10)
4. **Amino acid substitution schemes**: `--aa-subst blosum62|group|size` and per-residue `:*`, with substitution-aware scoring (§14)
5. **Index-time expansion**: `folddisco index --expand-radius/--expand-distance/--expand-angle/--aa-subst` (§14)

Also: superposition-free deformation metrics (`drmsd`, `min_drmsd`, §9) and a fixed `-q F204-F215` parse (§12).

Sensitivity was measured with the repository author's protocol (commands, answer sets,
Sens@1FP and F1; §2) and on the M-CSA catalytic-site benchmark over a rebuilt 62,122-entry
PDB index (§3). Conclusions that changed the code: the rare-hash IDF filter was removed (§8),
and `--sensitive` needs `--max-node` to pay off (§2.1) and costs more at scale (§2.3).
Largest open question: §3.2.

## 1. Method

- **Index**: human proteome, 23,391 AlphaFold structures (`index/h_sapiens_folddisco`,
  `PDBTrRosetta`, default bins) → 20,504 UniProt accessions.
  `aria2c https://opendata.mmseqs.org/folddisco/h_sapiens_folddisco.tar.lz4`
- **Commands.** A/B/D scored against the author's zinc-finger set (761 accessions), C against
  their MEROPS S01 serine-peptidase set:

```
A1  query -p query/1G2F.pdb -q F207,F212,F225,F229 -i $IDX -t 12 --covered-node 3 --skip-match
A2  query -p query/1G2F.pdb -q F207,F212,F225,F229 -i $IDX -t 12 --covered-node 3 --max-node 4 --rmsd 1.0 --per-structure
B1  query -p query/1G2F.pdb -q F207,F225,F229      -i $IDX -t 12 --covered-node 3 --skip-match
B2  query -p query/1G2F.pdb -q F207,F225,F229      -i $IDX -t 12 --covered-node 3 --max-node 3 --rmsd 1.0 --per-structure
C1  query -p query/4CHA.pdb -q B57,B102,C195       -i $IDX -t 12 --covered-node 3 --skip-match
C2  query -p query/4CHA.pdb -q B57,B102,C195       -i $IDX -t 12 --max-node 3 --rmsd 1.0 --per-structure
D1  query -p query/1G2F.pdb -q F204-215,F222-232   -i $IDX -t 1 --top 800 --skip-match
D2  query -p query/1G2F.pdb -q F204-215,F222-232   -i $IDX -t 1 --top 800 --per-structure --max-node 15
```

- **Metrics.** `--fp 1`: field 10 = TP@1FP, field 15 = Sens@1FP. Without `--fp`: field 17 is a
  set F1 (no rank dependence). F1 decides; Sens@1FP is reported (§6).
- **Answer sets.** Zinc: 761. Serine: 130 lines, `answer_len` 124 (P20231, P48740, Q15661,
  Q2TV78, Q5K4E3, Q7RTY7 carry two S01 ids and are deduplicated). `--afdb-to-uniprot`
  throughout; `hits` counts raw result lines.
- **Robustness.** 200 paired replicates dropping 5% of answer accessions. Sign consistency =
  fraction of replicates keeping the unperturbed sign (1.00 robust, 0.5 noise).
- **Machine.** 20 cores, warm cache, timing = median of 13 interleaved repeats. Searches are
  deterministic.

```bash
IDX=index/h_sapiens_folddisco
folddisco query -i $IDX -p query/1G2F.pdb -q F207,F212,F225,F229 \
  -t 12 --covered-node 3 --max-node 4 --rmsd 1.0 --per-structure > result.tsv
folddisco benchmark -r result.tsv -a <answers>.tsv -i $IDX --afdb-to-uniprot --fp 1   # Sens@1FP
folddisco benchmark -r result.tsv -a <answers>.tsv -i $IDX --afdb-to-uniprot          # P/R/F1
python3 scripts/eval_metrics.py result.tsv <answers>.tsv $IDX.lookup                  # + TP@kFP
```

## 2. `--sensitive`

`default` = this branch without flags, byte-identical to master (2a756d9) on all eight commands (§4).

| cmd | config | hits | TP@1FP | Sens@1FP | precision | recall | F1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **A1** 4-res, prefilter | master = default | 814 | 31 | 0.0407 | 0.9100 | 0.9435 | **0.9265** |
| | `--sensitive` | 878 | **42** | 0.0552 | 0.8787 | **0.9711** | 0.9226 |
| **A2** 4-res, matched | master = default | 746 | 31 | 0.0407 | 0.9655 | 0.9198 | 0.9421 |
| | `--sensitive` | 802 | **42** | 0.0552 | 0.9584 | **0.9698** | **0.9641** |
| **B1** 3-res, prefilter | master = default | 767 | **43** | 0.0565 | 0.9489 | 0.9277 | 0.9382 |
| | `--sensitive` | 826 | 22 | 0.0289 | 0.9328 | **0.9671** | **0.9497** |
| **B2** 3-res, matched | master = default | 756 | **43** | 0.0565 | 0.9591 | 0.9251 | 0.9418 |
| | `--sensitive` | 813 | 22 | 0.0289 | 0.9485 | **0.9671** | **0.9577** |
| **C1** Ser-His-Asp, prefilter | master = default | 107 | 29 | 0.2339 | 0.9533 | 0.8226 | 0.8831 |
| | `--sensitive` | 114 | **32** | 0.2581 | **0.9561** | **0.8790** | **0.9160** |
| **C2** Ser-His-Asp, matched | master = default | 113 | 29 | 0.2339 | 0.9640 | 0.8629 | 0.9106 |
| | `--sensitive` | 113 | **32** | 0.2581 | 0.9640 | 0.8629 | 0.9106 |
| **D1** 23-res segments, prefilter | master = default | 800 | **26** | 0.0342 | 0.8819 | 0.8830 | **0.8825** |
| | `--sensitive` | 800 | 25 | 0.0329 | 0.8729 | 0.8752 | 0.8740 |
| **D2** 23-res segments, matched | master = default | 697 | **26** | 0.0342 | 0.9722 | 0.8739 | **0.9204** |
| | `--sensitive` | 705 | 25 | 0.0329 | 0.9594 | 0.8686 | 0.9117 |

Recall always rises (+0.026 to +0.056) and precision always falls:

| cmd | Δ F1 | sign consistency | Δ Sens@1FP (TP) | sign consistency |
| --- | --- | --- | --- | --- |
| A1 4-res prefilter | −0.0039 | **1.00** | +11 | 0.53 |
| A2 4-res matched | **+0.0219** | **1.00** | +11 | 0.53 |
| B1 3-res prefilter | **+0.0115** | **1.00** | −21 | 0.56 |
| B2 3-res matched | **+0.0159** | **1.00** | −21 | 0.56 |
| C1 Ser-His-Asp prefilter | **+0.0328** | **1.00** | +3 | 0.85 |
| C2 Ser-His-Asp matched | 0.0000 | — | +3 | 0.85 |
| D1 23-res prefilter | **−0.0085** | **1.00** | −1 | 0.50 |
| D2 23-res matched | **−0.0087** | **1.00** | −1 | 0.50 |

Every F1 delta is robust; no Sens@1FP delta is.

### 2.1 Short motifs only, with `--max-node`

4-residue query, author's filters added one at a time (F1, 761 set):

| variant | default | `--sensitive` | Δ |
| --- | --- | --- | --- |
| `--skip-match`, no filters | 0.6511 | 0.5979 | **−0.0532** |
| `--covered-node 3 --skip-match` | 0.9265 | 0.9226 | −0.0039 |
| `--covered-node 3 --rmsd 1.0` | 0.9312 | 0.9360 | +0.0048 |
| `--covered-node 3 --max-node 4` | 0.9397 | 0.9609 | **+0.0212** |
| `--covered-node 3 --max-node 4 --rmsd 1.0` (A2) | 0.9421 | 0.9641 | **+0.0220** |

- `--max-node`, not `--rmsd`, turns the extra recall into F1.
- `--covered-node 3` dominates the protocol: dropping it costs the default 0.28 F1.
- Long segment queries lose (D1/D2 −0.009 F1, sign consistency 1.00): already saturated.

### 2.2 Not tolerance widening

Across a `-d`/`-a` grid on all eight commands, widening never reaches `--sensitive`'s F1 (best
cell +0.0008 over default, worst −0.14) and loses recall (A2: default 0.920, `-d 2.0 -a 20`
0.894, `--sensitive` 0.970). `-d`/`-a` set how far one feature moves; `--sensitive` how many move
at once. Wide single-dimension sub-steps also exhaust `MAX_HASHES_PER_PAIR = 4096` before any
joint variant is emitted. §3 agrees on a different index and metric.

### 2.3 Runtime

13 repeats after 2 warm-ups, wall clock, `-t 12`:

| config | A1 prefilter | A2 matched |
| --- | --- | --- |
| master | 9.1 ms (8.5–9.8) | 158.9 ms (154–177) |
| default | 9.1 ms (8.6–10.3) | 157.4 ms (148–174) |
| `--sensitive` | 9.9 ms (8.6–10.4) | 174.7 ms (168–183) |

+11% matched here; the prefilter difference is within spread. On M-CSA (§3) it is ~73%
slower: cost scales with candidate-pool inflation, so with motif and index size.

## 3. M-CSA on a PDB-derived index

M-CSA catalytic-site benchmark, rebuilt `mcsa_subset` index (62,122 PDB entries / 131,384
rows), paper metric: mean Sens@1FP at `--top 6000`.

**Rebuild fidelity.** All 131,384 lookup paths match the author's lookup; 99.90% of rows have
identical residue counts (413 of 99,719,731 residues differ, from 69 re-versioned entries);
`.type` matches on all fields. Note: `-y default` gives `PDBTrRosetta`; `-y pdb` gives
`PDBMotifSinCos`.

**Gate** (upstream master, `--top 6000`):

| | n | mean Sens@1FP | median |
| --- | --- | --- | --- |
| author's published figure | 752 | 0.4323 | 0.3940 |
| this rebuild, master | 752 | **0.4507** | **0.4000** |

Per-query direction is balanced (295 higher, 258 lower, 26.5% identical): corpus scatter.

**Grid** (frozen seeded 250-query subset, paired):

| config | mean Sens@1FP | median | zero-TP queries | win / loss / tie vs master |
| --- | --- | --- | --- | --- |
| master `--top 6000` | 0.4475 | **0.4032** | 10 | — |
| `--sensitive` | **0.4536** | 0.3810 | 12 | **94 / 42 / 114** |
| wide `-d`/`-a`, same cap | 0.3599 | 0.2608 | 15 | 38 / 137 / 75 |
| wide `-d`/`-a` + `--sensitive` | 0.3538 | 0.2500 | **43** | 55 / 129 / 66 |

Branch default was byte-identical to master on 25 of 25 M-CSA queries (before 5ae7321, §14).

### 3.1 `--sensitive` here

Mean and paired count improve, median falls: a modest mixed win. This protocol has no
`--max-node`. Wall time 5.53 → 9.56 s per query (250 queries, 5 processes × 4 threads).

### 3.2 Open: the uncapped sensitive preset

The wide rows are our reconstruction under `--top 6000`, not the author's sensitive preset,
which ran uncapped (per-query `result_len` up to 44,412, mean 5,671 vs 1,623 default).
Supported claim: **under a fixed top-N budget, `--sensitive` beats widening and the two do not
compose.** Against the uncapped preset: untested (>10 min per query uncapped vs 10–23 s capped).

### 3.3 IDF alone does not rank at PDB scale

On the full PDB index (230,655 structures, 20 GB mapped), a matched query sorted by raw IDF
puts sparse partial matches first (`2n8a`, node_count 2, IDF 26,829); the prefilter ranking of
the same query is correct. Use a node-count floor. Prefilter 0.11 s warm; matched 11.3 s
(peak RSS 770 MB).

## 4. Default path

With the rare-hash filter removed, the default is byte-identical (`cmp`) to master on all
eight commands (A1 814 lines, A2 746, B1 767, B2 756, C1 107, C2 113, D1 800, D2 697).
Per-match output (M-CSA) differs from master since 5ae7321 by design (§14).

Not a no-op at wider tolerance (A1):

| tolerance | master vs default | master TP@1FP / recall / F1 | default TP@1FP / recall / F1 |
| --- | --- | --- | --- |
| `-d 0.5 -a 5` | byte-identical | 31 / 0.9435 / 0.9265 | 31 / 0.9435 / 0.9265 |
| `-d 2.0 -a 20` | **differs** | 22 / 0.8476 / 0.8764 | **29 / 0.9290 / 0.8927** |

Master queries only the two extreme bins of a tolerance wider than one bin; the branch
sub-steps so covered bins are contiguous. Other fixes in `expand.rs`:

- Angles are wrapped (torsions) or reflected (`acos` angles) into their domain; sensitivity no
  longer falls past `-a 30` (was 12 → 10 → 4 → 2 matches at `-a 30/45/60/90`).
- Substitutions compose with the geometric tolerance.
- Distances are not clamped, so query encoding mirrors the index.

## 5. Binary lookup cache

Written on the first text parse, used afterwards, falls back to text on any mismatch. No
effect on ranking. See UPDATE_REPORT.md §3 for the current mmap layout.

PDB index (230,655 entries): ~10 ms saved (0.12 s parse vs 0.11 s decode, 6 runs), ~8% of a
0.11 s prefilter, negligible for an 11.3 s matched query. Cache 14,069,999 B (61 B/entry vs 48
text). Synthetic AFDB-style lookups, medians of 5 (v2 decode cache):

| entries | file size | threads | text | cache | speedup |
| --- | --- | --- | --- | --- | --- |
| 23,391 (human index) | 1.19 MB | 1 | 5.28 ms | 1.23 ms | 4.3x |
| | | 8 | 1.23 ms | 0.97 ms | 1.27x |
| | | 20 | 1.94 ms | 1.87 ms | 1.04x |
| 10^5 | 4.9 MB | 1 | 13.16 ms | 3.41 ms | 3.9x |
| | | 8 | 3.70 ms | 2.19 ms | 1.69x |
| | | 20 | 3.35 ms | 2.84 ms | 1.18x |
| 10^6 | 50.8 MB | 1 | 145.9 ms | 67.7 ms | 2.2x |
| | | 8 | 40.7 ms | 31.4 ms | 1.30x |
| | | 20 | 33.1 ms | 30.6 ms | 1.08x |
| 10^7 | 528 MB | 1 | 1498.6 ms | 666.8 ms | 2.25x |
| | | 8 | 346.8 ms | 262.6 ms | 1.32x |
| | | 20 | 277.7 ms | 235.7 ms | 1.18x |

- Never slower; 10^6–10^7 rows are projections (no shipping index is that large).
- First load pays the write (+0.63 ms at 23,391 / 8 threads, +38 ms at 10^6, +312 ms at 10^7);
  breaks even after ~4 loads at 8 threads, 9 at 20.
- Validated on source length, mtime (ns) and blob length. A same-length rewrite with a forged
  mtime is not detected (no content hash, by design); delete the cache in that case.
- Concurrent first loads are safe: content is deterministic and a partial file is rejected
  (2 threads, 40 trials, 0 corrupt).

## 6. Metrics: F1 decides, Sens@1FP is reported

1. F1 deltas hold sign in 1.00 of replicates; Sens@1FP in 0.50–0.85.
2. At k=1 a change is one accession (an earlier A1 31-vs-27 gap was Q8NB15 moving rank 32 → 28,
   a true positive in the other zinc set).
3. TP@1FP cannot see `--covered-node`, `--max-node`, `--rmsd` or `--skip-match` (31 default /
   42 `--sensitive` in all five §2.1 variants).

F1 without `--fp` is rank-blind: on C2 `--sensitive` returns the same 111 accessions reordered
(P49862, P08246, Q9P0G3 promoted); F1 stays 0.9106, Sens@1FP 29 → 32. For ranking use
`--fp 10` or `--fp 100` (TP@100FP sign consistency 1.00 vs 0.28–0.49 at TP@5FP).

## 7. The answer set can flip both metrics

Same lists scored against the older `data/zinc_answer.tsv` (1817 accessions):

| cmd | config | 761-set TP@1FP / F1 | 1817-set TP@1FP / F1 |
| --- | --- | --- | --- |
| A1 | default | 31 / 0.9265 | 90 / 0.5902 |
| A1 | `--sensitive` | **42** / 0.9226 | 61 / **0.6005** |
| A2 | default | 31 / 0.9421 | 89 / 0.5673 |
| A2 | `--sensitive` | **42** / **0.9641** | 61 / **0.5922** |

On A1 both metrics reverse with the answer set alone: F1 is precision-limited on the 761 set
(recall 0.87–0.97) and recall-limited on the 1817 set (0.40–0.44). Always name the answer set.

## 8. Rejected alternatives

| feature | verdict |
| --- | --- |
| **Rare-hash IDF filter** (`EXPANDED_HASH_IDF_MAX_EXCESS`) | Removed. Margins 0/2/3/5/8/off on eight commands at both radii: shipped 5.0 worse-or-equal in 16/16 cells on both metrics. F1 harm 0.0001–0.0020 (sign consistency 0.91–1.00); no runtime gain (158.7 vs 159.8 ms). Branch `rare-hash-idf-filter`. |
| **Torsion-angle ENM sampling (`--enm-sample`)** | Removed. vs `--sensitive`: matched F1 −0.0013 / +0.0006 / 0.0000; C1 prefilter −0.1206 F1 (precision 0.9561 → 0.7630). Runtime 3.34x (zinc) to 46.7x (triad). Conformer count 2/5/12 byte-identical on zinc; on the triad F1 falls 0.8240 → 0.7675 from 2 to 40. Added 17 crates. Branch `enm-torsion-sampling`. |
| ANM/NMA sampling | Rejected. Never beat torsion-ENM; ~110x slower (1.4 s vs 9 ms). |
| DMS / PAS / SOS composite ranking | Rejected. Same 1635 candidates; none beat RMSD as tiebreaker after `node_count`. |
| `--dist-ratio` | Removed. More total recall, fewer true positives at every early-precision depth. |
| `--max-idf` | Not implemented; moot (see §3.3 for the node-count floor). |
| `analyze-perturb` | Diagnostics only; search output byte-identical. |

The last five rows predate the author's protocol; their figures in other units were dropped.

## 9. Deformation metrics

`drmsd` and `max_dist_deviation` compare internal distances of a match (no superposition);
available in `--sort-by`, `--format-output`, `--drmsd`, and per structure as `min_drmsd`.

Synthetic hinge (4CHA chains C and G rotated about the chain-C centroid, motif
`B55-58,C193-196`): at 20° the default returns a 3-residue mis-assignment at 2.81 Å while
`--sensitive` recovers all 8 residues at 1.31 Å. Below 10° the default suffices; at 30° radius 2
is not enough. Sorting by dRMSD did not beat RMSD after `node_count`.

## 10. Novelty evidence mode

`--novelty-mode` prints one tab-separated row per query and makes no novelty call. The earlier
verdict (KNOWN / PARTIAL_MATCH / NOVEL / FILTERED_OUT / NO_HASHES, with `--novelty-coverage`
and `--novelty-rmsd`) was removed. Columns:

`query_id, status, candidates, hits, index_coverage, best_hit, best_coverage, best_rmsd, query_residues`

- `status`: `ok`, `no_candidates`, `filtered_out` (candidates existed, filters kept none; also
  warned on stderr), `no_hashes` (residues farther apart than the index cutoff).
- `candidates` / `index_coverage`: from the inverted index before any filter.
- `hits` / `best_*`: after filters and matching; `best_rmsd` is `NA` with `--skip-match`.

Measured facts that still apply (PDB index, human chymotrypsin C triad
`data/AF-P17538-F1-model_v4.pdb -q A57,A102,A195`, `--top 6000`):

| flags | best coverage | best RMSD | after filters |
| --- | --- | --- | --- |
| `--skip-match` | 0.6667 | NA | kept |
| `--max-node 2` | 0.6667 | 0.0205 | kept |
| `--max-node 3` | — | — | all removed (index held 0.6667) |

- A high `--max-node` discards partial matches, here a 2-of-3 match at 0.0205 Å to the query's
  own family. Screen with `--skip-match` or a low `--max-node`.
- Duplicate residues count twice in coverage denominators (`A1-A5,A3-A7`: 7 distinct, 10
  entries; a perfect match reads 0.70). Warned, not deduplicated, to keep the default path
  byte-identical.
- Earlier verdict data (100 random motifs, 32 of 96 "known" best hits above 2.0 Å, max 10.4 Å)
  shows why RMSD must be read alongside coverage.

## 11. Independent serine answer set

`data/serine_answer.tsv` (89 accessions) from `scripts/build_serine_answer.py`: S1 clan by
PROSITE PS00134/PS00135, sequence-defined and so independent of Folddisco. Tied to
`index/h_sapiens`; regenerate for other indices. Agreed with MEROPS on the `--sensitive`
direction. The author's `serinepeptidase_answer.tsv` (124 accessions) and
`zincfinger_answer.tsv` (761) are the sets of record for §2–§7 and are not committed.

## 12. Known issues

- `CompactStructure::build` does not reset `c` between residues, so a missing-Cβ fallback can
  use an earlier residue's carbonyl carbon. Fixing it changes every hash and published index.
- A 3-column Foldcomp-native lookup panics the text lookup parser (columns 3–4 unwrapped),
  e.g. `data/foldcomp/example_db.lookup`.

Fixed: `-q F204-F215` (range end repeating the chain) panicked; it now parses and a different
chain on the end is rejected with a message.

## 13. Limits

- Two indices, three answer sets; the full PDB index was used only for §3.3, §5 and §10.
- Uncapped sensitive-preset comparison open (§3.2).
- M-CSA grid is a 250-query subset; only master ran all 755.
- Fragility analysis perturbs answer sets, not indices or queries.
- `-d 2.0 -a 20` in §4 is outside the author's protocol.
- §9, §10 and the last five rows of §8 predate the author's protocol.
- Substitution and index-time expansion were measured on M-CSA only (§14).
- Single machine, 20 cores, local NVMe, warm cache.

## 14. Amino acid substitution

**Benchmarks** (M-CSA index, `--top 6000`, Sens@1FP; paired over queries, 95% bootstrap CI):

- *exact*: the M-CSA queries as published. Substitution can only cost here.
- *mutant*: each query with one residue renamed to its best BLOSUM62 alternative (His→Tyr,
  Asp→Glu, …; `fd-branchbench/scripts/make_mutant_queries.py`), same answer sets. Renaming
  drops q250 from 0.4465 to 0.2851; substitution should win it back.
- M-CSA answer sets hold aligned homologues, not identical residues: 4,468 q250 answers
  differ from the query at one catalytic position.

**Why the first scheme failed.** The per-pair 4096-hash cap never binds at default
tolerances (radius 2: 50 bins × ≤36 residue pairs). The loss came from scoring: every
substituted hash added its full IDF, a rare substitution outscored the exact pair, and
common residues (the aliphatic group) grew large matched components. q60, exact:
Sens@1FP 0.4680 → 0.3537 (blosum62), 0.3982 (group), 0.3097 (size), at 37–60 s per query.

**Selection** (q250; Δ vs exact search on the same queries):

| scoring | mutant `:*` | mutant `--aa-subst blosum62` | exact `--aa-subst blosum62`, Δ | exact `:*`, Δ |
| --- | --- | --- | --- | --- |
| aa5542d (full IDF) | 0.4241 | 0.3703¹ | −0.1143¹ | — |
| weight 1.0 + cap + edge + matched | **0.4237** | 0.3968 | −0.0198 (−0.034, −0.007) | −0.0050 |
| weight 0.75 + edge + matched | 0.4149 | **0.4021** | −0.0094 (−0.022, +0.002) | −0.0045 |
| **weight 0.75 + cap + edge + matched** (shipped) | 0.4123 | 0.4015 | **−0.0017** (−0.011, +0.008) | **−0.0034** |

¹ q60. *cap*: substituted IDF ≤ the query pair's; *edge*: one substituted hit per query edge
and none with an exact hit; *matched*: per-match IDF counts substituted edges only between
matched residues. On q60, weight 0.5 without *matched* kept exact at −0.0055 but recovered
only 0.3879 on mutants; weight 0.25 was worse on both; weight 1.0 with *edge* alone left exact
at −0.085.

Shipped scoring against exact search on mutants: `:*` +0.1273 (159 wins / 11 losses),
`--aa-subst blosum62` +0.1164 (158 / 18); mean 3.4 s and 5.7 s per query vs 2.8 s exact
and 10.6 s for aa5542d `:*`.

Matching above 200 query hashes now prefilters target residues by amino acid code; outputs
were byte-identical on q60 (blosum62) and runtime fell 17.1 → 8.3 s per query.

**Neighbouring bins** (5ae7321). Geometry-only expanded hashes were scored as exact: own IDF,
every hit summed, all component edges in per-match IDF. Variants (q250; motif = human index):

| variant | M-CSA default | M-CSA `--sensitive` | motif |
| --- | --- | --- | --- |
| unchanged | 0.4465 | 0.4536 | — |
| like substitutions (cap, one hit per edge, matched-only) | 0.4584 | 0.4659 | A1 Sens@50FP 0.895 → 0.369, B1 Sens@10FP 0.720 → 0.189 |
| same, weight 0.75 | 0.4516 | 0.4650 | similar losses |
| IDF cap in candidate scoring only | 0.4464 | 0.4536 | — |
| **matched-only in per-match IDF** (shipped) | **0.4573** | **0.4647** | byte-identical |

Capping hits per edge removes the advantage of repeat proteins (C2H2 arrays), which are the
zinc-finger answers. With substitution the shipped variant also helps: mutant `:*` 0.4123 →
0.4216, mutant blosum62 0.4015 → 0.4116, exact blosum62 0.4447 → 0.4561.
Data: `fd-branchbench/selection/geometric_sweep.txt`.

**Shipped (5ae7321), q250** (Sens@1FP; seconds per query):

| flags | exact | Δ vs master (win/loss) | mutant | runtime |
| --- | --- | --- | --- | --- |
| master | 0.4475 | — | — | 5.34 |
| (default) | 0.4591 | +0.0116 (27/10) | 0.2936 | 2.67 |
| `--sensitive` | 0.4690 | +0.0216 (107/35) | — | 3.04 |
| `:*` on the varying residue | — | — | 0.4216 | 3.40 |
| `--aa-subst blosum62` | 0.4575 | +0.0100 (61/42) | 0.4119 | 5.72 |
| `--aa-subst group` | 0.4547 | +0.0073 (58/42) | — | 6.58 |
| `--aa-subst size` | 0.4532 | +0.0058 (56/57) | — | 6.80 |
| `--sensitive --aa-subst blosum62` | 0.4670 | +0.0196 (107/56) | — | 5.52 |

- Queries with no true positive: 10 (master) → 3. Motif commands: default byte-identical to master.
- Answers ranked before the first FP: conservative single-substitution 473 → 560 with
  blosum62, exact 5,875 → 5,606.
- Motif commands: substitution keeps default Sens@1FP on C1/C2 (0.2339; aa5542d 0.0000);
  D1/D2 F1 0.66–0.84 (aa5542d 0.005–0.34); D2 runtime 51–77 s → 17–23 s. Set F1 on the
  prefilter commands still falls (C1 0.883 → 0.370), since more structures pass `--covered-node`.
- Use substitution when the residues may differ; on exact-residue answer sets it is neutral
  at best.

**Index-time expansion** (M-CSA motif-only index, 24,762 motifs; site queries within 12 Å;
5ae7321): index-side expansion is not a substitute for query-side expansion. Radius 1 on the
index vs on the query: Sens@1FP −0.025 [−0.040, −0.011], returned motifs overlap (Jaccard) 0.58;
the index grows 3.5× (r1), 6× (r2), 9× (blosum62), 32× (both). An expanded index is looked up
exactly, so candidate scoring cannot tell its substituted entries from exact ones. On mutant
sites, `:*` on a plain index gains +0.101 while a blosum62 index loses −0.031 (idf ranking).
Details: `fd-branchbench/index_expansion/README.md`.

## 15. Default sort order and `--confident`

Selected from one run per benchmark, re-sorted and re-filtered offline (a build that skips
sorting and `--top` truncation prints every match with every metric column; re-sorting it
reproduced the shipped ranking on 250/250 M-CSA queries). Data:
`fd-branchbench/defaults/`.

### 15.1 Sort order

172 lexicographic orders over 10 per-match keys (1-2 keys, plus 3 keys when the first is a
count) and 9 over the structure keys, scored by Sens@1FP and average precision (AP) on M-CSA
q250 (default, `--sensitive`, mutant `:*`, `--aa-subst blosum62`) and the four motif queries.

| per-match order | M-CSA Sens@1FP | M-CSA AP | motif Sens@1FP |
| --- | --- | --- | --- |
| `idf,rmsd` (2.x default) | 0.4518 | 0.5951 | 0.1817 |
| best plain alternative (`idf,tm_score`) | 0.4523 | 0.5951 | 0.1837 |
| `node_count,idf,rmsd` | 0.4078 | 0.5645 | 0.1958 |
| **`coverage_idf,rmsd`** (shipped) | **0.4583** | **0.6024** | **0.2172** |

`coverage_idf` is `idf` times the fraction of query residues the match covers. No lexicographic
order beat `idf,rmsd` by more than 0.0005 — putting `node_count` first trades M-CSA for the motif
sets — while the product wins on both. Exponents 0.5 / 2 / 3 on the coverage term scored 0.4559 /
0.4514 / 0.4445; `idf` divided by node count was much worse (0.3666). `drmsd` as the second key
was slightly worse than `rmsd` (0.4502 vs 0.4518 with `idf` first).

**Held out** (the 492 runnable q755 queries outside q250, paired, 95% bootstrap CI): AP rises on
all three configs — default +0.0095 [+0.006, +0.013] (262 wins / 137 losses), `--sensitive`
+0.0046 [+0.001, +0.008], mutant `:*` +0.0100 [+0.007, +0.013] — while Sens@1FP is unchanged
within noise (+0.0015, −0.0029, +0.0002; all CIs span 0). The Sens@1FP gain on the selection set
did not reproduce, so the shipped claim is a better overall ranking, not a better first-FP cut.

Per structure, plain coverage wins and the product does not:

| per-structure order | M-CSA Sens@1FP | M-CSA AP | motif Sens@1FP |
| --- | --- | --- | --- |
| `idf,min_rmsd` (2.x default) | 0.3535 | 0.5237 | 0.0913 |
| `idf*coverage^2,min_rmsd` | 0.4045 | 0.5782 | 0.1192 |
| **`max_node_count,min_rmsd`** (shipped) | **0.4301** | **0.5920** | **0.1671** |

### 15.2 `--confident`

Grid over coverage ratio × RMSD × dRMSD, scored as the precision and recall of the returned
structure set (M-CSA q250, per match):

| filter | precision | recall | F1 | queries with a hit | median hits |
| --- | --- | --- | --- | --- | --- |
| none (default output) | 0.055 | 0.741 | 0.081 | 1.00 | 1,804 |
| coverage ≥ 0.8 only | 0.407 | 0.568 | 0.344 | 0.94 | 94 |
| coverage = 1.0, RMSD ≤ 1.0 | 0.795 | 0.405 | 0.421 | 0.88 | 13 |
| **coverage ≥ 0.8, RMSD ≤ 1.0** (shipped) | **0.756** | **0.478** | 0.471 | 0.92 | 21 |
| coverage ≥ 0.8, RMSD ≤ 0.75 | 0.823 | 0.466 | 0.489 | 0.92 | 16 |

- Coverage 0.8 keeps every residue of a 3-4 residue motif and allows one missing residue from
  five up; it beats exact-full coverage on F1 and leaves more queries with an answer. Ratios
  0.65 / 0.70 / 0.75 / 0.80 / 0.90 score F1 0.372 / 0.439 / 0.435 / **0.474** / 0.431 on M-CSA
  (precision 0.50 / 0.61 / 0.61 / 0.76 / 0.80), so 0.8 is the best single value there. On the
  four motif queries a looser 0.70 scores higher (F1 0.898 vs 0.827) because the 23-residue
  query recovers recall; `--confident --max-node-ratio 0.65` reproduces the published protocol
  for such segment queries.
- A dRMSD cap adds nothing once RMSD is capped (identical rows across `drmsd` 0.75-∞).
- RMSD 0.75 scores marginally better than 1.0 on M-CSA; 1.0 ships because it keeps more answers
  (recall 0.478 vs 0.466) and matches the threshold the published protocol uses.
- Zinc/serine motif queries: precision 0.332 → 0.943, recall 0.966 → 0.687.
- Both cutoffs are needed up to 12 residues: coverage alone gives precision 0.33 (≤4 residues)
  to 0.61 (9-12), and the RMSD cap raises those to 0.64-0.95. Past 12 residues coverage alone is
  already precise (M-CSA 13-21 residues: 0.991 at recall 0.627, versus 1.000 at recall 0.319 with
  the cap and one of three queries emptied), and on the 23-residue two-segment zinc query true
  hits are assembled across neighbouring fingers at 8-9 Å RMSD, so the cap empties the list
  (precision 0.996, recall 0.356 with coverage only). Hence `CONFIDENT_RMSD_MAX_RESIDUES = 12`,
  supported by four queries only.
- With `--skip-match` only hash coverage applies, and that alone is weak (precision 0.25 on
  M-CSA, 0.71 on the motif set).

**Held out** (492 M-CSA queries outside q250): coverage ≥ 0.8 with RMSD ≤ 1.0 gives precision
0.789-0.810, recall 0.464-0.474, median 17-18 hits and 92.3-92.8% of queries with a hit across
default, `--sensitive` and mutant `:*` — matching the selection set (0.756 / 0.478 / 21 / 92.4%).
Unfiltered on the same queries: precision 0.045-0.053, median 1,679-1,802 hits.
