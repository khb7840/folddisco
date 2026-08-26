# Feature evaluation: what shipped, what it measured, what was rejected

This branch adds three things to Folddisco:

1. **Non-rigid motif search** — joint bin expansion, `--nonrigid` / `--expand-radius`
2. **A binary lookup cache**, automatic, no flag
3. **A novelty verdict mode** for `query` — `--novelty-mode`

plus superposition-free deformation metrics (`drmsd`, `min_drmsd`) for filtering,
sorting and output, and a fix for a residue-range syntax that used to panic.

This document is the measurement record. Everything in the first half was measured with
the **repository author's own benchmark protocol** — their commands, their answer sets
and their two metrics, Sens@1FP and F1. Where an earlier round of this document reported
average precision, those tables are gone: not because AP is wrong, but because it is not
the metric this project uses, and a measurement nobody on the project reads is not
evidence.

Two conclusions from the earlier rounds were overturned by that protocol and the code
changed with them: the rare-hash IDF filter is **gone** (§7), and `--nonrigid` turns out
to need the author's `--max-node` filter to pay off at all (§2.3).

## 1. How the measurements were made

- **Index**: human proteome, 23,391 AlphaFold structures (`index/h_sapiens_folddisco`,
  hash type `PDBTrRosetta`, default bins) → 20,504 distinct UniProt accessions.
  Download: `aria2c https://opendata.mmseqs.org/folddisco/h_sapiens_folddisco.tar.lz4`
- **The author's eight commands.** A/B/D are scored against their zinc-finger answer set
  (761 accessions), C against their MEROPS S01 serine-peptidase set:

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

- **Metrics.** `--fp 1` walks the ranked list and stops after the first false positive, so
  field 10 is TP@1FP and field 15 is **Sens@1FP** = TP@1FP / answer_len. Without `--fp`,
  all four cells are set operations over the whole list, so field 17 is a **set F1** with
  no dependence on rank order. **F1 is the decision metric here and Sens@1FP is reported
  beside it**; §5 gives the three reasons, and the one flaw in F1 that goes with them.
- **Answer sets.** The zinc set holds 761 accessions. The serine set has 130 lines but
  `answer_len` prints **124**: six proteins carry two MEROPS S01 identifiers (P20231,
  P48740, Q15661, Q2TV78, Q5K4E3, Q7RTY7) and the harness deduplicates. A reader
  recomputing these numbers will not match without that. `--afdb-to-uniprot` throughout,
  so AlphaFold fragments of one protein collapse to one accession — `hits` below is raw
  result lines and is always the larger number.
- **Robustness.** Every delta was re-measured over 200 paired replicates with 5% of the
  answer accessions dropped at random, the same replicate applied to both lists so
  common-mode annotation noise cancels. "Sign consistency" is the fraction of replicates
  keeping the unperturbed sign: 1.00 is robust, 0.5 is noise.
- **Machine**: 20 cores, warm page cache, medians of 13 interleaved repeats for timing.
  Searches are deterministic — repeated runs are byte-identical, so all the variance
  below is answer-set variance.

Reproduce a row:

```bash
IDX=index/h_sapiens_folddisco
folddisco query -i $IDX -p query/1G2F.pdb -q F207,F212,F225,F229 \
  -t 12 --covered-node 3 --max-node 4 --rmsd 1.0 --per-structure > result.tsv

# field 10 = TP@1FP, field 15 = Sens@1FP
folddisco benchmark -r result.tsv -a <answers>.tsv -i $IDX --afdb-to-uniprot --fp 1
# fields 14, 15, 17 = precision, recall, F1 over the whole list
folddisco benchmark -r result.tsv -a <answers>.tsv -i $IDX --afdb-to-uniprot

# the same numbers plus TP@kFP at other depths, from the repository
python3 scripts/eval_metrics.py result.tsv <answers>.tsv $IDX.lookup
```

## 2. `--nonrigid`

`default` below is this branch with no extra flag. It is **byte-identical to master**
(2a756d9) on all eight commands (§3), so one column serves for both.

| cmd | config | hits | TP@1FP | Sens@1FP | precision | recall | F1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **A1** 4-res, prefilter | master = default | 814 | 31 | 0.0407 | 0.9100 | 0.9435 | **0.9265** |
| | `--nonrigid` | 878 | **42** | 0.0552 | 0.8787 | **0.9711** | 0.9226 |
| **A2** 4-res, matched | master = default | 746 | 31 | 0.0407 | 0.9655 | 0.9198 | 0.9421 |
| | `--nonrigid` | 802 | **42** | 0.0552 | 0.9584 | **0.9698** | **0.9641** |
| **B1** 3-res, prefilter | master = default | 767 | **43** | 0.0565 | 0.9489 | 0.9277 | 0.9382 |
| | `--nonrigid` | 826 | 22 | 0.0289 | 0.9328 | **0.9671** | **0.9497** |
| **B2** 3-res, matched | master = default | 756 | **43** | 0.0565 | 0.9591 | 0.9251 | 0.9418 |
| | `--nonrigid` | 813 | 22 | 0.0289 | 0.9485 | **0.9671** | **0.9577** |
| **C1** Ser-His-Asp, prefilter | master = default | 107 | 29 | 0.2339 | 0.9533 | 0.8226 | 0.8831 |
| | `--nonrigid` | 114 | **32** | 0.2581 | **0.9561** | **0.8790** | **0.9160** |
| **C2** Ser-His-Asp, matched | master = default | 113 | 29 | 0.2339 | 0.9640 | 0.8629 | 0.9106 |
| | `--nonrigid` | 113 | **32** | 0.2581 | 0.9640 | 0.8629 | 0.9106 |
| **D1** 23-res segments, prefilter | master = default | 800 | **26** | 0.0342 | 0.8819 | 0.8830 | **0.8825** |
| | `--nonrigid` | 800 | 25 | 0.0329 | 0.8729 | 0.8752 | 0.8740 |
| **D2** 23-res segments, matched | master = default | 697 | **26** | 0.0342 | 0.9722 | 0.8739 | **0.9204** |
| | `--nonrigid` | 705 | 25 | 0.0329 | 0.9594 | 0.8686 | 0.9117 |

`--nonrigid` always raises deep recall (+0.026 to +0.056) and always lowers precision.
Whether that trade is a win depends on the command:

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

Every F1 delta is robust (1.00). Not one Sens@1FP delta is (0.50–0.85) — including the
two largest, A1's +11 and B1's −21.

### 2.1 It is a short-motif feature and it needs `--max-node`

Two conditions, both measured, both worth stating before anyone turns the flag on:

**Use it with `--max-node <n_residues>`.** The same 4-residue query, adding the author's
filters one at a time, F1 (whole list, 761-accession zinc set):

| variant | default | `--nonrigid` | Δ |
| --- | --- | --- | --- |
| `--skip-match`, no filters | 0.6511 | 0.5979 | **−0.0532** |
| `--covered-node 3 --skip-match` | 0.9265 | 0.9226 | −0.0039 |
| `--covered-node 3 --rmsd 1.0` | 0.9312 | 0.9360 | +0.0048 |
| `--covered-node 3 --max-node 4` | 0.9397 | 0.9609 | **+0.0212** |
| `--covered-node 3 --max-node 4 --rmsd 1.0` (author's A2) | 0.9421 | 0.9641 | **+0.0220** |

`--covered-node 3` is the dominant filter of the whole protocol — dropping it costs the
default search 0.28 F1, so any evaluation run without it is measuring a different search.
And it is `--max-node`, not `--rmsd`, that converts the extra recall into F1: requiring
the full motif to be covered within one structure discards the low-coverage material bin
expansion pulls in while keeping the fully-covered true positives it finds. An earlier
round of this document measured `--nonrigid` with none of these filters — that is,
in the one configuration where it is at its worst (−0.053 F1).

**Do not use it for long segment queries.** On the author's 23-residue segment query it
loses on both metrics and both commands: F1 −0.0085 (D1) and −0.0087 (D2) at sign
consistency 1.00, TP@1FP 26 → 25. Long queries are already saturated (recall 0.87–0.88
before the flag) and the extra candidates only dilute the list. This reproduces, on the
author's own query, the regression an earlier round found on a 16-residue motif.

### 2.2 Runtime

Interleaved, 13 timed repeats after 2 warm-ups, whole-process wall clock, `-t 12`:

| config | A1 prefilter | A2 matched |
| --- | --- | --- |
| master | 9.1 ms (8.5–9.8) | 158.9 ms (154–177) |
| default | 9.1 ms (8.6–10.3) | 157.4 ms (148–174) |
| `--nonrigid` | 9.9 ms (8.6–10.4) | 174.7 ms (168–183) |

`--nonrigid` costs **+11% on the matched command** (157 → 175 ms). Its prefilter
difference, 0.8 ms, is inside the 1.7 ms spread of a single configuration and should not
be quoted as a cost. `--covered-node 3` shrinks the candidate set before matching, which
is why this is cheaper than the +16% an earlier round measured without it.

## 3. The default path: byte-identical to master, and one real repair

With the rare-hash filter removed (§7) the branch default is **byte-identical (`cmp`) to
master on all eight of the author's commands** — A1 814 lines, A2 746, B1 767, B2 756,
C1 107, C2 113, D1 800, D2 697. So the `expand.rs` rewrite, `MAX_HASHES_PER_PAIR`, the
substitution/expansion composition and the angle-domain fix-up are a provable no-op at
the author's tolerances, and the small default-path regression an earlier round of this
document reported was the filter, not this code.

They are not a no-op in general. Widening the tolerance on the A1 query:

| tolerance | master vs default | master TP@1FP / recall / F1 | default TP@1FP / recall / F1 |
| --- | --- | --- | --- |
| `-d 0.5 -a 5` | byte-identical | 31 / 0.9435 / 0.9265 | 31 / 0.9435 / 0.9265 |
| `-d 2.0 -a 20` | **differs** | 22 / 0.8476 / 0.8764 | **29 / 0.9290 / 0.8927** |

At a tolerance wider than one bin, master queries only the two extreme bins and skips
those in between — visible as recall *falling* as the tolerance widens. The branch
sub-steps the offsets so the covered bins are contiguous, which repairs it. That is the
standing justification for the default-path changes: they are a correctness fix that the
author's protocol cannot see because it uses the default tolerances.

The other corrections in the same module, all active at every setting:

- **Angles are pulled back into their domain.** Torsions wrap at ±π, `acos`-derived
  angles reflect at their ends. Previously a wide `-a` pushed the bin index past the end
  of its bit field on hash types that discretize radians directly, and on sin-cos types
  it asked for a `sin` sign no real structure can produce. Sensitivity used to *fall*
  past `-a 30` (12 → 10 → 4 → 2 matches at `-a 30/45/60/90`); it is now monotone.
- **Substitutions compose with tolerance.** `164:H` used to be hashed only at the
  observed geometry, so the alternative residue got no distance or angle slack.

Distances are deliberately *not* clamped. A structure can hold a Cβ–Cβ distance beyond
`MAX_DIST` — only the Cα distance is checked against the cutoff — so the index itself
contains those out-of-window encodings. Indexing and querying share one discretizer, so
mirroring it exactly is what makes a perturbed query hash mean the same thing as an index
hash.

## 4. Binary lookup cache

Loading the lookup file dominates index startup at AFDB scale. The parsed result is now
cached in a binary sidecar next to the lookup, written on the first text parse and used by
every load after that. There is no flag: it is on, and it falls back to text parsing
whenever anything is wrong. It does not touch ranking, so the author's protocol has
nothing to say about it and the measurements below stand from the previous round.

The cache is laid out for **parallel** decoding — a 44-byte header, then a fixed-stride
40-byte record table, then a names blob — so the reader mmaps it and decodes records with
rayon, matching the existing mmap + `par_lines` text path. A serial reader over
variable-length records was measured *slower* than text parsing at 8 threads, which would
have defeated the purpose.

Medians of 5, warm page cache; `text` is the pre-cache path with neither cache read nor
cache write:

| entries | file size | threads | text | cache | speedup |
| --- | --- | --- | --- | --- | --- |
| 23,391 (human index) | 1.19 MB | 1 | 5.28 ms | 1.23 ms | **4.3x** |
| | | 8 | 1.23 ms | 0.97 ms | 1.27x |
| | | 20 | 1.94 ms | 1.87 ms | 1.04x |
| 10^5 | 4.9 MB | 1 | 13.16 ms | 3.41 ms | **3.9x** |
| | | 8 | 3.70 ms | 2.19 ms | 1.69x |
| | | 20 | 3.35 ms | 2.84 ms | 1.18x |
| 10^6 | 50.8 MB | 1 | 145.9 ms | 67.7 ms | **2.2x** |
| | | 8 | 40.7 ms | 31.4 ms | 1.30x |
| | | 20 | 33.1 ms | 30.6 ms | 1.08x |
| 10^7 | 528 MB | 1 | 1498.6 ms | 666.8 ms | **2.25x** |
| | | 8 | 346.8 ms | 262.6 ms | 1.32x |
| | | 20 | 277.7 ms | 235.7 ms | 1.18x |

Never slower at any size or thread count. The speedup shrinks as threads rise because
both paths allocate one `String` per entry and the cache only removes the integer and
float parsing; at 20 threads allocation dominates. The next win at AFDB scale would be
interning names or returning slices into the mmap, not a different binary layout.

The first load pays for the write — roughly one extra text parse (+0.63 ms at 23,391
entries and 8 threads, +38 ms at 10^6, +312 ms at 10^7) — so a fresh index pays the cache
back after about **4 loads at 8 threads**, 9 at 20. The cache file is 1.23–1.33x the size
of the text lookup.

### What invalidates it, and the one case it cannot see

The header carries the source lookup's **length**, its **mtime in nanoseconds** and the
**length of the names blob**; all three are validated on every load, and a record with an
empty name is rejected. That covers every realistic edit — including the coarse-timestamp
collisions an mtime comparison alone gets wrong (a 1-second-granularity filesystem, or a
rebuild inside one tick, previously served the *old* records with every field wrong and no
signal), a corrupted entry count that would otherwise decode as an empty index, and the
all-zero region a torn write leaves.

**It cannot see a same-length edit whose mtime is forged back to its original value.**
That is by design: there is no content hash, because hashing hundreds of megabytes would
cost more than the load the cache exists to speed up. Delete the `.cache` sidecar if you
have rewritten a lookup that way.

Concurrent first-loads are safe without a temp-and-rename: cache content is a pure
function of the lookup, so two writers emit identical bytes, and a reader that catches a
write in progress sees a short file and rejects it (tested: two threads, 40 trials, 0
corrupt, 0 validated-but-wrong). That argument holds only while the bytes stay
deterministic, which is recorded in the code.

## 5. Which metric governs, and what each one is blind to

**F1 (no `--fp`) decides; Sens@1FP is reported.** Three reasons, all from this round:

1. **Reproducibility.** F1 deltas hold their sign in 1.00 of dropout replicates
   throughout §2; Sens@1FP deltas hold theirs in 0.50–0.85. Every Sens@1FP conclusion in
   the headline table is inside its own replicate spread.
2. **At k=1 it is one accession.** The 31-vs-27 gap that an earlier build showed on A1 was
   a single accession (Q8NB15) moving from rank 32 to rank 28 in a 788-accession list —
   and that accession is a *true positive* in the repository's other zinc answer set, so
   the metric was partly reading annotation coverage.
3. **It is blind to the author's own filters.** TP@1FP is identical across every
   combination of `--covered-node`, `--max-node`, `--rmsd` and `--skip-match` — 31 for
   the default and 42 for `--nonrigid` in all five variants of §2.1 — because those
   filters only remove structures ranked below the first false positive. A metric that
   cannot see `--max-node 4` cannot be used to tune a protocol built around it.

**And the flaw in F1, which this round exposed: without `--fp` it is a pure set metric
with no rank dependence at all.** C2 is the proof. `--nonrigid` returns *the same 111
accessions* as the default there, verified, in a different order — promoting genuine
proteases (kallikrein 7 P49862, neutrophil elastase P08246, HGF activator Q9P0G3 among
them) toward the top. F1 does not move by 0.0001 (0.9106 both ways); only Sens@1FP sees it
(29 → 32). If ranking quality matters, F1 needs a rank-aware companion, and Sens@1FP at
k=1 is too fragile to be it: **`--fp 10` or `--fp 100` is the better one** (an earlier
round measured sign consistency 1.00 for TP@100FP against 0.28–0.49 for TP@5FP).

Where the two metrics disagree, they disagree because they measure different things:

* `--nonrigid` on A1: Sens@1FP better (31 → 42), F1 worse (−0.0039). The list is worse.
* `--nonrigid` on B1/B2: Sens@1FP much worse (43 → 22), F1 better (+0.0115/+0.0159). One
  extra early false positive buys ~29 extra deep true positives.

## 6. The answer set can flip both metrics on its own

Same ranked lists, scored against the repository's older `data/zinc_answer.tsv` (1817
accessions, of which the author's 761 are essentially a subset) instead of the author's set:

| cmd | config | 761-set TP@1FP / F1 | 1817-set TP@1FP / F1 |
| --- | --- | --- | --- |
| A1 | default | 31 / 0.9265 | 90 / 0.5902 |
| A1 | `--nonrigid` | **42** / 0.9226 | 61 / **0.6005** |
| A2 | default | 31 / 0.9421 | 89 / 0.5673 |
| A2 | `--nonrigid` | **42** / **0.9641** | 61 / **0.5922** |

On A1 **both metrics reverse their verdict from the answer set alone**, with the ranked
lists untouched: `--nonrigid` wins Sens@1FP and loses F1 on the 761-accession set, and
loses Sens@1FP and wins F1 on the 1817-accession set. The reason is what each set limits:
these searches reach recall 0.87–0.97 against the 761 set, so F1 there is
precision-limited and the harder test; against the 1817 set recall is 0.40–0.44, so F1 is
recall-limited and rewards anything that adds recall — which is exactly what `--nonrigid`
does.

**Absolute F1 values are not comparable between answer sets, and neither is the sign of a
feature comparison.** Every number in this repository's documentation names its answer set
for that reason.

## 7. Alternatives measured and rejected

| feature | verdict |
| --- | --- |
| **Rare-hash IDF filter** (`EXPANDED_HASH_IDF_MAX_EXCESS`) | **Rejected and removed.** Dropped a tolerance-expanded hash whose IDF exceeded its observed parent's by more than a margin. Swept at margins 0/2/3/5/8/off across the author's eight commands at both radii: the shipped 5.0 was **worse-or-equal in 16 of 16 cells on both metrics and better in none**. F1 harm 0.0001–0.0020 with sign consistency 0.91–1.00; margins 8.0 and off indistinguishable, so no margin beat no filter; no runtime to trade against it (158.7 vs 159.8 ms matched). It was also the entire difference between this branch's default path and master. Preserved on branch `rare-hash-idf-filter`. |
| **Torsion-angle ENM query sampling (`--enm-sample`)** | **Rejected and removed.** Against `--nonrigid`: −0.047 AP on an independent serine family with 0/200 replicates positive, −0.004 on a 3-residue zinc motif, against +0.019 on the 4-residue motif it was developed on; 3.1x prefilter runtime; 17 crates (`nalgebra`) for one eigendecomposition. Its conformer count could not be chosen by measurement — 2 through 40 conformers produced byte-identical output. Three real defects in it were found and fixed first (rotations crossing chain boundaries, a `--nma-rmsd` whose Ångström units were fictional, trivial modes crowding out real ones); the fixes took a cross-chain motif from 0.778 to 0.850 AP and still lost. Preserved on branch `enm-torsion-sampling`. |
| ANM/NMA query sampling | Rejected. Never beat torsion-ENM on quality and ran at 1.4 s against a 9 ms baseline (~110x), most likely a full 3N×3N eigendecomposition where torsion space is far smaller. |
| DMS / PAS / SOS composite ranking metrics | Rejected. With an identical 1635-structure candidate set none beats plain RMSD as the tiebreaker after `node_count`. DMS is a normalised dRMSD, and that does not improve ranking either. |
| `--dist-ratio` (distance-proportional tolerance) | Rejected and removed. Raised total recall while losing true positives at every early-precision point: widening the tolerance most for long pairs adds exactly the least specific matches. |
| An absolute IDF ceiling (`--max-idf`) instead of the relative filter | Not implemented, and now moot: the relative filter it would have replaced earns nothing at all (above). |
| Perturbation analytics (`analyze-perturb`) | Not a search feature — byte-identical search output. Judge on the diagnostics alone. |

The `--enm-sample`, ANM, DMS/PAS/SOS and `--dist-ratio` rows were measured in an earlier
round against average precision on this repository's own answer sets, not under the
author's protocol. Their verdicts are unaffected by the metric change — none of them was a
close call, and three of them are also slower — but the AP figures quoted there are not
comparable with the F1 figures above.

## 8. Deformation metrics

`drmsd` (distance-matrix RMSD) and `max_dist_deviation` compare the internal distances of
a match instead of superposing it, so a motif whose halves swung apart on a hinge keeps a
low dRMSD where its RMSD is large. Available to `--sort-by`, `--format-output` and
`--drmsd`, plus `min_drmsd` per structure. They are read straight off the matched
coordinates, so no superposition is involved.

On a synthetic hinge (chains C and G of `4CHA` rotated about the chain-C centroid, motif
`B55-58,C193-196` spanning the hinge, searched against an index of the unrotated
structure) there is a window where the sensitivity work decides the result: at a 20° hinge
the default search collapses to a 3-residue mis-assignment at 2.81 Å while `--nonrigid`
recovers all 8 residues at 1.31 Å. Below 10° the default already suffices; by 30°
radius-2 expansion is not enough either.

Ranking honestly did not improve: sorting by dRMSD does not beat sorting by RMSD after
`node_count`. These are filtering and reporting metrics, and the help text says so.

## 9. Novelty verdict mode

Replaces the result listing with one line per query: `query_id, verdict, best hit,
coverage, RMSD, query residues`. KNOWN requires both coverage ≥ `--novelty-coverage`
(default 0.8) **and** a best-hit RMSD ≤ `--novelty-rmsd` (default 2.0 Å); covered but
failing either is PARTIAL_MATCH; nothing covered is NOVEL. Like the cache, it does not
touch ranking, so these numbers stand from the previous round.

The RMSD gate exists because of a measurement. Over 100 randomly drawn, spatially compact,
indexable motifs (sizes 3–8, five source structures), **96 came back KNOWN and 32 of those
96 had a best hit worse than 2.0 Å** — median 0.85 Å, p90 5.41 Å, max 10.4 Å. An
8-residue motif was called KNOWN at coverage 1.0 with a 6.15 Å best hit. The verdict was
ignoring the RMSD printed on its own line. A 2.0 Å gate separates the measured cases
cleanly: the chymotrypsin triad (best hit 0.06 Å, correctly identifying human chymotrypsin
C) and the zinc finger (0.17 Å) stay KNOWN even at 0.5 Å, while the 4.3–6.1 Å lookalikes
drop to PARTIAL_MATCH.

Two honest limits:

- **Under `--skip-match` no RMSD is computed**, so the verdict falls back to coverage
  alone and is correspondingly weaker.
- **A query whose residues are further apart than the index distance cutoff produces no
  hashes**, and therefore no hits. That is reported as `NO_HASHES` with a warning on
  stderr, not as NOVEL — an unsearchable query is not a discovery.

`--nonrigid` composes with the verdict, but on this evidence the effect is theoretical:
across 200 comparisons it changed exactly one verdict, and coverage rose twice and never
fell.

## 10. A third answer set, kept for independence

`data/serine_answer.tsv` (89 accessions) and `scripts/build_serine_answer.py` are **not**
the author's serine set. They are an independent S1-clan set built from two PROSITE
sequence patterns (PS00134 His, PS00135 Ser, His before Ser), defined by sequence and so
never circular with anything Folddisco computes. They are kept because a third,
independently derived answer set is worth having, and because the generator matters more
than the file: the set lists entries of `index/h_sapiens` and must be regenerated for any
other index. It agreed with the author's MEROPS set on the direction of the `--nonrigid`
result.

The author's `serinepeptidase_answer.tsv` (130 lines / 124 accessions, MEROPS S01) and
`zincfinger_answer.tsv` (761 accessions) are the sets of record for everything in §2–§6
and are not committed here.

## 11. Known issues found but not fixed

- **`CompactStructure::build` never resets its `c` variable between residues**, so the Cβ
  approximation for a residue lacking Cβ can fall back on an *earlier residue's* carbonyl
  carbon. Fixing it changes every hash and would invalidate published indices, so the
  existing behaviour is preserved verbatim. Worth addressing on a future index-format
  change.
- **`FoldcompDbReader::new` sorts its lookup by database key**, while
  `get_foldcomp_db_entry_by_name` binary-searches the same lookup by *name*, so
  `read_single_structure(name)` fails for most entries unless `sort_lookup_by_name()` is
  called first. `examples/seqdump.rs` calls it; other callers may not.
- **A 3-column foldcomp-native lookup panics `load_lookup_from_file`**, which indexes
  columns 2 and 3 unconditionally (`data/foldcomp/example_db.lookup` is such a file).
  Pre-existing, and it affects the text and cache paths equally.

A fourth was found and **fixed**: `-q F204-F215` — a residue range whose end repeats the
chain letter, the form used in the author's own benchmark commands — panicked with
`Invalid end residue: ParseIntError`, on this branch and on master. It now parses, a
mismatched chain is rejected with a message, and a malformed `-q` exits with a diagnostic
instead of a backtrace.

## 12. Limits of these measurements

- **One index (23,391 human structures) and two answer sets** for the protocol above. An
  886-query M-CSA benchmark on the PDB index is a separate exercise and is not covered
  here.
- **The fragility analysis holds the ranked lists fixed and perturbs the answer set.** It
  models annotation noise, not index or query variation.
- **The `-d 2.0 -a 20` comparison in §3 is outside the author's protocol.** It is included
  only to bound the scope of the byte-identity claim.
- **Sections 4, 8, 9 and the last four rows of §7 predate the author's protocol** and were
  measured against average precision on this repository's own answer sets. None of them is
  a ranking-quality claim that the metric change touches, but their numbers are not
  comparable with §2's.
- **Single machine, warm cache.** 20 cores, local NVMe.
