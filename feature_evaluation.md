# Feature evaluation: what shipped, what it measured, what was rejected

This branch adds four things to Folddisco:

1. **Non-rigid motif search** — joint bin expansion, `--nonrigid` / `--expand-radius`
2. **A rare-hash floor** on tolerance-expanded hashes, always on
3. **A binary lookup cache**, automatic, no flag
4. **A novelty verdict mode** for `query` — `--novelty-mode`

plus superposition-free deformation metrics (`drmsd`, `min_drmsd`) for filtering,
sorting and output.

This document is the measurement record: what each change is, what it measured, and
what was measured and rejected along the way. Every number here was produced against
the shipped build; where a claim could not survive measurement it has been removed
rather than softened.

## How the measurements were made

- **Index**: human proteome, 23,391 AlphaFold structures (`index/h_sapiens_folddisco`,
  hash type `PDBTrRosetta`, default bins). Download:
  `aria2c https://opendata.mmseqs.org/folddisco/h_sapiens_folddisco.tar.lz4`
- **Two motif families**:
  - zinc finger from `query/1G2F.pdb`, three residue selections — `F207,F212,F225,F229`
    (4 residues), `F207,F212,F225` (3), `F205-212,F223-230` (16). Answers:
    `data/zinc_answer.tsv`, 1,817 UniProt accessions.
  - serine-peptidase catalytic triad from `query/4CHA.pdb`, `B57,B102,C195`. Answers:
    `data/serine_answer.tsv`, 89 human S1-clan peptidases, built from two PROSITE
    sequence patterns by `scripts/build_serine_answer.py` — defined independently of
    anything Folddisco computes, so it cannot be circular.
- **Metrics**: **average precision** over the ranked list, **TP@100FP** (true positives
  found while walking the ranked list until the 100th false positive) and **recall**.
  Matching is on UniProt accessions after `--afdb-to-uniprot`; the 23,391 index entries
  dedupe to 20,504 accessions, so "hits" (raw result lines, AFDB fragments) is always
  larger than the number of accessions scored.
- **Robustness**: each delta was re-measured over 200 paired replicates with 5% of the
  answer set dropped at random. "Sign consistency" is the fraction of replicates in
  which the delta keeps its unperturbed sign: 1.00 or 0.00 means robust, near 0.5 means
  noise.
- **Machine**: 8 threads, warm page cache, medians of 13 interleaved repeats for timing.
  All searches are deterministic — repeated runs are byte-identical.

**TP@5FP and TP@10FP are not used as evidence here.** They are the rank of the *k*-th
false positive minus *k*, so a single false positive changing position moves them by
more than 100, and their sign flips in 50–70% of dropout replicates (sign consistency
0.34 on the 4-residue zinc motif). Earlier versions of this document led with them; the
tables below do not report them at all.

Reproduce a row:

```bash
IDX=index/h_sapiens_folddisco
folddisco query -i $IDX -p query/1G2F.pdb -q F207,F212,F225,F229 \
  -t 8 --skip-match --per-structure --format-output tid > result.tsv

# every column of the tables below, AP included
python3 scripts/eval_metrics.py result.tsv data/zinc_answer.tsv $IDX.lookup

# the same TP@kFP and recall from the shipped binary. Note that with --fp the
# recall column is TP@kFP / answer_len, so deep recall needs a second run without it
folddisco benchmark -r result.tsv -a data/zinc_answer.tsv -i $IDX --afdb-to-uniprot --fp 100
folddisco benchmark -r result.tsv -a data/zinc_answer.tsv -i $IDX --afdb-to-uniprot
```

`folddisco benchmark` does not compute average precision, so `scripts/eval_metrics.py`
does: it replicates `--afdb-to-uniprot`'s identifier parsing and deduplication exactly
and adds AP, TP@kFP and precision at fixed depths. AP is the standard one — accumulate
precision@i at every rank holding a true positive, divide by the size of the full answer
set, so answers that never appear contribute zero and two configurations returning
different numbers of hits stay comparable. Six rows of the tables below were recomputed
with the committed script as a check, and all six agree to four decimals.

## 1. Non-rigid search — `--nonrigid`, `--expand-radius`

Motifs are rarely rigid. The same catalytic site in two homologs, or the same site
before and after a conformational change, keeps its residues in the same arrangement
while the distances and angles between them drift. Folddisco discretizes that geometry
into bins, so a pair whose features land on the far side of a bin boundary produces a
different hash and is missed. The expansion used to perturb one feature dimension at a
time, so a target pair whose distance *and* angle both drifted across a boundary was
unreachable. `--expand-radius 2` searches pairs of dimensions together, and `--nonrigid`
is a preset for it.

`master` is 2a756d9; `default` is this branch with no new flag; all rows
`-t 8 --skip-match --per-structure`.

| benchmark | config | hits | AP | TP@100FP | recall |
| --- | --- | --- | --- | --- | --- |
| serine triad, 89 answers | master | 581 | 0.8580 | 88 | 1.0000 |
| | default | 581 | 0.8580 | 88 | 1.0000 |
| | `--nonrigid` | 587 | **0.8947** | 89 | 1.0000 |
| zinc, 4 residues | master | 1645 | 0.4716 | 782 | 0.4959 |
| | default | 1640 | 0.4715 | 782 | 0.4959 |
| | `--nonrigid` | 1901 | **0.4893** | 799 | **0.5190** |
| zinc, 3 residues | master | 595 | 0.1600 | 268 | 0.1970 |
| | default | 591 | 0.1604 | 270 | 0.1970 |
| | `--nonrigid` | 1252 | **0.3866** | 304 | **0.4843** |
| zinc, 16 residues | master | 22491 | 0.5026 | 663 | 0.9923 |
| | default | 22489 | 0.5027 | 664 | 0.9923 |
| | `--nonrigid` | 22776 | *0.4988* | 665 | 0.9934 |

Robustness of `--nonrigid` versus master, 200 paired replicates:

| benchmark | AP Δ | median Δ | p10..p90 | sign consistency |
| --- | --- | --- | --- | --- |
| serine triad | +0.0367 | +0.0415 | +0.0026..+0.0480 | 0.91 |
| zinc 4 res | +0.0177 | +0.0168 | +0.0147..+0.0184 | 1.00 |
| zinc 3 res | +0.2266 | +0.2151 | +0.2116..+0.2180 | 1.00 |
| zinc 16 res | **−0.0038** | −0.0036 | −0.0053..−0.0020 | **0.01** |

**`--nonrigid` is a short-motif feature, and it costs something on long ones.** The
3-residue zinc motif more than doubles its average precision and its recall; the
independent serine family gains 0.037 AP with 0.91 sign consistency, which is the
generality signal that matters most here — it is a family the feature was not developed
against. The 16-residue motif *loses* 0.0038 AP, negative in 99% of replicates: a motif
that long is already at recall 0.99 and the extra candidates only dilute the ranking.
That regression is small but it is real and reproducible, and it belongs next to the
recommendation rather than in a footnote.

On a synthetic hinge (chains C and G of `4CHA` rotated about the chain-C centroid, motif
`B55-58,C193-196` spanning the hinge, searched against an index of the unrotated
structure), there is a window where this decides the result: at a 20° hinge the default
search collapses to a 3-residue mis-assignment at 2.81 Å while `--nonrigid` recovers all
8 residues at 1.31 Å. Below 10° the default already suffices; by 30° radius-2 expansion
is not enough either.

### Runtime

Interleaved medians of 13 repeats, `-t 8`, 4-residue zinc motif:

| config | prefilter (`--skip-match`) | with residue matching |
| --- | --- | --- |
| master | 8.7 ms (7.1–9.0) | 410.9 ms |
| default | 8.7 ms (7.6–9.2) | 411.4 ms |
| `--nonrigid` | 8.6 ms (7.9–9.5) | 478.1 ms (**1.16x**) |

The prefilter spread *within* one configuration is 1.9 ms, larger than any difference
between the configurations, so **there is no prefilter cost to quote for `--nonrigid`**
at this precision. An earlier version of this document claimed 1.04–1.05x; that was
inside the noise. The full-match cost is real and comes from the larger candidate set
that has to be superposed.

### Corrections that ship at radius 1, on correctness grounds

Three further corrections in `src/controller/expand.rs` are active at every setting,
including the default. **They carry no quality claim**: at the shipped rare-hash margin
the default path is byte-identical to master on the serine benchmark and within 0.0005
AP of master on all three zinc motifs. They ship because the old behaviour was wrong,
not because the benchmark improved.

- **Wide tolerances are sub-stepped.** A single jump of 2.5 bins skipped the bins in
  between; offsets are now spaced within one bin so the covered bins are contiguous.
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
mirroring it exactly is what makes a perturbed query hash mean the same thing as an
index hash.

## 2. Rare-hash floor (`src/controller/count_query.rs`)

A tolerance-expanded hash that turns out far rarer in the database than the observed
hash it came from is dropped. A rare hash carries a large IDF, so one spurious hit on it
can outrank several real ones. The threshold is a ratio between two IDFs from the same
index, so it needs no tuning per database size. The shipped margin is 5.0 in log2 units
— "up to 32x rarer than its parent is fine".

**This is a floor, not an improvement, and not a false-positive filter.** Measured in AP
against *no filter at all*, over 200 paired replicates on three motifs:

| margin | effect |
| --- | --- |
| 0 | **−0.0116 AP** on the 4-residue motif, 0/200 replicates positive — robustly harmful |
| 3 | −0.0011 AP on the 16-residue motif, worse in 97% of replicates |
| **5 (shipped)** | never robustly worse on any of the three motifs |
| 2 … off | all within ±0.004 AP of each other |

Margin 0 being measurably harmful is the only robust result in the sweep, and it is the
entire justification for running a filter. 5.0 is then the smallest margin that is never
robustly worse; 3.0, which this branch shipped earlier, was robustly worse on the
16-residue motif. Five of the eight shipped benchmark rows are marginally *lower* at 5.0
than at 3.0 — the retune buys robustness, not quality.

The filter does not reduce false positives: FP in the top 100 is 0–1 at every setting
including off. It is not a cost control either: the prefilter runs at 9.1 ms with it and
9.3 ms without.

## 3. Binary lookup cache (`src/index/lookup.rs`)

Loading the lookup file dominates index startup at AFDB scale. The parsed result is now
cached in a binary sidecar next to the lookup, written on the first text parse and used
by every load after that. There is no flag: it is on, and it falls back to text parsing
whenever anything is wrong.

The cache is laid out for **parallel** decoding — a 44-byte header, then a fixed-stride
40-byte record table, then a names blob — so the reader mmaps it and decodes records
with rayon, matching the existing mmap + `par_lines` text path. A serial reader over
variable-length records was measured *slower* than text parsing at 8 threads, which
would have defeated the purpose.

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

The cache is never slower than text at any size or thread count. The speedup shrinks as
threads rise because both paths allocate one `String` per entry and the cache only
removes the integer and float parsing; at 20 threads allocation dominates. The next win
at AFDB scale would be interning names or returning slices into the mmap, not a
different binary layout.

The first load pays for the write — roughly one extra text parse (+0.63 ms at 23,391
entries and 8 threads, +38 ms at 10^6, +312 ms at 10^7) — so a fresh index pays the
cache back after about **4 loads at 8 threads**, 9 at 20. The cache file is 1.23–1.33x
the size of the text lookup.

### What invalidates it, and the one case it cannot see

The header carries the source lookup's **length**, its **mtime in nanoseconds** and the
**length of the names blob**; all three are validated on every load, and a record with
an empty name is rejected. That covers every realistic edit — including the
coarse-timestamp collisions an mtime comparison alone gets wrong (a 1-second-granularity
filesystem, or a rebuild inside one tick, previously served the *old* records with every
field wrong and no signal), a corrupted entry count that would otherwise decode as an
empty index, and the all-zero region a torn write leaves.

**It cannot see a same-length edit whose mtime is forged back to its original value.**
That is by design: there is no content hash, because hashing hundreds of megabytes would
cost more than the load the cache exists to speed up. Delete the `.cache` sidecar if you
have rewritten a lookup that way.

Concurrent first-loads are safe without a temp-and-rename: cache content is a pure
function of the lookup, so two writers emit identical bytes, and a reader that catches a
write in progress sees a short file and rejects it (tested: two threads, 40 trials, 0
corrupt, 0 validated-but-wrong). That argument holds only while the bytes stay
deterministic, which is recorded in the code.

## 4. Novelty verdict mode — `--novelty-mode`

Replaces the result listing with one line per query:
`query_id, verdict, best hit, coverage, RMSD, query residues`. KNOWN requires both
coverage ≥ `--novelty-coverage` (default 0.8) **and** a best-hit RMSD ≤ `--novelty-rmsd`
(default 2.0 Å); covered but failing either is PARTIAL_MATCH; nothing covered is NOVEL.

The RMSD gate exists because of a measurement. Over 100 randomly drawn, spatially
compact, indexable motifs (sizes 3–8, five source structures), **96 came back KNOWN and
32 of those 96 had a best hit worse than 2.0 Å** — median 0.85 Å, p90 5.41 Å, max
10.4 Å. An 8-residue motif was called KNOWN at coverage 1.0 with a 6.15 Å best hit. The
verdict was ignoring the RMSD printed on its own line. A 2.0 Å gate separates the
measured cases cleanly: the chymotrypsin triad (best hit 0.06 Å, correctly identifying
human chymotrypsin C) and the zinc finger (0.17 Å) stay KNOWN even at 0.5 Å, while the
4.3–6.1 Å lookalikes drop to PARTIAL_MATCH.

Two honest limits:

- **Under `--skip-match` no RMSD is computed**, so the verdict falls back to coverage
  alone and is correspondingly weaker.
- **A query whose residues are further apart than the index distance cutoff produces no
  hashes**, and therefore no hits. That is now reported as `NO_HASHES` with a warning on
  stderr, not as NOVEL — an unsearchable query is not a discovery. Before this change,
  three of the degenerate queries tried returned NOVEL.

`--nonrigid` composes with the verdict, but on this evidence the effect is theoretical:
across 200 comparisons it changed exactly one verdict, and coverage rose twice and never
fell.

## 5. Deformation metrics (`src/structure/metrics.rs`)

`drmsd` (distance-matrix RMSD) and `max_dist_deviation` compare the internal distances
of a match instead of superposing it, so a motif whose halves swung apart on a hinge
keeps a low dRMSD where its RMSD is large. Available to `--sort-by`, `--format-output`
and `--drmsd`, plus `min_drmsd` per structure. They are read straight off the matched
coordinates, so no superposition is involved.

Ranking honestly did not improve: sorting by dRMSD does not beat sorting by RMSD after
`node_count` on either benchmark. These are reporting and filtering metrics, and the
help text says so.

## Alternatives measured and rejected

| feature | verdict |
| --- | --- |
| **Torsion-angle ENM query sampling (`--enm-sample`)** | **Rejected and removed.** Against `--nonrigid`, post-fix: **−0.047 AP** on the independent serine family with 0/200 replicates positive, −0.004 AP on the 3-residue zinc motif (0/200 positive), against +0.019 AP on the 4-residue zinc motif it was developed on. 3.1x prefilter runtime and 17 crates (`nalgebra`) for one eigendecomposition. Its conformer count could not be chosen by measurement: `ENSEMBLE_CONFORMERS` from 2 through 40 produces byte-identical output, i.e. the ensemble saturates at one deformation direction and its sign. Three real defects in it were found and fixed first (rotations propagating across chain boundaries, a `--nma-rmsd` flag whose Ångström units were fictional, trivial modes crowding out real ones) — those fixes took a cross-chain motif from 0.778 to 0.850 AP but still lost to `--nonrigid` at 0.897. The work is preserved on branch `enm-torsion-sampling`. |
| ANM/NMA query sampling | Rejected. Never beat torsion-ENM on quality and ran at 1.4 s against a 9 ms baseline (~110x), most likely a full 3N×3N eigendecomposition where torsion space is far smaller. |
| DMS / PAS / SOS composite ranking metrics | Rejected. With an identical 1635-structure candidate set none beats plain RMSD as the tiebreaker after `node_count` (RMSD 730/742/843/898 against DMS 729/741/834/896, PAS 730/740/837/880, SOS 668/741/844/895). DMS is a normalised dRMSD, and that does not improve ranking either. |
| `--dist-ratio` (elastic, distance-proportional tolerance) | Rejected and removed. Raised total recall while losing true positives at every early-precision point: widening the tolerance most for long pairs adds exactly the least specific matches. A knob documented as making results worse is a knob nobody can set correctly. |
| An absolute IDF ceiling (`--max-idf`) instead of the relative floor | Not implemented. It measured better at the top of the ranking in an early sweep, but on the metric that sweep used (TP@5FP), which the robustness analysis later showed cannot support a choice between thresholds. Worth revisiting with AP if short-motif precision becomes the priority. |
| Perturbation analytics (`analyze-perturb`) | Not a search feature — byte-identical search output. Judge on the diagnostics alone. |

## Known issues found but not fixed

- **`CompactStructure::build` never resets its `c` variable between residues**, so the
  Cβ approximation for a residue lacking Cβ can fall back on an *earlier residue's*
  carbonyl carbon. Fixing it changes every hash and would invalidate published indices,
  so the existing behaviour is preserved verbatim. Worth addressing on a future
  index-format change.
- **`FoldcompDbReader::new` sorts its lookup by database key**, while
  `get_foldcomp_db_entry_by_name` binary-searches the same lookup by *name*, so
  `read_single_structure(name)` fails for most entries unless `sort_lookup_by_name()` is
  called first. `examples/seqdump.rs` calls it; other callers may not.
- **A 3-column foldcomp-native lookup panics `load_lookup_from_file`**, which indexes
  columns 2 and 3 unconditionally (`data/foldcomp/example_db.lookup` is such a file).
  Pre-existing, and it affects the text and cache paths equally.

## Limits of these measurements

- **Two motif families, one index.** The serine set is a genuinely independent second
  family, which is what makes the `--nonrigid` result credible, but two families on one
  proteome is still a narrow base. Nothing here has been measured at AFDB scale.
- **The serine answer set covers only the S1 clan.** A geometric Ser-His-Asp query
  legitimately also hits S8/S9/S28 peptidases and α/β-hydrolases, which count as false
  positives here, so absolute precision is *understated*. Comparisons between
  configurations are unaffected.
- **The answer sets are index-specific.** `data/serine_answer.tsv` lists entries of
  `index/h_sapiens`; regenerate it with `scripts/build_serine_answer.py` for any other
  index rather than reusing the file.
- **Single machine, warm cache.** All timings are 20-core, 8-thread, warm page cache,
  local NVMe.
