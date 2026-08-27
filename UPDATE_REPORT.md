# folddisco — update report

Everything below lives in **one repository**: `/home/user1/projects/afdbv6/folddisco`.
The two new branches are checked out in separate worktrees so your own checkout was
never touched, but they are ordinary branches in that repo's ref store — `git branch`
from the main checkout lists them.

## 0. What to commit, and where

| branch | worktree | base | contents |
|---|---|---|---|
| `lookup-mmap-cache` | `/home/user1/projects/afdbv6/fd-mmap` | `feature-integration` (`fcd8742`) | mmap-backed caches for both lookups (§3) |
| `widen-chain-ids-compat` | `/home/user1/projects/afdbv6/fd-chainid` | `origin/widen-chain-ids` (`2a756d9`) | multi-character chain IDs, output format unchanged by default (§4) |

Both worktrees are clean; nothing is uncommitted. The only untracked file in your
own checkout is `index/iotest`, which was already there.

**`lookup-mmap-cache` → `feature-integration`.** A plain fast-forward:

```bash
cd /home/user1/projects/afdbv6/folddisco     # already on feature-integration
git merge --ff-only lookup-mmap-cache
```

**`widen-chain-ids-compat` → a PR against `steineggerlab/folddisco:widen-chain-ids`.**
This one is based on master, not on `feature-integration`, because that is the branch
you asked the PR to target. It could not be pushed from here — this environment has no
git credentials (no token, no SSH key, no credential helper), so the push and the PR
are yours to run:

```bash
cd /home/user1/projects/afdbv6/fd-chainid
git push origin widen-chain-ids-compat       # or: git push fork widen-chain-ids-compat
gh pr create --repo steineggerlab/folddisco \
  --base widen-chain-ids --head widen-chain-ids-compat \
  --title "feat: multi-character chain IDs, without changing the default output" \
  --body-file /home/user1/projects/afdbv6/CHAINID_PR_BODY.md
```

Supporting documents outside the repo, none of which need committing:

* `/home/user1/projects/afdbv6/CHAINID_PR_BODY.md` — the PR body for the chain-ID branch.
* `/home/user1/projects/afdbv6/bench_data/LOOKUP_LOAD_BENCH.md` — every measurement
  behind §3, including the machine spec and the commands used.
* `/home/user1/projects/afdbv6/bench_data/lookup/` — 45 GB of downloaded lookups and
  built caches, reclaimable when you are done re-checking numbers.

## 1. Non-rigid motif search

Live at `fcd8742` as `--expand-radius` (default **1**) and `--nonrigid`
(= `--expand-radius 2`). `src/controller/expand.rs`.

### What is hashed, and why a rigid lookup misses a deformed motif

folddisco indexes ordered **residue pairs**, not residues. The default hash type
`PDBTrRosetta` takes the feature vector
`[aa1, aa2, ca_dist, cb_dist, ca_cb_angle, theta1, theta2]` and packs it into 30
bits (`geometry/pdb_tr.rs:73-75`):

* two residue identities, 5 bits each — **exact**, never expanded;
* **two distances** (Cα–Cα, Cβ–Cβ), 16 bins over [2, 20] Å → **1.2 Å per bin**
  (`utils/convert.rs:3-4`, width `(max-min)/(nbin-1)` at `geometry/core.rs:209-216`);
* **three angles** (the Cα–Cβ angle and two torsions), each stored as `sin` **and**
  `cos` separately, 4 bins over [-1, 1] each (`pdb_tr.rs:17,46-70`).

Five continuous features become ten binned fields. Lookup is exact equality on the
packed `u32` (`retrieve.rs:137-141`), so a motif that is present but slightly
deformed puts one feature on the far side of a bin boundary and the hash simply
differs. The test at `expand.rs:400-423` drifts Cα and Cβ distance by 0.1 Å across
the 7.4 Å boundary and the hash changes.

### The mechanism: a Hamming ball in bin-offset space

"Joint" means distances and angles go into **one** list of tolerant dimensions —
indices `[2,3]` for distances and `[4,5,6]` for angles (`controller/feature.rs:269-286`),
built by `FeatureExpander::new` (`expand.rs:93-130`). A radius-2 combination can
therefore be (dist, dist), (dist, angle) or (angle, angle).

Each dimension gets signed sub-step offsets covering `[-tol, +tol]`, nearest first
(`substep_offsets`, `expand.rs:208-224`):

```rust
steps = ceil(tol / bin_width).clamp(1, MAX_SUBSTEPS = 8)
magnitude = tolerance * step / steps
```

Only the **widest** threshold in `-d`/`-a` is used (`widest_tolerance`,
`expand.rs:199-201`). `for_each_neighbor` then walks levels `1..=radius`, and
`visit_level` picks `level` *distinct* dimensions and takes the full cartesian
product **of that subset only** (`expand.rs:150-193`). So it is a Hamming ball of
radius r, not the product over all dimensions.

Every variant is domain-repaired before hashing — `sanitize_perturbed_feature`
(`geometry/core.rs:308-319`) wraps torsions at ±π and reflects `acos`-derived
angles. Distances are deliberately **not** clamped (`core.rs:300-307`) so the query
encodes exactly as the index does.

### Cost

With D tolerant dimensions and k offsets each, hashes per pair is
`Σ_{l=1..r} C(D,l)·k^l`, times `1 + aa_variants.len()` insert attempts
(`query.rs:274-292`). At defaults (D = 5; `-d 0.5` → `ceil(0.5/1.2)` = 1 → k = 2;
`-a 5` → 0.087 rad / 0.667 → k = 2):

| radius | hashes per pair | note |
|---|---|---|
| 0 | 1 | expansion off (`ToleranceConfig::none`) |
| **1** | **10** | the default; asserted at `expand.rs:257-263` |
| 2 | 50 | `--nonrigid` |
| 3 | 130 | |

Bounds: `MAX_SUBSTEPS = 8` per dimension (`expand.rs:36`) and
`MAX_HASHES_PER_PAIR = 4096` as an early exit (`query.rs:24,293`). Dedup is
implicit — `insert_binned_hash` skips a key that is already present
(`query.rs:99-103`). Note the 4096 cap counts *attempts*, not unique hashes.

### What keeps the extra candidates honest

Expanded hashes enter the inverted-index prefilter identically
(`count_query.rs:120-160`). Survivors are then re-derived from the target structure
and re-hashed, so a candidate must reproduce one of the query hashes **exactly**
(`retrieve.rs:120-141`), gated by the Cα-distance agreement test
`(curr_dist - dist).abs() < ca_distance_cutoff` (`retrieve.rs:104-107`, default
`--ca-distance 1.0`). Connected components then give node counts, and Kabsch/LMS
give RMSD and dRMSD (`retrieve.rs:766-800`).

**Every quality filter is off by default.** CLI defaults are 0 for
`--covered-node`, `--max-node`, `--rmsd`, `--connected-node` and `--drmsd`
(`main.rs:62-79`), and the filters treat 0 as "skip" (`filter.rs:81-124,207-248`).
The `StructureFilter::default()` / `MatchFilter::default()` constructors that carry
sensible values (`covered_node_ratio` 0.8, `rmsd` 1.0) are **never called by the
CLI** — dead code.

That matters, because `feature_evaluation.md:122-141` is explicit that expansion is
only a win once `--max-node` is on: on the matched 4-residue zinc query, F1 goes
0.9421 → **0.9641** with `--max-node`, but 0.9265 → **0.9226** without it. Cost:
`+11%` on that command (157 → 175 ms) and `~73%` slower on M-CSA
(`feature_evaluation.md:182-190`).

### Reverted, and confirmed gone

`--enm-sample` built a torsion-space elastic network model of the query, sampled
low-frequency modes into conformers and searched the union of their hashes.
`3cb2946` removed it: at `fcd8742` a grep for `enm|nma|nalgebra` finds **only prose
in `feature_evaluation.md:453-460`** — no code path, no flag, no dependency.
`d29c1f4` removed the rare-hash IDF filter; `is_primary` survives as a vestige with
no production consumer, documented honestly at `query.rs:62-73`.

**Do not cite `fd-bench/nonrigid_search.md` for current behaviour.** It predates
both reverts, documents `--enm-sample` as shipped, and its `--nonrigid` row is
numerically the margin-3 *filtered* row — so even its non-ENM figures were measured
with the removed filter, against a different answer set and metric.
`feature_evaluation.md` was written after both reverts and says so at line 21.

### Concerns worth a look

1. **`MAX_SUBSTEPS = 8` silently reinstates the gap it exists to prevent.** Beyond
   8 bins of tolerance (`-d > 9.6 Å` at defaults) `steps` clamps, spacing exceeds
   one bin, and bins are skipped again — the exact failure `expand.rs:203-207` was
   written against. The test only goes to 3.0 Å (`expand.rs:287-307`).
2. **Under the 4096 cap, which dimension pairs survive is order-dependent, not
   geometry-dependent.** Dimensions are pushed distance-first (`expand.rs:100-121`)
   and levels iterate lexicographically, so distance–distance joints systematically
   outlive angle–angle ones.
3. **First-writer-wins dedup can shadow another pair's observed hash.**
   `insert_binned_hash` returns early if the key exists (`query.rs:99-103`) and
   `map_query_and_retrieved_residues` reads the single stored pair
   (`retrieve.rs:645`), so an expanded variant of pair A inserted first steals the
   vote pair B's *observed* hash would have cast. Expansion enlarges that collision
   surface; RMSD is the only backstop.
4. **Out-of-window Cβ distances still alias into the residue-identity field.**
   `expand.rs:481-482` deliberately exempts distances from the bit-field bound test.
   Pre-existing (`core.rs:300-307`), but exercised more often once offsets move
   distances.
5. **Doc mismatch:** the README and the `query` help say results are sorted by node
   count then RMSD; the code default is **IDF then RMSD** (`sort.rs:227-231`).

## 2. IDF scoring and filtering

This section answers the question you opened with: why the EF-hand query
(`1RFJ`, `-q A21,A23,A25,A27,A32`) put an 852,272 hit above the 110.688 self-match.
**It is not a scaling artefact of the IDF formula. The `idf` column in per-match
output is summed over the wrong set of edges.**

### Two different quantities, both printed in a column called `idf`

**Per-structure** (`--per-structure`), `count_query.rs`:

```rust
131:  (lookup.len() as f32 / hash_count as f32).log2()   // idf of one query hash
154:  entry.idf_sum += idf;                              // once per matching hash
203:  merged_entry.idf_sum *= (lookup_entry.2 as f32).powf(-lp);   // nres^-lp
```

i.e. `nres^(-lp) · Σ log2(N / df(h))` over matched query hashes, where
`N = lookup.len()` and `df(h)` is the posting-list length. Posting lists are
deduplicated per structure at build time (`controller/mod.rs:344-345`), so each
distinct hash contributes at most once per structure. It is a **plain sum** —
never a mean, never divided by node or hash count.

**Per-match** (the default mode — `mode.rs` maps no flags to `PerMatch`), which is
what your command produced: `retrieve.rs:231` / `:421` call
`calculate_subgraph_idf` and that value is passed as `avg_idf` into
`MatchResult::new` (`result.rs:135,148`) despite no averaging happening.

### Why 852,272

Follow the per-match path:

```rust
// retrieve.rs:212-214
let graph = create_index_graph(&indices_found);
let connected = connected_components_with_given_node_count(&graph, node_count);
```

`connected_components_with_given_node_count` retains components with **`len() >=
node_count`** (`graph.rs`), i.e. *at least* the motif size, with no upper bound.
Then, per component:

```rust
// retrieve.rs:217-231
let subgraph = graph.filter_map(
    |node, _| if component.contains(&graph[node]) { Some(graph[node]) } else { None },
    |_, edge| Some(*edge)          // <- every edge of the component is kept
);
let subgraph_idf = calculate_subgraph_idf(&subgraph, query_map);
```

and

```rust
// retrieve.rs:715-728
for edge in subgraph.edge_indices() {
    if let Some(&(_, _, idf)) = query_map.get(&hash) { total_idf += idf; }
}
```

So the reported `idf` is the IDF sum over **every edge of the whole connected
component the match was drawn from**, not over the ~10 edges of the five matched
residues. Meanwhile the `node_count` printed in the same row is the number of
`Some` entries in `matching_residues` — five. **The two columns describe different
objects.**

A single hash's term is bounded by `log2(N)` ≈ 17.8 at PDB scale, so nothing about
one rare hash can reach 852,272. What scales is the **number of edges in the
component**, which grows roughly with the square of its node count. A small,
highly repetitive helical protein whose residue pairs densely hit the expanded
query hash set produces a large component and a large score; `1RFJ` matching its
own EF-hand sits in a small component and scores 110.688. Both rows are the same
column, and neither number is comparable to the other.

Two further defects compound it, found independently in the same read:

* **Expanded variants inherit the observed parent's IDF.** `calculate_idf_for_hash`
  is called once for the observed hash (`query.rs:243`) and that value is stored on
  every expanded variant (`query.rs:268,280,288`). A hit on a *common* neighbouring
  bin is therefore credited with the *rare* parent bin's rarity. `count_query.rs`
  recomputes IDF per hash from the index instead — so the two code paths use two
  different conventions in one query.
* **The per-structure sum is not deduplicated per edge.** `idf_sum += idf` fires
  once per matching hash (`count_query.rs:154`) while `edge_count` fires once per
  edge (`:165-168`). At radius 2, a target matching many neighbour bins of *one*
  edge can outrank a target matching one bin of *every* edge — and IDF is the
  default first sort key.

I could not reproduce the exact 852,272 here: it needs your PDB-wide index, which
is not on this machine (`/data/xray2/...` does not exist here). The mechanism above
is read directly off the code; the specific value is not verified. Point me at an
index and I will confirm it.

### Length penalty

`lp = length_penalty_power.unwrap_or(0.5)` (`count_query.rs:91`), flag
`--length-penalty`, default **0.5**. Multiplies the sum by `nres^-0.5` of the
*target*, damping long structures and boosting short ones. `lp = 0` disables it.
Note it damps only with the square root of length while the sum grows with the
number of matched hashes, so it does not offset the effect above.

### Every live knob that touches the score or the candidate set

Hash admission, before scoring (`count_query.rs`):

| flag | default | effect |
|---|---|---|
| `--freq-filter <f>` | **off** | skips hashes with `df/N > f` (`:124-128`). The help says `[0.0]`, but passing `0` skips *every* hash and empties the result. |
| `--sampling-ratio` / `--sampling-count` | off | keep only the rarest hashes (`:226-257`) |
| `--length-penalty` | 0.5 | as above |

Per-structure, before matching (`filter.rs:81-105`) — all off unless noted:

| flag | default | compares |
|---|---|---|
| `--total-match` | 0 = off | `total_match_count >= X` |
| `--covered-node` | 0 = off | `node_count >= X` |
| `--covered-node-ratio` | 0.0 = off | `node_count / residue_count >= r` |
| `--score` | 0.0 = off | `result.idf >= X` — the length-penalised **sum** |
| `--num-residue` | 50000 — **on** | `nres <= X` |
| `--plddt` | 0.0 = off | `plddt >= X` |

Per-structure, after matching (`filter.rs:108-124`): `--max-node`,
`--max-node-ratio`, `--rmsd`, `--drmsd` — all 0 = off.

Per-match, applied only in `--per-match` / `--web` (`filter.rs:207-251`):
`--connected-node`, `--connected-node-ratio`, `--rmsd`, `--tm-score`, `--gdt-ts`,
`--gdt-ha`, `--chamfer`, `--hausdorff`, `--drmsd` — all off by default — plus
**`--score` again, here against the subgraph IDF** (`:217`). One flag, two
incomparable scales. The e-value cutoff is hardwired to `f64::MAX`
(`query_pdb.rs:616`, "Currently not used").

`--top` truncates *after* sorting by IDF desc and *before* residue matching
(`query_pdb.rs:559-566`). Default sort is IDF desc → RMSD asc (`sort.rs:227,467`).

### E-value

`result.rs:364-384`. `mu = 4.2161·exp(0.0489 l) + 3.6661`,
`lam = 0.2894·exp(-0.0762 l) + 0.0316`, `k = exp(lam·mu)/10546`,
`E = k·m·l·exp(-lam·x)`, then squashed by `E·m/(E+m)`. `x` is the per-match
subgraph IDF — so it inherits everything above. Its own comments read
`// TODO: Finalize fitting.` and `// WARNING: This is not validated yet.` It is
**not printed by default** (`//"e_value",` is commented out of
`MATCH_RESULT_DEFAULT_COLUMNS`, `result.rs:342`); reachable via
`--format-output e_value`.

### The reverted rare-hash filter

`7bac6ea` raised `EXPANDED_HASH_IDF_MAX_EXCESS` from 3.0 to 5.0: it dropped a
tolerance-expanded hash whose own IDF exceeded its observed parent's by more than
5 bits. `d29c1f4` removed the constant and its 49-line block; it appears nowhere
in `src/` now. The stated reason: under the author's protocol (Sens@1FP, F1) the
shipped margin was "worse-or-equal in 16 of 16 cells on both metrics and better in
none", with no runtime justification (158.7 ms vs 159.8 ms) — and it had been "the
*entire* difference between this branch's default path and master". Preserved on
branch `rare-hash-idf-filter`.

`foldcomp-shard-eval/FOLDDISCO-IDF-FILTER.md` is entirely about that removed
filter, is self-marked "Superseded in part", and its `sed` recipe will no-op on
`fcd8742`. Nothing in it describes the live scoring path.

### What the repo already says to do about it

`feature_evaluation.md:266-272`, on the live code and the PDB index (230,655
structures):

> a matched query sorted by raw IDF puts sparse partial matches on rare NMR
> entries at the top (`2n8a`, node_count 2, IDF 26,829) rather than the correct
> family… On a database that size, IDF needs a node-count floor to be usable.

So the only in-repo remedy today is **gate on node count** (`--max-node`,
`--covered-node`) and **sort by node count** (`--sort-by node_count,rmsd`, which
`sort.rs:473-486` supports), rather than trusting the score. Line 457 notes
`--max-idf` was never implemented. No per-node normalisation exists anywhere in
the tree — the only trace is a disabled `idf_max_per_edge` field
(`count_query.rs:65,76`).

### What I would fix, in order

1. **Score the matched subgraph, not the component.** `calculate_subgraph_idf`
   should be given the edges induced by `retrieved_indices` — the residues actually
   reported — not the whole component. That alone makes per-match `idf` comparable
   across hits and bounded by `C(node_count, 2) · log2(N)`.
2. **Pick one IDF convention.** Either recompute per hash as `count_query` does, or
   inherit the parent's as `query.rs` does — not both in one query.
3. **Rename one of the two columns.** Printing a component-wide sum and a
   per-hash-sum both as `idf` is what made this hard to see from the output.
## 3. Cache data structure

Two lookup loads used to dominate the start of every query. Both are now mapped
rather than parsed, from caches whose **on-disk layout is the in-memory layout**.

### The problem, as a number

`load_lookup_from_file` returned `Vec<(String, usize, usize, f32, usize)>`. On the
3.07 GB `afdb50_v4` lookup (53,665,860 entries) that is 53.7 M simultaneously
live `String`s: a counting allocator measured **53,779,018 allocations and 8.08 GB
of live heap for 1.32 GB of actual name bytes**. Of the 1.79 s decode, 1.33 s
(74%) was the allocator, not the data — the field decoding itself is at most
0.21 s. Because every one of those allocations stays live, the kernel had to zero
and map 6.02 GB of fresh anonymous pages: 1.47 M extra minor faults and ~8.8 s of
*system* CPU. That is also why it did not scale — best at 8 threads for 2.0x, then
*slower* at 20 (1.58 s → 1.88 s) as system CPU climbed 5.95 s → 10.87 s on glibc
arena and kernel `mmap_lock` contention.

The Foldcomp DB path was worse: `FoldcompDbReader::new` parsed **and sorted** both
the `.lookup` and the `.index` on every query, with no cache at all. On
`afdb_uniprot_v6` (10.78 GB, 241,070,489 entries) that is 6.14 s and 27.9 GB
resident before a single structure is read.

### The layout

`src/utils/pod_cache.rs`. One header, then a plain array of a `#[repr(C)]` POD
record, then a single blob of names.

```text
offset  0   magic         [u8; 8]   FDLOOKUP | FDFCLKUP | FDFCIDXT
offset  8   version       u32
offset 12   pad           u32       keeps the counts 8-aligned
offset 16   count         u64       number of records
offset 24   source_len    u64       byte length of the file this was built from
offset 32   source_mtime  u64       mtime of that file, in nanoseconds
offset 40   names_len     u64       byte length of the name blob
offset 48   records       count * size_of::<R>()      <- 8-aligned, cast as &[R]
            names         names_len bytes of UTF-8    <- sliced, never copied
```

`HEADER_SIZE` is 48 — a multiple of 8 — so the record table starts 8-aligned and
is handed to callers as a `&[R]` cast straight out of the mapping. **This one
constant is the whole trick.** The previous v2 cache already had a fixed-stride
40-byte POD record; its 44-byte header put the table at a 4-aligned offset, which
is why it had to be decoded rather than cast.

| record | size | fields |
|---|---|---|
| `LookupRecord` | 40 B | `id u64`, `nres u64`, `db_key u64`, `name_offset u64`, `plddt f32`, `name_len u32` |
| `FoldcompLookupRecord` | 24 B | `key u64`, `name_offset u64`, `name_len u32`, `_pad u32` |
| `FoldcompIndexRecord` | 24 B | `key u64`, `offset u64`, `length u64` |

Field order is chosen so each record is padding-free and 8-aligned: the 8-byte
fields first, the 4-byte ones paired. Padding-free matters because it means the
bytes on disk are exactly the value, in both directions.

### What makes the cast sound

`CacheRecord` is an `unsafe trait` whose contract is the safety argument, stated
once where it can be checked:

> `#[repr(C)]`, no padding bytes, alignment at most 8, size a multiple of 8 (so
> consecutive records stay aligned), and valid for every bit pattern — i.e. built
> only from integer and floating-point fields, with no enums, references, `bool`s
> or `NonZero` types.

The alignment is additionally asserted at load rather than assumed, because a
misaligned cast is undefined behaviour rather than a wrong answer.

When the index directory is read-only — a normal deployment — the cache is built
in memory instead. That image is held as a `Vec<u64>`, not a `Vec<u8>`: a `Vec<u8>`
is only guaranteed 1-aligned, and most allocators happening to return 16 is not
something an alignment invariant should rest on.

### What is validated, and why that set

Every load checks, in O(1): magic, version, the **length and mtime** of the source
file, and that the file size is exactly what the header's counts imply.

* Length *and* mtime, not mtime alone: an mtime comparison by itself accepts a
  stale cache whenever the source is rewritten inside one filesystem timestamp
  tick (1 s on ext3, HFS+, most NFS servers), which serves wrong names and pLDDTs
  with no signal at all.
* The size check pins the counts, so a corrupted `count` cannot decode as a
  short-but-well-formed cache — a zeroed count would otherwise look like an empty
  index and make a search silently find nothing.

Then one **parallel pass over the record table**, which is the part that is not a
bare cast and is deliberate:

| table | checked |
|---|---|
| index lookup | non-empty name, in-bounds name range, `id < len()` |
| Foldcomp lookup | non-empty name, in-bounds name range, strictly ascending keys |
| Foldcomp index | non-zero entry length, strictly ascending keys |

Every field of an all-zero record is *individually* in range, so a torn write
would otherwise decode as a real entry. For the index lookup its `id` of 0
attributes those hits to the first structure in the index — and `id` indexes a
vector of exactly `len()` entries, so one past the end panics mid-search. For the
Foldcomp tables a zeroed record serves another entry's structure for key 0.
Ascending keys are also exactly what the per-hit binary search relies on.

Anything that fails sends the load back to parsing the text, which is always
still there. That pass is most of the load time; a bare cast with no pass is
~20 µs, if that trade is ever wanted.

Two things are deliberately absent, with the reasoning inherited from the v2
cache: **no checksum**, because computing one costs a full pass, which is the
whole thing the cache exists to avoid, and after the size and identity checks no
reachable write path leaves a table that is corrupt yet passes; and **no
temp-file-and-rename**, because cache content is a pure function of the source, so
two processes writing the same cache emit identical bytes from offset 0 and a
reader that catches a write in progress sees a short file and rejects it. Both
arguments hold only while the bytes stay deterministic — putting a timestamp or a
thread id in the body would break them and bring atomicity back into scope.

### Names

Names live in one contiguous blob and come back as `&str` slices of the mapping.
Nothing is copied, so the entries a query never reports are never touched — a
top-1000 query resolves 1000 names, not 53.7 M. Bounds and UTF-8 are checked at
resolution, not up front: validating a multi-gigabyte blob on every load would
defeat the purpose, and a name that is not UTF-8 reads back empty rather than
panicking.

### Ordering

The Foldcomp tables are sorted by key **when the cache is built**, so no load ever
sorts. That removed `sort_lookup_by_id` and `sort_lookup_by_name` entirely, along
with the re-sort `read_compact_structure` did on every `db:name` read.

Name resolution is therefore a scan, not a search — and that is the deliberate
trade. Keys are resolved **once per hit** (binary search, ~5 µs); names at most
**twice per query**, when the query structure is given as `db:name`. A
name-ordered index would have added ~0.96 GB to the cache and a sort of 241 M
names to its build to save microseconds on a per-query operation. `keys_of_names`
resolves any number of names in one pass, for callers that ask for several.

### Measured (warm, `/usr/bin/time`)

**Index lookup** — `afdb50_v4`, 53,665,860 entries, 3.07 GB:

| | before | after | |
|---|---|---|---|
| load | 2.98 s | **17.5–17.9 ms** | **170x** |
| peak RSS | 8.55 GB private anon | **2.10 GB page-cache** | 4.1x less, and evictable |
| allocations | 53,779,018 | ~250 | |
| cold load (parse + build + map) | 4.64 s | 3.34–3.45 s | once |
| cache size | 3,462,598,388 B (v2) | 3,462,598,392 B (v3) | +4 B of header padding |

**Foldcomp DB lookup** — `afdb_uniprot_v6`, 241,070,489 entries, 10.78 GB:

| | before | after | |
|---|---|---|---|
| load | 6.14 s | **45.8–50.3 ms** | **~130x** |
| peak RSS | 27.9 GB private anon | **5.66 GB page-cache** | 4.9x less |
| key → name (per hit) | — | ~4.4–5.1 µs | binary search |
| name → key (≤2x per query) | — | ~98–102 ms | parallel scan |
| cold load (parse + sort + build + map) | 6.14 s *every query* | 10.7 s | once |
| cache size | — | 11,725,383,048 B | |

For scale, on that same 241 M-entry file: Foldcomp's own `read_lookup` with the
`stable_sort` its `USE_LOOKUP` triggers takes **141.6 s** (~3000x), and MMseqs2's
`readLookup` plus the sort `open()` always does takes **29.5 s** (~620x).

Fusing the two Foldcomp validation checks into one pass rather than two matters at
this size: 90 ms → 47 ms, because the table is 5.8 GB and reading it twice doubled
the only per-entry work left.

### Two things to decide

1. **The Foldcomp caches are large** — 11.7 GB for a 10.78 GB lookup, being
   24-byte records plus a copy of the names. Narrowing `key` to the `u32` that
   MMseqs2-format DBs actually use takes it to about 8 GB. The format is
   versioned, so this needs no migration.
2. **The validation pass is most of the remaining time.** Dropping it reaches
   ~20 µs and 1.9 MB, at the cost of the torn-write protection above. It is a
   two-line change if that is the trade you want.

### Related finding, not fixed here

MMseqs2 and Foldseek feel instant on lookups because **their search path never
loads the lookup**: it is gated on `USE_LOOKUP` (`DBReader.cpp:138`, flag `= 8` in
`DBReader.h:310-315`), and that flag appears nowhere under `workflow/`,
`prefiltering/` or `alignment/`. When they do read one, `readLookup`
(`DBReader.cpp:1138-1158`) is a serial loop with one `std::string` per entry at
53.4 ns/entry — slower per entry than what folddisco had. Foldcomp *does* have a
`.cache`, but it caches the **index only**: on a cache hit `make_reader` sets
`reader->lookup = NULL` and returns early (`database_reader.cpp:74-86`), so
`--use-cache` silently gives up name→key resolution and every `--id-mode 1` id
warns "not found". Its `load_cache` also validates nothing, and its `NULL` return
is dereferenced one line later — `chmod 000` on the cache segfaults.

## 4. Multi-character chain IDs

Branch `widen-chain-ids-compat`, two commits, targeting `widen-chain-ids`.

### What was broken

mmCIF `auth_asym_id` can be more than one character — large cryo-EM entries use
`AA`, `10` and so on — but chain IDs were a single `u8`. The CIF reader worked
around that by *rejecting* any multi-character `auth_asym_id` and silently falling
back to `label_asym_id`, so folddisco reported chain labels that do not exist in
the file's own numbering. On a copy of `1G2F.cif` with the auth IDs renamed to
`AA`/`BB`/`10`/`DD`/`EE`/`FF`:

```text
before:  1G2F_multichain.cif  4  0.0000  F207,F212,F225,F229   <- no chain F in that file
         1G2F_multichain.cif  4  0.1633  E107,E112,E125,E129   <- no chain E either
after:   1G2F_multichain.cif  4  0.0000  FF_207,FF_212,FF_225,FF_229
         1G2F_multichain.cif  4  0.1633  10_107,10_112,10_125,10_129
```

### Output format — the part that had to stay compatible

`matching_residues` and `query_residues` keep the `A21,A23` spelling. A field
switches to `AA_250` / `10_250` **only** when a chain ID in it would otherwise be
unreadable, i.e. when it is multi-character or numeric: `A21` re-parses because a
leading letter can only be the chain, while `AA250` and `10250` cannot be split at
all. The decision is taken **per field**, so no field is ever half in one spelling
and half in the other. `--chain-sep` forces the separated form everywhere, for
callers who would rather parse one fixed format.

`-q` accepts both spellings (`B57` == `B_57`), so any `matching_residues` field
folddisco prints can be pasted straight back in as a query. A unit test round-trips
every shape the formatter emits back through `parse_query_string`.

Evidence that the default output is unchanged: nine query variants over
`data/serine_peptidases` — per-match, per-structure, `--header`, `--superpose`,
`--format-output`, `--serial-index`, chain-less queries, and the Foldcomp FCZ DB
path — are byte-identical to the `widen-chain-ids` tip, as are the index table,
offset file and lookup that `folddisco index` produces. **Existing indices stay
valid**: chain IDs are never serialised into them, only derived when a structure is
read.

### Representation

`ChainId` holds up to 8 printable ASCII bytes inline and stays `Copy`. This is the
one deliberate departure from PR #44, which used `String`: one chain ID is carried
per residue of every candidate structure, so `Copy` keeps that path's byte
comparison and allocation-free behaviour and adds no `.clone()` calls.

`Atom::chain` stays a single byte on purpose — `Atom::from_c` transmutes Foldcomp's
C `atom_t` into it (48 bytes, `chain: c_char` at offset 20), so that layout is
pinned. Widening happens one level up in `AtomVector` / `Structure` /
`CompactStructure`, which is where the CIF parser can supply the full identifier.
The FCZ path is covered by the byte-identical comparison above.

### A pre-existing bug this exposed (second commit, separable)

`CompactStructure::build` writes a residue out when the *following* one starts, but
took the chain ID from `atom_vector.chain[idx]` — and at that moment `idx` is the
first atom of the next residue. **Every residue that ended a chain was filed under
the chain that came after it.**

On `data/serine_peptidases/4cha.pdb`, whose chain A holds residues 1–11 and whose
chain B starts at 16, `chain_per_residue` reported chain A as residues 1–10 and put
residue 11 into chain B: `-q A11` found nothing, and `-q B11` matched a residue of
chain A. The same at all six chain boundaries of that file, and at every chain
boundary of every multi-chain structure.

`folddisco analyze` confirms it from a different code path — five enriched positions
move to the chain they belong to, and two that had been colliding with a real
residue of the next chain reappear:

```text
4cha.pdb B11  -> A11     (residue 11 ends chain A)
4cha.pdb C146 -> B146    (residue 146 ends chain B)
4cha.pdb E245 -> C245    (residue 245 ends chain C)
4cha.pdb F10  -> E10     (residue 10 ends chain E)
4cha.pdb G146 -> F146    (residue 146 ends chain F)
1azw.pdb         A313    (new: was reported as B313, colliding with the real one)
1l7a.pdb         A318    (new: was reported as B318, likewise)
```

**This changes results for multi-chain structures** — correctly. Indices do not need
rebuilding. It is a separate commit so it can be reviewed or dropped on its own.

`b_factors` in the same loop has the same one-residue shift. Left alone: its only
consumer is `get_avg_bfactor`, where an off-by-one in the summands barely moves the
mean, and correcting it would change reported pLDDT for reasons unrelated to chain
IDs.

## 5. Test status

`cargo test --release`:

| branch | result |
|---|---|
| `lookup-mmap-cache` | 142 passed, 0 failed, 7 pre-existing `#[ignore]` |
| `widen-chain-ids-compat` | 117 passed, 0 failed, 7 pre-existing `#[ignore]` |

New tests worth knowing about, because each encodes a failure mode rather than a
happy path:

* `pod_cache`: a source rewritten inside one mtime tick is rejected; truncated,
  wrong-magic, wrong-version, zeroed-count and over-large-count caches are all
  rejected; an unwritable cache path degrades to an aligned in-memory image.
* `index::lookup`: a zeroed record and an out-of-range `id` are rejected and the
  text is reparsed; a 4-column legacy lookup reuses the id as `db_key`.
* `structure::io::fcz`: the mapped tables match the text files entry for entry in
  both passes (built, then mapped); the cache path does **not** collide with
  Foldcomp's own `<index>.cache`; a zeroed record is rejected.
* `structure::core`: the last residue of a chain keeps its own chain, with ground
  truth from `4cha.pdb`. Verified to fail on the first chain-ID commit and pass on
  the second.
* `controller::result`: the legacy spelling is preserved, `_` placeholders survive,
  an ambiguous chain switches the whole field, and every emitted shape round-trips
  through the query parser.
