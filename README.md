# Folddisco

<p align="center">
<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/steineggerlab/folddisco/master/.github/img/folddisco_logo_with_light.png">
<img src="https://raw.githubusercontent.com/steineggerlab/folddisco/master/.github/img/folddisco_logo_without_light.png" max-height="300px" height="300" display="block" margin-left="auto" margin-right="auto" display="block"/>
</picture>

Folddisco is tool for searching discontinuous motifs in protein structures.
It is designed to handle large-scale protein databases with efficiency, enabling the detection of structural motifs across thousands of proteomes or millions of structures.

## Publications
[Kim H, Kim RS, Mirdita M, Yoon J, Steinegger M. Structural motif search across the protein-universe with Folddisco. Nature Biotechnology, (2026)](https://www.nature.com/articles/s41587-026-03162-9)

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/folddisco.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/folddisco) [![Github All Releases](https://img.shields.io/github/downloads/steineggerlab/folddisco/total.svg)](https://github.com/steineggerlab/folddisco/releases/latest) 

## Webserver 
Search protein structures motifs against the [AlphaFoldDB](https://alphafold.ebi.ac.uk/) and [PDB](https://www.rcsb.org/) in seconds using the Folddisco webserver ([code](https://github.com/soedinglab/mmseqs2-app)): [search.foldseek.com/folddisco](https://search.foldseek.com/folddisco) 🚀

## Installation
```bash
# Install from Bioconda
conda create -n folddisco -c conda-forge -c bioconda folddisco

# Install through docker
docker pull ghcr.io/steineggerlab/folddisco:master

# Precompiled binary for Linux x86-64
wget https://mmseqs.com/folddisco/folddisco-linux-x86_64.tar.gz; tar xvfz folddisco-linux-x86_64.tar.gz; export PATH=$(pwd)/folddisco/bin/:$PATH

# Precompiled binary for Linux ARM64
wget https://mmseqs.com/folddisco/folddisco-linux-arm64.tar.gz; tar xvfz folddisco-linux-arm64.tar.gz; export PATH=$(pwd)/folddisco/bin/:$PATH

# macOS (universal, works on Apple Silicon and Intel Macs)
wget https://mmseqs.com/folddisco/folddisco-macos-universal.tar.gz; tar xvfz folddisco-macos-universal.tar.gz; export PATH=$(pwd)/folddisco/bin/:$PATH
```

**Compile from source**

Compiling from source requires the Rust toolchain (Cargo). Installation instructions are available [here](https://www.rust-lang.org/tools/install).

```bash
git clone https://github.com/steineggerlab/folddisco.git
cd folddisco
cargo install --features foldcomp --path .
```

## Quick start
Folddisco queries a database of precomputed geometric hashes computed from structures. 

### Download pre-build database 
You can download the pre-built human proteome index and use it to search for a common motif, like a zinc finger.

This example is fully self-contained. You can copy and paste the entire block into your terminal.

```bash
# Download human proteome index. Use wget or aria2 to download the index.
cd index
aria2c https://opendata.mmseqs.org/folddisco/h_sapiens_folddisco.tar.lz4

# Extract the index
lz4 -dc h_sapiens_folddisco.tar.lz4 | tar -xvf -
cd ..
```

#### Pre-built Indices

Download pre-built index files:
- [Human proteome](https://opendata.mmseqs.org/folddisco/h_sapiens_folddisco.tar.lz4)
- [E. coli proteome](https://opendata.mmseqs.org/folddisco/e_coli_folddisco.tar.lz4)
- [AFDB proteome of 16 model organisms](https://opendata.mmseqs.org/folddisco/afdb_proteome_v4_folddisco.tar.lz4)
- [Swiss-Prot](https://opendata.mmseqs.org/folddisco/afdb_swissprot_v4_folddisco.tar.lz4)
- [AFDB50](https://opendata.mmseqs.org/folddisco/afdb50_v4_folddisco.tar.lz4) 
- [ESM30](https://opendata.mmseqs.org/folddisco/highquality_clust30_folddisco.tar.lz4)
- [PDB](https://opendata.mmseqs.org/folddisco/pdb_folddisco.tar.lz4)
- To get the old version of Folddisco indices, please **visit** https://opendata.mmseqs.org/folddisco/
  - `*.tar.gz` indices are legacy indices (version 1.0), which are not compatible with version 2.0. 
    Please use `*.tar.lz4` indices for version 2.0.
  - **AFDB50** (`afdb50_v4_folddisco*` + `afdb50_v4*`)
  - **ESM30** (`highquality_clust30_folddisco*` + `highquality_clust30*`)

### Build an custom index 
The command below will read all PDB or mmCIF from `serine_peptidases` folder and generate an index `serine_peptidases_folddisco`.
```bash
folddisco index -p data/serine_peptidases -i index/serine_peptidases_folddisco
```

### Querying a Single Motif
To search for a specific structural motif, you'll use three main flags:
-   **`-p`**: Provides the query protein's structure file (PDB/mmCIF).
-   **`-q`**: Specifies the comma-separated list of residues that form your motif.
-   **`-i`**: Points to the target database index you want to search against.

If you omit the **`-q`** flag, `folddisco` defaults to a "whole structure" search. It will find all possible motifs from your entire query protein and search for them in the index.

```bash
# Search for the catalytic triad from 4CHA.pdb against the indexed peptidases.
folddisco query -i index/serine_peptidases_folddisco -p query/4CHA.pdb -q B57,B102,C195
```
#### Residue & motif syntax
We allow to customize the query motif using some motif syntax.
* **Residues:** `B57` = chain `B`, residue number `57`. Ranges are inclusive and may
  repeat the chain on the end: `1-10`, `F204-215` and `F204-F215` all work. A range
  cannot span two chains.
* **Lists:** comma-separated: `B57,B102,C195`.
* **Substitutions:** `:<ALT>` allows alternatives:
  * Single amino acid: `164:H`
  * Set: `247:ND` (Asp or Asn)
  * Wildcard/categories:
    * `X`: any amino acid
    * `p`: positively charged (Arg, His, Lys)
    * `n`: negatively charged (Asp, Glu)
    * `h`: polar (Asn, Gln, Ser, Thr, Tyr)
    * `b`: hydrophobic (Ala, Cys, Gly, Ile, Leu, Met, Phe, Pro, Val)
    * `a`: aromatic (His, Phe, Trp, Ty)

### Searching Multiple Motifs (Batch Mode)
To search for many motifs at once, you can provide a single query file to the **`-q`** flag (and omit the `-p` flag).

This file must be a **tab-separated** text file with these columns:
1.  **Column 1:** Path to the query structure (PDB/mmCIF).
2.  **Column 2:** Comma-separated list of motif residues.
3.  **Column 3:** (Optional) path to the output file (default: `stdout`).

```bash
# Search a zinc finger motif against pre-downloaded human proteome (see Download pre-build database)
folddisco query -i index/h_sapiens_folddisco -q query/serine_peptidase.txt
```

## Commands

### Usage of Query Module
```bash
folddisco query -i <INDEX> -p <QUERY_PDB> [-q <QUERY_RESIDUES> -d <DISTANCE_THRESHOLD> -a <ANGLE_THRESHOLD> --skip-match -t <THREADS>]
```

**Important parameter:**
- `-d`: Distance tolerance in Å, increase sensitivity during the prefilter (default: 0.5)
- `-a`: Angle tolerance in degrees, increase sensitivity during the prefilter (default: 5)
- `--nonrigid`: Preset for deformed motifs (see [Non-rigid search](#non-rigid-search))
- `--expand-radius`: How many geometric features may fall in a neighbouring bin at once (default: 1)
- `--novelty-mode`: One verdict line per query instead of a hit list (see [Novelty screening](#novelty-screening))
- `--skip-match`: Skips residue matching and RMSD calculation (prefilter only, much faster with same ranking)
- `--top`: Only report top N hits from the prefilter (controls speed and size of result)
- `-t`: Threads used for search
- `-v`: Verbose output

#### Example Querying
```bash
# Search with default settings. This will print out matching motifs with sorting by RMSD.
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6
folddisco query -p query/1G2F.pdb -q F207,F212,F225,F229 -i index/h_sapiens_folddisco -d 0.5 -a 5 -t 6
folddisco query -p query/1LAP.pdb -q 250,255,273,332,334 -i index/h_sapiens_folddisco --skip-match -t 6 # Skip residue matching

# Query file given as separate text file
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 -d 0.5 -a 5

# Querying a whole structure
folddisco query -i index/h_sapiens_folddisco -p query/1G2F.pdb -t 6 --skip-match
# For a long query, low `--sampling-ratio` can be used to speed up the search
folddisco query -i index/h_sapiens_folddisco -p query/1G2F.pdb -t 6  --skip-match --sampling-ratio 0.3

# Using a query file with distance and angle thresholds
folddisco query -i index/h_sapiens_folddisco -q query/knottin.txt -d 0.5 -a 5 --skip-match -t 6

# Query with amino-acid substitutions and range. 
# Alternative amino acids can be given after colon. 
# X: substitute to any amino acid, p: positive-charged, n: negative-charged, h: hydrophilic, b: hydrophobic, a: aromatic
# Here's enolase query with 3 substitutions; Allow His at 164, Asp & Asn at 247, and His at 297. (Install e_coli_folddisco index first)
folddisco query -p query/2MNR.pdb -q 164:H,195,221,247:ND,297:H -i index/e_coli_folddisco -d 0.5 -a 5 --top 10 --header --per-structure
# Range can be given with dash. This will query first 10 residues and 11th residue with subsitution to any amino acid.
folddisco query -p query/4CHA.pdb -q 1-10,11:X -i index/h_sapiens_folddisco -t 6 --serial-index

# Advanced query with filtering and sorting
## Based on connected node and rmsd
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --connected-node 0.75 --rmsd 1.0

## Coverage based filtering & top N filtering without residue matching
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --covered-node 3 --top 1000 --per-structure --skip-match

# Print top 100 structures with sorting by score
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --top 100 --per-structure --sort-by idf
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --covered-node 4 --top 100 --sort-by idf --per-structure --skip-match

# Comprehensive filtering with multiple criteria
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 -d 0.5 -a 10.0 --ca-distance 1.0 --covered-node-ratio 0.3 --max-node-ratio 0.35 --rmsd 5.0 --tm-score 0.2 --gdt-ts 0.25 --gdt-ha 0.15 --chamfer 5.5 --hausdorff 12.0 --sort-by node_count,gdt_ts,rmsd,idf --format-output tid,node_count,gdt_ts,rmsd,idf,matching_residues,query_residues
```

### Non-rigid search

Motifs are rarely rigid. The same catalytic site in two homologs, or the same site
before and after a conformational change, keeps its residues in the same arrangement
while the distances and angles between them drift by a few tenths of an Ångström —
and a residue pair whose geometry lands on the far side of a bin boundary produces a
different hash and is missed.

| flag | what it does |
| --- | --- |
| `--nonrigid` | Preset: `--expand-radius 2`. The flag to reach for |
| `--expand-radius <INT>` | How many of the geometric features of a residue pair may sit in a neighbouring bin *at the same time*. The default of 1 searches one feature at a time and misses a pair whose distance and angle both drift across a boundary; 2 covers those. 0 searches the observed bins only |

```bash
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --nonrigid
```

**What it buys**, measured on the human proteome index with the benchmark protocol,
answer sets and metrics used by this project (full tables, robustness analysis and
runtime in [feature_evaluation.md](feature_evaluation.md)). F1 over the whole result
list, against the same command without the flag:

| query | F1 without | F1 with `--nonrigid` |
| --- | --- | --- |
| 4-residue zinc finger, matched | 0.9421 | **0.9641** |
| 3-residue zinc finger, matched | 0.9418 | **0.9577** |
| Ser-His-Asp triad, prefilter (independent MEROPS S01 set) | 0.8831 | **0.9160** |
| 23-residue two-segment query, matched | **0.9204** | 0.9117 |

Two conditions come with it, both measured:

- **Use it with `--max-node <n_residues>`.** It raises recall and lowers precision, and
  it is `--max-node` — requiring the whole motif to be covered inside one structure —
  that converts that trade into a win. On the 4-residue query the F1 delta runs from
  −0.053 with no filters, to −0.004 with `--covered-node 3`, to **+0.022** with
  `--covered-node 3 --max-node 4 --rmsd 1.0`.
- **Not for long segment queries.** On the 23-residue query above it loses 0.009 F1
  consistently: long queries are already saturated and the extra candidates only dilute
  the list.

Residue matching runs about 11% longer with the flag on that command (157 → 175 ms), and
about 73% longer on the M-CSA benchmark below. The cost tracks how much the expansion
inflates the candidate pool that then has to be matched, so it grows with motif size and
index size — expect anything in that range. The prefilter difference is smaller than the
run-to-run spread and is not worth quoting.

It is not a repackaging of `-d`/`-a`: across a grid of wider tolerances, widening never
reaches `--nonrigid`'s F1 and *loses* recall rather than gaining it (0.920 → 0.894 against
`--nonrigid`'s 0.970 on the matched 4-residue query). `-d`/`-a` set how far one geometric
feature may move; `--nonrigid` sets how many may move at once, and a deformed motif is
usually two features crossing a bin boundary together.

On a second benchmark — **M-CSA catalytic sites over a 62,122-entry PDB index**, 250
queries, the paper's own mean-sensitivity-at-first-false-positive metric — the picture is
more mixed: it improves 94 queries and degrades 42, raising the mean (0.4475 → 0.4536)
while *lowering* the median (0.4032 → 0.3810). That protocol has no `--max-node`, so it
is the flag measured outside the configuration that makes it work. On the same benchmark,
widening tolerances under a fixed `--top` budget is far worse (0.3599), and stacking
widening with `--nonrigid` triples the queries that return no true positive at all.

Rank the results with **dRMSD** rather than RMSD when the motif may be deformed. dRMSD
compares the internal distances of the match instead of superposing it, so a motif whose
halves swung apart on a hinge keeps a low dRMSD where its superposition RMSD is large:

```bash
# Deformation-ranked non-rigid search
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --nonrigid \
  --sort-by node_count,drmsd --format-output tid,node_count,idf,rmsd,drmsd,matching_residues

# Maximum recall, accepting a worse ranking: widen everything and rescore by deformation
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 \
  -d 1.0 -a 10 --expand-radius 2 --ca-distance 2.0 --sort-by drmsd

# Per-structure output reports the deformation of the best match as min_drmsd
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --nonrigid \
  --per-structure --sort-by max_node_count,min_drmsd --format-output tid,idf,max_node_cov,min_rmsd,min_drmsd
```

`drmsd` and `max_dist_deviation` (the worst single internal-distance change) are
available to `--format-output`, `--sort-by` and `--drmsd`, and work in every output
mode. Both are read straight off the matched coordinates, so no superposition is
involved. They are filtering and reporting metrics: sorting by dRMSD did not beat
sorting by RMSD on the benchmark.

### Novelty screening

`--novelty-mode` replaces the hit list with a single verdict line per query, so a batch
of designed motifs can be screened against a reference database in one pass:

```bash
folddisco query -p design.pdb -q A10,A20,A30 -i index/pdb_folddisco --novelty-mode
# design.pdb	KNOWN	4cha.pdb	1.0000	0.0599	A10,A20,A30
```

Columns are `query_id`, verdict, best hit, residue coverage, best-hit RMSD, query
residues. The verdict is one of:

| verdict | meaning |
| --- | --- |
| `KNOWN` | Coverage ≥ `--novelty-coverage` (default 0.8) **and** best-hit RMSD ≤ `--novelty-rmsd` (default 2.0 Å) |
| `PARTIAL_MATCH` | Something matched, but it failed the coverage or the RMSD bar |
| `NOVEL` | Nothing in the reference database covered a single residue of the motif |
| `FILTERED_OUT` | The database had candidates and this run's search filters kept none of them. The coverage column reports what the database held; there is no hit and no RMSD. Also warned about on stderr |
| `NO_HASHES` | The query could not be searched at all — its residues are further apart than the index distance cutoff, so it produced no hashes. Also warned about on stderr |

The coverage column is not one quantity down the whole file: on `KNOWN` and
`PARTIAL_MATCH` rows it is geometric match coverage, on `FILTERED_OUT` rows the hash-level
coverage the index held, which is generally larger. And a residue named twice in the query
counts twice in the denominator — `A1-A5,A3-A7` has seven distinct residues but ten
entries, so a perfect match reads 0.70 — which is warned about on stderr but not yet
deduplicated, because the query length also feeds `--covered-node-ratio` and
`--max-node-ratio`.

Both bars matter. Coverage alone is not enough: the same residues in a different
arrangement cover everything and are still a different motif — over 100 arbitrary
motifs, a third of the KNOWN verdicts had a best hit worse than 2.0 Å before the RMSD
bar was added. Under `--skip-match` no RMSD is computed, so the verdict falls back to
coverage alone and is weaker. Verdicts append to `-o`, so a batch of queries sharing one
output file accumulates, while rerunning the same command replaces it.

**Screen with the prefilter, or with a low `--max-node`.** `--max-node <n>` drops
structures whose best match covers fewer than *n* residues, so on a 3-residue motif
`--max-node 3` throws away every partial match — including, measured against the PDB
index, a 2-of-3 match at 0.02 Å to the query's own family. A high `--max-node` answers
"is there a full-coverage match"; a novelty screen is asking "does anything like this
exist", and partial matches are most of the answer. Note this pulls the opposite way from
the `--max-node` advice for [`--nonrigid`](#non-rigid-search), which is about ranking
precision: the two recommendations serve different questions and should not be stacked
without thinking about which one you are asking.

### Lookup cache

The first time an index is loaded, Folddisco writes a binary cache of its parsed lookup
file next to it as `<index>.lookup.cache`, and every later load decodes that instead of
re-parsing the text. There is no flag and nothing to manage:

- **On the largest index that ships — the PDB index, 230,655 entries — it saves about
  10 ms** (0.12 s of text parsing against 0.11 s of decoding). That is ~8% of a
  prefilter query on that index and negligible against a matched one, which takes
  ~11 s. A real saving, and a small one.
- It scales with the file: 2.2x at one thread on a 10^6-entry lookup and 2.25x at 10^7,
  which is AFDB-v6 territory and roughly 43x larger than any index that exists today.
  Read those as a projection, not as something you get on a current index. It is
  **never slower** at any size or thread count.
- The first load pays for the write, roughly one extra parse, so a fresh index breaks
  even after about four loads at 8 threads.
- The cache file is ~1.3x the size of the text lookup, and is safe to delete at any time.
- It is invalidated automatically when the lookup file changes: its length, modification
  time and internal sizes are all validated on every load, and anything unexpected falls
  back to parsing the text. The one case it cannot detect is a rewrite that keeps the
  same length *and* restores the original mtime — there is no content hash, because
  hashing hundreds of megabytes would cost more than the load it protects. Delete the
  `.cache` file if you have done that.

### Indexing

### Usage of Index Module
```bash
folddisco index -p <PDB_DIR|FOLDCOMP_DB> -i <INDEX_PATH> -t <THREADS> [-d <DISTANCE_BINS> -a <ANGLE_BINS> -y <FEATURE_TYPE>]
```

**Important parameter:**
- `-d`: Distance threshold in Å for pairs to be included (default: 16)
- `-a`: Bin size of Angle (default: 4)
- `-m`: For big databases (>65k structures) enable -m big for efficiency. Mode `big`, generates an 8GB fixed-size offset.
- `-t`: Threads used for search
- `-v`: Verbose output
- `--type`: Define which features sets are stored in the index; `default` (Folddisco), `pdb` (RCSB feature sets), or `tr` (trRosetta).

#### Examples
```bash
# Default indexing for a small dataset
# h_sapiens directory or foldcomp database is indexed with default parameters
folddisco index -p h_sapiens -i index/h_sapiens_folddisco -t 12

# Indexing big protein dataset
folddisco index -p swissprot -i index/swissprot_folddisco -t 64 -m big -v

# Indexing with custom hash type and parameters
folddisco index -p h_sapiens -i index/h_sapiens_folddisco -t 12 --type default -d 16 -a 4 # Default
folddisco index -p h_sapiens -i index/h_sapiens_pdbtype -t 12 --type pdb -d 8 -a 3 # PDB
```

## Output
### Match Result
Default output which prints out one matching motif per line
```
tid	node_count	idf	rmsd	matching_residues	query_residues
data/serine_peptidases/4cha.pdb	3	8.7616	0.0000	B57,B102,C195	B57,B102,C195
data/serine_peptidases/4cha.pdb	3	8.7616	0.0874	F57,F102,G195	B57,B102,C195
data/serine_peptidases/1pq5.pdb	3	4.1178	0.2609	A56,A99,A195	B57,B102,C195
data/serine_peptidases/1ju3.pdb	2	1.4739	0.7792	_,A223,A234	B57,B102,C195
data/serine_peptidases/1l7a.pdb	2	1.4739	0.7883	_,A146,A127	B57,B102,C195
data/serine_peptidases/1l7a.pdb	2	1.4739	0.8078	_,B146,B127	B57,B102,C195
data/serine_peptidases/1azw.pdb	2	4.6439	0.9234	A179,_,B176	B57,B102,C195
```
- `tid`: Identifier of the target protein structure
- `node_count`: Number of nodes in the match
- `idf`: Inverse document frequency score of matched structure
- `rmsd`: Root mean square deviation
- `matching_residues`: Residue indices in the match (comma-separated, _ for no match)
- `query_residues`: Residue indices in the query (comma-separated)

### Structure Result
Output with one structure per line (`--per-structure`)
```
tid	idf	total_match_count	node_count	edge_count	max_node_cov	min_rmsd	nres	plddt	matching_residues	db_key	query_residues
data/serine_peptidases/4cha.pdb	0.6138	8	3	6	3	0.0000	477	13.5404	B57,B102,C195:0.0000;F57,F102,G195:0.0874	4	B57,B102,C195
data/serine_peptidases/1pq5.pdb	0.4869	4	3	4	3	0.2609	224	5.1340	A56,A99,A195:0.2609	3	B57,B102,C195
data/serine_peptidases/1ju3.pdb	0.0617	2	2	2	2	0.7792	570	19.4881	_,A223,A234:0.7792	1	B57,B102,C195
data/serine_peptidases/1l7a.pdb	0.0584	2	2	2	2	0.7883	636	11.7037	_,A146,A127:0.7883;_,B146,B127:0.8078	2	B57,B102,C195
data/serine_peptidases/1azw.pdb	0.1856	2	2	2	2	0.9234	626	34.2399	A179,_,B176:0.9234	0	B57,B102,C195
```
- `tid`: Identifier of the target protein structure
- `idf`: Inverse document frequency score with length penalty; Higher score indicates more matches within smaller structures
- `total_match_count`: Total number of matches
- `node_count`: Number of nodes in the structure
- `edge_count`: Number of edges in the structure
- `max_node_cov`: Maximum node coverage
- `min_rmsd`: Minimum root mean square deviation
- `nres`: Number of residues
- `plddt`: Predicted local distance difference test score
- `matching_residues`: Residue indices in the match (comma-separated, _ for no match, semicolon-separated for multiple matches with RMSD)
- `key`: Numeric identifier of the structure
- `query_residues`: Residue indices in the query (comma-separated)

### Display Options
- `--per-structure`: Outputs results per structure.
- `--per-match`: Outputs results per match.
- `--sort-by`: Sorts results by given columns (comma-separated).
- `--format-output`: Custom output format using column names.
- `--top <N>`: Outputs top N results.
- `--header`: Outputs header for the result.

## Contributions

<a href="https://github.com/steineggerlab/folddisco/graphs/contributors">
  <img src="https://contributors-img.firebaseapp.com/image?repo=steineggerlab/folddisco" />
</a>
