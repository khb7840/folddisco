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
* **Residues:** `B57` = chain `B`, residue number `57`. Ranges are inclusive: `1-10`.
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
- `--dist-ratio`: Distance tolerance proportional to the pair distance (default: 0)
- `--expand-radius`: How many geometric features may fall in a neighbouring bin at once (default: 1)
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
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --top 100 --per-structure --sort-by-score
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --covered-node 4 --top 100 --sort-by-score --per-structure --skip-match

# Comprehensive filtering with multiple criteria
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 -d 0.5 -a 10.0 --ca-distance 1.0 --covered-node-ratio 0.3 --max-node-ratio 0.35 --rmsd 5.0 --tm-score 0.2 --gdt-ts 0.25 --gdt-ha 0.15 --chamfer-distance 5.5 --hausdorff-distance 12.0 --sort-by node_count,gdt_ts,rmsd,idf --format-output tid,node_count,gdt_ts,rmsd,idf,matching_residues,query_residues
```

### Non-rigid search

Design notes, the full benchmark and the alternatives that were measured and rejected
are in [nonrigid_search.md](nonrigid_search.md).

Motifs are rarely rigid. The same catalytic site in two homologs, or the same site
before and after a conformational change, keeps its residues in the same arrangement
while the distances and angles between them drift by a few tenths of an Ångström —
and a residue pair whose geometry lands on the far side of a bin boundary produces a
different hash and is missed. Three flags widen the search for that case:

| flag | what it does |
| --- | --- |
| `--nonrigid` | Preset: `--expand-radius 2`. An explicit `--expand-radius` is kept when it is looser |
| `--expand-radius <INT>` | How many of the geometric features of a residue pair may sit in a neighbouring bin *at the same time*. The default of 1 searches one feature at a time and misses a pair whose distance and angle both drift across a boundary; 2 covers those. 0 searches the observed bins only |
| `--enm-sample` | Wiggle the query along its low-frequency **torsional normal modes** and search the union of the ensemble's hashes. Tune with `--num-confs`, `--nma-rmsd`, `--nma-modes`. Best deep recall measured, at ~4x runtime |
| `--dist-ratio <FLOAT>` | Adds a tolerance proportional to the pair distance, on top of `-d`, on the theory that elastic deformation scales with distance. **Measured as a recall-for-ranking trade** (see below), so it is off by default and not part of `--nonrigid` |

Independently of the flags, a tolerance-expanded hash that turns out far rarer in the
database than the observed hash it came from is now dropped during search. A rare hash
carries a large IDF, so one spurious hit on such a hash can outrank several real ones.
The threshold is a ratio between two IDFs from the same index, so it needs no tuning
per database.

#### Benchmark

Zinc-finger motif from `1G2F` against the human proteome index (23,391 proteins,
1,816 annotated zinc-finger answers in `data/zinc_answer.tsv`). `TP@k FP` is the
number of true positives found walking the ranked list until *k* false positives:

| motif | setting | hits | TP@5FP | TP@10FP | TP@100FP | recall | sec |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `F207,F212,F225,F229` | before these changes | 1645 | 464 | 649 | 782 | 0.496 | 0.011 |
| | default | 1607 | 578 | 691 | 789 | 0.495 | 0.010 |
| | `--nonrigid` | 1850 | 684 | 731 | 795 | 0.518 | 0.011 |
| | `--nonrigid --enm-sample --num-confs 10` | 2021 | **709** | **746** | **823** | **0.542** | 0.040 |
| `F207,F212,F225` | before these changes | 595 | 94 | 113 | 268 | 0.197 | 0.009 |
| | default | 590 | 94 | 113 | 269 | 0.197 | 0.009 |
| | `--nonrigid` | 1248 | **108** | **145** | 301 | 0.484 | 0.009 |
| | `--nonrigid --enm-sample --num-confs 10` | 1417 | 98 | 118 | **304** | **0.509** | 0.056 |

The rare-hash filter alone lifts the 4-residue motif by 25% at 5 false positives with
no runtime cost. `--nonrigid` adds most for short motifs. `--enm-sample` gives the best
deep recall but can lose ground at the very top of the ranking for a 3-residue motif,
which is why it is opt-in. `--expand-radius 3` measured no better than 2.

`--dist-ratio` behaves differently: it raises total recall but loses true positives at
every early-precision point, because widening the tolerance most for long pairs adds
exactly the least specific matches. Reach for it only when total recall is what
matters and you intend to rescore the hits yourself.

Rank the results with **dRMSD** rather than RMSD. dRMSD compares the internal
distances of the match instead of superposing it, so a motif whose halves swung apart
on a hinge keeps a low dRMSD where its superposition RMSD is large:

```bash
# Deformation-ranked non-rigid search
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --nonrigid \
  --sort-by node_count,drmsd --format-output tid,node_count,idf,rmsd,drmsd,matching_residues

# Maximum recall, accepting a worse ranking: widen everything and rescore by deformation
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 \
  -d 1.0 -a 10 --dist-ratio 0.08 --expand-radius 2 --ca-distance 2.0 --sort-by drmsd

# Per-structure output reports the deformation of the best match as min_drmsd
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --nonrigid \
  --per-structure --sort-by max_node_count,min_drmsd --format-output tid,idf,max_node_cov,min_rmsd,min_drmsd
```

`drmsd` and `max_dist_deviation` (the worst single internal-distance change) are
available to `--format-output`, `--sort-by` and `--drmsd`, and work in every output
mode. Both are read straight off the matched coordinates, so no superposition is
involved.

#### Runtime

Same index, 8 threads, median of 15 runs, warm cache. `--nonrigid` roughly doubles the
number of query hashes, which shows up as single-digit percent for a motif query and
more for a whole-structure query:

| query | prefilter (`--skip-match`) | with residue matching |
| --- | --- | --- |
| 3–8 residue motif | 8–30 ms, `--nonrigid` +3–5% | 39–44 ms, `--nonrigid` +5–6% |
| 16 residue motif | 43 ms, `--nonrigid` +4% | 82 ms, `--nonrigid` +6% |
| whole structure (1G2F) | 129 ms, `--nonrigid` +18% | 302 ms, `--nonrigid` +14% |
| whole structure (1LAP) | 510 ms, `--nonrigid` +37% | 6.2 s, `--nonrigid` ~0% |

For a motif query — the usual case — the sensitivity gain costs a few percent. Pair
`--nonrigid` with `--top`, `--sampling-ratio` or `--skip-match` for whole-structure
queries, where the query is already thousands of hashes before any expansion.

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
