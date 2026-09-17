# Changelog

## 3.0.0

Indices built with 2.x are read unchanged.

### New
- `query --sensitive` (= `--expand-radius 2`) and `--expand-radius`: residue pairs may fall
  in neighbouring distance/angle bins, several features at once. `-d`/`-a` are sub-stepped.
- `query --aa-subst blosum62|group|size` and per-residue `:*`: substitution schemes.
  Substituted residues score below exact ones, so exact matches rank first.
- `index --expand-radius/--expand-distance/--expand-angle/--aa-subst`: index-time expansion.
- `query --novelty-mode`: one evidence row per query (coverage, best hit, RMSD).
- Multi-character and numeric chain IDs (`AA_250`, `10_250`); `--chain-sep`.
- Superposition-free `drmsd` and `max_dist_deviation` (columns, sort keys, `--drmsd`).
- Memory-mapped caches for index and Foldcomp lookups (`*.lookup.cache`, `*.fdcache`),
  built on first use.

### Changed
- Matching prefilters target residues by amino acid code for large queries (same output, faster).
- Per-match IDF counts edges from neighbouring bins and substitutions only between the matched
  residues, so large components of similar residues no longer rank first.
- A residue that ends a chain is no longer labelled with the next chain's ID.
- `-d`/`-a` with several values use only the widest one.
- `-q` ranges may repeat the chain (`F204-F215`); malformed queries are reported, not panicked on.
- Help states the actual default sort, `idf:desc,rmsd:asc`.

### Fixed
- Substituted residues (`:H`) were found by the index but dropped at matching.
