// File: novelty.rs
// Description: Novelty checking for structural motifs.
//   Classifies a query motif as NOVEL, PARTIAL_MATCH, or KNOWN based on
//   how well known structures cover the query residues.

use std::fmt;

use crate::controller::result::{MatchResult, StructureResult, evalue_fitting};

// ---------------------------------------------------------------------------
// Formatting helpers
// ---------------------------------------------------------------------------

/// Format an f32 metric as a string, returning "NA" when the value equals the
/// sentinel (typically `f32::MAX` meaning "not computed").
#[inline]
fn format_optional_f32(value: f32, sentinel: f32, fmt_str: &str) -> String {
    if value == sentinel {
        "NA".to_string()
    } else if fmt_str == "{:.4}" {
        format!("{:.4}", value)
    } else {
        format!("{:.4e}", value)
    }
}

/// Format an f64 metric as a string, returning "NA" when the value equals the
/// sentinel (typically `f64::MAX` meaning "not computed").
#[inline]
fn format_optional_f64(value: f64, sentinel: f64, fmt_str: &str) -> String {
    if value == sentinel {
        "NA".to_string()
    } else if fmt_str == "{:.4}" {
        format!("{:.4}", value)
    } else {
        format!("{:.4e}", value)
    }
}

/// Three-tier classification of motif novelty.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NoveltyTier {
    /// No structure in the reference database covers the motif above the
    /// coverage threshold – the motif is genuinely novel.
    Novel,
    /// At least one structure partially covers the motif but none reaches
    /// the full-coverage threshold.
    PartialMatch,
    /// At least one structure covers the motif at or above the coverage
    /// threshold – the motif is known.
    Known,
}

impl fmt::Display for NoveltyTier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NoveltyTier::Novel => write!(f, "NOVEL"),
            NoveltyTier::PartialMatch => write!(f, "PARTIAL_MATCH"),
            NoveltyTier::Known => write!(f, "KNOWN"),
        }
    }
}

/// Summary of sub-motif (pair/triplet) overlap with the reference database.
#[derive(Debug, Clone)]
pub struct SubMotifSummary {
    /// Number of residue pairs present in at least one reference structure.
    pub known_pairs: usize,
    /// Total residue pairs in the query motif.
    pub total_pairs: usize,
}

impl SubMotifSummary {
    pub fn new(known_pairs: usize, total_pairs: usize) -> Self {
        Self { known_pairs, total_pairs }
    }
}

impl fmt::Display for SubMotifSummary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}/{}_pairs_known", self.known_pairs, self.total_pairs)
    }
}

/// Per-design novelty report produced by `folddisco novelty`.
#[derive(Debug, Clone)]
pub struct NoveltyReport {
    /// Identifier for the design (typically the PDB file path or basename).
    pub design_id: String,
    /// Query residue specification string.
    pub query_residues: String,
    /// Novelty classification.
    pub novelty_tier: NoveltyTier,
    /// Total number of residues in the query motif.
    pub query_residue_count: usize,
    /// Target ID of the best-coverage hit, if any.
    pub best_hit: Option<String>,
    /// Fraction of query residues covered by the best hit (0.0–1.0).
    pub best_hit_coverage: f32,
    /// RMSD of the best hit match, or f32::MAX when there is no match.
    pub best_hit_rmsd: f32,
    /// IDF score of the best hit, or 0.0 when there is no match.
    pub best_hit_idf: f32,
    /// E-value of the best hit match.
    pub best_hit_evalue: f64,
    /// Summary of sub-motif (pair-level) overlap, if computed.
    pub sub_motif_summary: Option<SubMotifSummary>,
    /// Name of the reference database index that produced the best hit.
    pub source_index: Option<String>,
}

impl NoveltyReport {
    /// Build a report from the output of the count+retrieve pipeline.
    ///
    /// # Arguments
    /// * `design_id`          – identifier for this design
    /// * `query_residues`     – residue query string
    /// * `query_residue_count`– number of residues in the query
    /// * `results`            – ranked `StructureResult` list (may be empty)
    /// * `coverage_threshold` – minimum coverage fraction to call a motif KNOWN
    /// * `index_size`         – number of structures in the reference database
    /// * `source_index`       – name/path of the reference index (optional)
    pub fn from_structure_results<'a>(
        design_id: String,
        query_residues: String,
        query_residue_count: usize,
        results: &[(usize, StructureResult<'a>)],
        coverage_threshold: f32,
        index_size: usize,
        source_index: Option<String>,
    ) -> Self {
        if results.is_empty() || query_residue_count == 0 {
            return Self {
                design_id,
                query_residues,
                novelty_tier: NoveltyTier::Novel,
                query_residue_count,
                best_hit: None,
                best_hit_coverage: 0.0,
                best_hit_rmsd: f32::MAX,
                best_hit_idf: 0.0,
                best_hit_evalue: f64::MAX,
                sub_motif_summary: None,
                source_index,
            };
        }

        // Find the result with the highest node (residue) coverage.
        // When coverage is equal, prefer lower RMSD (better geometric match),
        // so we reverse the RMSD ordering in the tiebreaker.
        let best = results.iter().max_by(|a, b| {
            let cov_a = a.1.max_matching_node_count;
            let cov_b = b.1.max_matching_node_count;
            cov_a.cmp(&cov_b)
                .then_with(|| a.1.min_rmsd_with_max_match.partial_cmp(&b.1.min_rmsd_with_max_match)
                    .map(|ord| ord.reverse()) // lower RMSD is better → reverse for max_by
                    .unwrap_or(std::cmp::Ordering::Equal))
        });

        match best {
            None => Self {
                design_id,
                query_residues,
                novelty_tier: NoveltyTier::Novel,
                query_residue_count,
                best_hit: None,
                best_hit_coverage: 0.0,
                best_hit_rmsd: f32::MAX,
                best_hit_idf: 0.0,
                best_hit_evalue: f64::MAX,
                sub_motif_summary: None,
                source_index,
            },
            Some((_, r)) => {
                let coverage = r.max_matching_node_count as f32 / query_residue_count as f32;
                let novelty_tier = if coverage >= coverage_threshold {
                    NoveltyTier::Known
                } else if r.max_matching_node_count > 0 {
                    NoveltyTier::PartialMatch
                } else {
                    NoveltyTier::Novel
                };
                let evalue = evalue_fitting(r.idf, index_size as f32, query_residue_count as f32);
                Self {
                    design_id,
                    query_residues,
                    novelty_tier,
                    query_residue_count,
                    best_hit: Some(r.tid.to_string()),
                    best_hit_coverage: coverage,
                    best_hit_rmsd: r.min_rmsd_with_max_match,
                    best_hit_idf: r.idf,
                    best_hit_evalue: evalue,
                    sub_motif_summary: None,
                    source_index,
                }
            }
        }
    }

    /// Build a report when only the prefilter (skip-match) results are available.
    /// In that case `max_matching_node_count` is zero, so we fall back to
    /// `node_count` (covered nodes from the inverted-index step).
    pub fn from_structure_results_skip_match<'a>(
        design_id: String,
        query_residues: String,
        query_residue_count: usize,
        results: &[(usize, StructureResult<'a>)],
        coverage_threshold: f32,
        index_size: usize,
        source_index: Option<String>,
    ) -> Self {
        if results.is_empty() || query_residue_count == 0 {
            return Self {
                design_id,
                query_residues,
                novelty_tier: NoveltyTier::Novel,
                query_residue_count,
                best_hit: None,
                best_hit_coverage: 0.0,
                best_hit_rmsd: f32::MAX,
                best_hit_idf: 0.0,
                best_hit_evalue: f64::MAX,
                sub_motif_summary: None,
                source_index,
            };
        }

        // Use node_count (inverted-index coverage) as best proxy
        let best = results.iter().max_by(|a, b| {
            a.1.node_count.cmp(&b.1.node_count)
                .then_with(|| a.1.idf.partial_cmp(&b.1.idf).unwrap_or(std::cmp::Ordering::Equal))
        });

        match best {
            None => Self {
                design_id,
                query_residues,
                novelty_tier: NoveltyTier::Novel,
                query_residue_count,
                best_hit: None,
                best_hit_coverage: 0.0,
                best_hit_rmsd: f32::MAX,
                best_hit_idf: 0.0,
                best_hit_evalue: f64::MAX,
                sub_motif_summary: None,
                source_index,
            },
            Some((_, r)) => {
                let coverage = r.node_count as f32 / query_residue_count as f32;
                let novelty_tier = if coverage >= coverage_threshold {
                    NoveltyTier::Known
                } else if r.node_count > 0 {
                    NoveltyTier::PartialMatch
                } else {
                    NoveltyTier::Novel
                };
                let evalue = evalue_fitting(r.idf, index_size as f32, query_residue_count as f32);
                Self {
                    design_id,
                    query_residues,
                    novelty_tier,
                    query_residue_count,
                    best_hit: Some(r.tid.to_string()),
                    best_hit_coverage: coverage,
                    best_hit_rmsd: r.min_rmsd_with_max_match,
                    best_hit_idf: r.idf,
                    best_hit_evalue: evalue,
                    sub_motif_summary: None,
                    source_index,
                }
            }
        }
    }

    /// Attach a sub-motif summary to this report.
    pub fn with_sub_motif_summary(mut self, summary: SubMotifSummary) -> Self {
        self.sub_motif_summary = Some(summary);
        self
    }

    /// TSV header line for a batch of novelty reports.
    pub fn tsv_header() -> &'static str {
        "design_id\tquery_residues\tnovelty_tier\tquery_residue_count\tbest_hit\tbest_hit_coverage\tbest_hit_rmsd\tbest_hit_idf\tbest_hit_evalue\tsub_motif_summary\tsource_index"
    }

    /// Format the report as a single TSV line (no trailing newline).
    pub fn to_tsv(&self) -> String {
        let best_hit = self.best_hit.as_deref().unwrap_or("NA");
        let rmsd_str = format_optional_f32(self.best_hit_rmsd, f32::MAX, "{:.4}");
        let evalue_str = format_optional_f64(self.best_hit_evalue, f64::MAX, "{:.4e}");
        let sub_motif = self.sub_motif_summary.as_ref()
            .map(|s| s.to_string())
            .unwrap_or_else(|| "NA".to_string());
        let source_index = self.source_index.as_deref().unwrap_or("NA");
        format!(
            "{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{}\t{}\t{}",
            self.design_id,
            self.query_residues,
            self.novelty_tier,
            self.query_residue_count,
            best_hit,
            self.best_hit_coverage,
            rmsd_str,
            self.best_hit_idf,
            evalue_str,
            sub_motif,
            source_index,
        )
    }
}

/// Determine the novelty tier of a single `MatchResult` (per-match output mode).
///
/// Returns `NoveltyTier::Known` when the match covers at least
/// `coverage_threshold` fraction of the query residues, `PartialMatch` when
/// there is at least one matched residue, and `Novel` otherwise.
pub fn novelty_tier_from_match(
    result: &MatchResult,
    query_residue_count: usize,
    coverage_threshold: f32,
) -> NoveltyTier {
    if query_residue_count == 0 {
        return NoveltyTier::Novel;
    }
    let coverage = result.node_count as f32 / query_residue_count as f32;
    if coverage >= coverage_threshold {
        NoveltyTier::Known
    } else if result.node_count > 0 {
        NoveltyTier::PartialMatch
    } else {
        NoveltyTier::Novel
    }
}

/// Count how many distinct residue pairs from the query appear in at least
/// one result structure (uses the inverted-index edge count as a proxy).
///
/// `total_pairs` = n*(n-1)/2 for a query of n residues.
pub fn count_known_pairs_from_results<'a>(
    results: &[(usize, StructureResult<'a>)],
    query_residue_count: usize,
) -> SubMotifSummary {
    let total_pairs = if query_residue_count >= 2 {
        query_residue_count * (query_residue_count - 1) / 2
    } else {
        0
    };
    if total_pairs == 0 || results.is_empty() {
        return SubMotifSummary::new(0, total_pairs);
    }
    // Use max edge_count across all results as a conservative lower bound on
    // how many distinct pairs are known.
    let max_edges = results.iter().map(|(_, r)| r.edge_count).max().unwrap_or(0);
    let known_pairs = max_edges.min(total_pairs);
    SubMotifSummary::new(known_pairs, total_pairs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_novelty_tier_display() {
        assert_eq!(NoveltyTier::Novel.to_string(), "NOVEL");
        assert_eq!(NoveltyTier::PartialMatch.to_string(), "PARTIAL_MATCH");
        assert_eq!(NoveltyTier::Known.to_string(), "KNOWN");
    }

    #[test]
    fn test_novelty_report_tsv_header() {
        let header = NoveltyReport::tsv_header();
        assert!(header.starts_with("design_id"));
        assert!(header.contains("novelty_tier"));
        assert!(header.contains("best_hit_coverage"));
    }

    #[test]
    fn test_sub_motif_summary_display() {
        let s = SubMotifSummary::new(3, 6);
        assert_eq!(s.to_string(), "3/6_pairs_known");
    }
}
