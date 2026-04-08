//! # Criterion Benchmarks for the GRAMS Metric
//!
//! Run with:
//!
//! ```text
//! cargo bench --bench metric_benchmark
//! ```
//!
//! Benchmark groups:
//!
//! | Group | What is measured |
//! |-------|-----------------|
//! | `grams_score_by_size` | End-to-end GRAMS computation for motifs of 3–10 residues |
//! | `grams_subcomponents` | Individual sub-scores (TM, DMS, PAS) for a 6-residue motif |
//! | `grams_10k_batch` | Throughput: 10 000 randomised 6-residue motif comparisons |

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use folddisco::ranking::{
    distance_matrix_score, grams_score, grams_score_detailed, pseudo_bond_angle_score,
    GramsWeights,
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Simple linear congruential step — deterministic, no dependency needed.
#[inline]
fn pseudo_rand(seed: u64, lo: f32, hi: f32) -> f32 {
    let x = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let frac = ((x >> 33) as f32) / (u32::MAX as f32);
    lo + frac * (hi - lo)
}

fn random_coords(n: usize, seed_base: u64) -> Vec<[f32; 3]> {
    (0..n)
        .map(|i| {
            let s = seed_base.wrapping_add(i as u64 * 3);
            [
                pseudo_rand(s, -20.0, 20.0),
                pseudo_rand(s.wrapping_add(1), -20.0, 20.0),
                pseudo_rand(s.wrapping_add(2), -20.0, 20.0),
            ]
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Benchmark: end-to-end GRAMS by motif size
// ---------------------------------------------------------------------------

fn bench_grams_by_size(c: &mut Criterion) {
    let mut group = c.benchmark_group("grams_score_by_size");

    for n in [3usize, 4, 5, 6, 7, 8, 9, 10] {
        let ref_coords = random_coords(n, 42);
        let model_coords = random_coords(n, 99);
        let weights = GramsWeights::default();

        group.throughput(Throughput::Elements(1));
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| {
                grams_score(
                    black_box(&ref_coords),
                    black_box(&model_coords),
                    black_box(weights),
                )
            })
        });
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// Benchmark: individual sub-scores for a 6-residue motif
// ---------------------------------------------------------------------------

fn bench_subcomponents(c: &mut Criterion) {
    let n = 6;
    let ref_coords = random_coords(n, 42);
    let model_coords = random_coords(n, 99);
    let weights = GramsWeights::default();

    let mut group = c.benchmark_group("grams_subcomponents");

    group.bench_function("distance_matrix_score", |b| {
        b.iter(|| distance_matrix_score(black_box(&ref_coords), black_box(&model_coords)))
    });

    group.bench_function("pseudo_bond_angle_score", |b| {
        b.iter(|| pseudo_bond_angle_score(black_box(&ref_coords), black_box(&model_coords)))
    });

    group.bench_function("grams_score_detailed", |b| {
        b.iter(|| {
            grams_score_detailed(
                black_box(&ref_coords),
                black_box(&model_coords),
                black_box(weights),
            )
        })
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Benchmark: 10 000 randomised comparisons (throughput test)
// ---------------------------------------------------------------------------

fn bench_10k_batch(c: &mut Criterion) {
    const BATCH: usize = 10_000;
    const N: usize = 6;

    let ref_coords: Vec<Vec<[f32; 3]>> = (0..BATCH)
        .map(|i| random_coords(N, i as u64 * 7 + 1))
        .collect();
    let model_coords: Vec<Vec<[f32; 3]>> = (0..BATCH)
        .map(|i| random_coords(N, i as u64 * 13 + 5))
        .collect();
    let weights = GramsWeights::default();

    let mut group = c.benchmark_group("grams_10k_batch");
    group.throughput(Throughput::Elements(BATCH as u64));

    group.bench_function("10k_6-residue_motifs", |b| {
        b.iter(|| {
            let mut total = 0.0_f32;
            for i in 0..BATCH {
                total += grams_score(
                    black_box(&ref_coords[i]),
                    black_box(&model_coords[i]),
                    black_box(weights),
                );
            }
            black_box(total)
        })
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Register groups
// ---------------------------------------------------------------------------

criterion_group!(
    benches,
    bench_grams_by_size,
    bench_subcomponents,
    bench_10k_batch
);
criterion_main!(benches);
