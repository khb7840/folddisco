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
//! | `grams_subcomponents` | Individual sub-scores (TM, ResComp, ClashFrac) |
//! | `grams_10k_batch` | Throughput: 10 000 randomised 6-residue motif comparisons |

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use folddisco::ranking::{
    clash_fraction, grams_score, grams_score_detailed, residue_compatibility, GramsWeights,
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Pseudo-random f32 from a seed, range [lo, hi).
#[inline]
fn pseudo_rand(seed: u64, lo: f32, hi: f32) -> f32 {
    // xorshift64 step
    let x = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
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

fn random_residues(n: usize, seed: u64) -> Vec<u8> {
    const AAS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    (0..n)
        .map(|i| {
            let x = seed.wrapping_add(i as u64).wrapping_mul(2862933555777941757);
            AAS[(x >> 58) as usize % 20]
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
        let ref_res = random_residues(n, 17);
        let model_res = random_residues(n, 31);
        let weights = GramsWeights::default();

        group.throughput(Throughput::Elements(1));
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| {
                grams_score(
                    black_box(&ref_coords),
                    black_box(&model_coords),
                    black_box(&ref_res),
                    black_box(&model_res),
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
    let ref_res = random_residues(n, 17);
    let model_res = random_residues(n, 31);
    let weights = GramsWeights::default();

    let mut group = c.benchmark_group("grams_subcomponents");

    group.bench_function("residue_compatibility", |b| {
        b.iter(|| residue_compatibility(black_box(&ref_res), black_box(&model_res)))
    });

    group.bench_function("clash_fraction", |b| {
        b.iter(|| clash_fraction(black_box(&model_coords)))
    });

    group.bench_function("grams_score_detailed", |b| {
        b.iter(|| {
            grams_score_detailed(
                black_box(&ref_coords),
                black_box(&model_coords),
                black_box(&ref_res),
                black_box(&model_res),
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

    // Pre-generate all data outside the benchmark loop
    let ref_coords: Vec<Vec<[f32; 3]>> = (0..BATCH)
        .map(|i| random_coords(N, i as u64 * 7 + 1))
        .collect();
    let model_coords: Vec<Vec<[f32; 3]>> = (0..BATCH)
        .map(|i| random_coords(N, i as u64 * 13 + 5))
        .collect();
    let ref_res: Vec<Vec<u8>> = (0..BATCH)
        .map(|i| random_residues(N, i as u64 * 3 + 2))
        .collect();
    let model_res: Vec<Vec<u8>> = (0..BATCH)
        .map(|i| random_residues(N, i as u64 * 11 + 8))
        .collect();
    let weights = GramsWeights::default();

    let mut group = c.benchmark_group("grams_10k_batch");
    group.throughput(Throughput::Elements(BATCH as u64));

    group.bench_function("10k_6-residue_motifs", |b| {
        b.iter(|| {
            let mut _total = 0.0_f32;
            for i in 0..BATCH {
                _total += grams_score(
                    black_box(&ref_coords[i]),
                    black_box(&model_coords[i]),
                    black_box(&ref_res[i]),
                    black_box(&model_res[i]),
                    black_box(weights),
                );
            }
            black_box(_total)
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
