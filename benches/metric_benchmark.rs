//! # Criterion Benchmarks for DMS/PAS/SOS
//!
//! Run with:
//!
//! ```text
//! cargo bench --bench metric_benchmark
//! ```

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use folddisco::structure::metrics::{
    distance_matrix_score, pseudo_bond_angle_score, side_chain_orientation_score,
};

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

fn bench_dms_by_size(c: &mut Criterion) {
    let mut group = c.benchmark_group("dms_by_size");

    for n in [3usize, 4, 5, 6, 7, 8, 9, 10] {
        let ref_coords = random_coords(n, 42);
        let model_coords = random_coords(n, 99);

        group.throughput(Throughput::Elements(1));
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| distance_matrix_score(black_box(&ref_coords), black_box(&model_coords)))
        });
    }

    group.finish();
}

fn bench_pas_by_size(c: &mut Criterion) {
    let mut group = c.benchmark_group("pas_by_size");

    for n in [3usize, 4, 5, 6, 7, 8, 9, 10] {
        let ref_coords = random_coords(n, 42);
        let model_coords = random_coords(n, 99);

        group.throughput(Throughput::Elements(1));
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| pseudo_bond_angle_score(black_box(&ref_coords), black_box(&model_coords)))
        });
    }

    group.finish();
}

fn bench_sos_by_size(c: &mut Criterion) {
    let mut group = c.benchmark_group("sos_by_size");

    for n in [3usize, 4, 5, 6, 7, 8, 9, 10] {
        let ref_ca = random_coords(n, 42);
        let model_ca = random_coords(n, 99);
        let ref_cb = random_coords(n, 7);
        let model_cb = random_coords(n, 13);

        group.throughput(Throughput::Elements(1));
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| {
                side_chain_orientation_score(
                    black_box(&ref_ca),
                    black_box(&model_ca),
                    black_box(&ref_cb),
                    black_box(&model_cb),
                )
            })
        });
    }

    group.finish();
}

criterion_group!(benches, bench_dms_by_size, bench_pas_by_size, bench_sos_by_size);
criterion_main!(benches);
