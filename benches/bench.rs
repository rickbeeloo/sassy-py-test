#![feature(portable_simd, array_chunks)]
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use sassy::*;
use std::time::Duration;

fn generate_test_data(size: usize) -> Vec<u8> {
    let mut rng = rand::rng();
    let bases = vec![b'A', b'T', b'G', b'C'];
    let mut seq = vec![b'A'; size];
    for i in 0..size {
        seq[i] = bases[rand::Rng::random_range(&mut rng, 0..4)];
    }
    seq
}

fn benchmark_base_lookup(c: &mut Criterion) {
    let mut group = c.benchmark_group("Base lookup");
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(1));
    group.sample_size(10);

    // Test different sequence sizes, for now just one

    #[allow(clippy::single_element_loop)]
    for size in [1024 * 1024] {
        let seq = generate_test_data(size);

        group.throughput(Throughput::Bytes(size as u64));

        let query_bases_defaults = b"NY";

        group.bench_with_input(BenchmarkId::new("profile u32", size), &seq, |b, seq| {
            b.iter(|| {
                let mut result = vec![0; 6];
                for chunk in seq.array_chunks() {
                    packed_nibbles_portable_32(chunk, black_box(query_bases_defaults), &mut result);
                    black_box(&mut result);
                }
            })
        });

        group.bench_with_input(BenchmarkId::new("profile u64", size), &seq, |b, seq| {
            b.iter(|| {
                let mut result = vec![0; 6];
                for chunk in seq.array_chunks() {
                    packed_nibbles_portable_64(chunk, black_box(query_bases_defaults), &mut result);
                    black_box(&mut result);
                }
            })
        });

        let query = b"ACTGCAACTGCAACGACGTA";
        group.bench_with_input(BenchmarkId::new("search_20", size), &seq, |b, seq| {
            let mut deltas = vec![];
            b.iter(|| {
                sassy::search::<Iupac>(black_box(query), seq, &mut deltas);
                black_box(&deltas);
            })
        });
        let query = b"ACTGCAANTGCAACGACGTA";
        group.bench_with_input(BenchmarkId::new("search_20_N", size), &seq, |b, seq| {
            let mut deltas = vec![];
            b.iter(|| {
                sassy::search::<Iupac>(black_box(query), seq, &mut deltas);
                black_box(&deltas);
            })
        });
        let query = b"ACTGCAACTGCAACGACGTAACACCTACTAAC";
        group.bench_with_input(BenchmarkId::new("search_32", size), &seq, |b, seq| {
            let mut deltas = vec![];
            b.iter(|| {
                sassy::search::<Iupac>(black_box(query), seq, &mut deltas);
                black_box(&deltas);
            })
        });
        let query = b"ACTGCAANTGCAACGAYGTAACARCTACTAAC";
        group.bench_with_input(BenchmarkId::new("search_32_NRY", size), &seq, |b, seq| {
            let mut deltas = vec![];
            b.iter(|| {
                sassy::search::<Iupac>(black_box(query), seq, &mut deltas);
                black_box(&deltas);
            })
        });
        let query = b"ACTGCAACTGCAACGACGTAACACCTACTAAC";
        group.bench_with_input(BenchmarkId::new("find_32", size), &seq, |b, seq| {
            let mut deltas = vec![];
            let mut positions = vec![];
            b.iter(|| {
                sassy::search::<Iupac>(black_box(query), seq, &mut deltas);
                positions.clear();
                sassy::find_below_threshold(black_box(query), 0, &deltas, &mut positions);
                eprintln!("len: {}", positions.len());
            })
        });
    }

    group.finish();
}

criterion_group!(benches, benchmark_base_lookup);
criterion_main!(benches);
