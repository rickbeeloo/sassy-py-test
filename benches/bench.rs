#![feature(portable_simd, array_chunks)]
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use sassy::profiles::*;
use std::time::Duration;

fn generate_dna_sequence(size: usize) -> Vec<u8> {
    let mut rng = rand::rng();
    let bases = [b'A', b'T', b'G', b'C'];
    let mut seq = vec![b'A'; size];
    for i in 0..size {
        seq[i] = bases[rand::Rng::random_range(&mut rng, 0..4)];
    }
    seq
}

fn generate_ascii_sequence(size: usize) -> Vec<u8> {
    let mut rng = rand::rng();
    let mut seq = vec![0; size];
    for i in 0..size {
        seq[i] = rand::Rng::random_range(&mut rng, 0..256) as u8;
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
        let dna_seq = generate_dna_sequence(size);
        let ascii_seq = generate_ascii_sequence(size);
        group.throughput(Throughput::Bytes(size as u64));

        group.bench_with_input(
            BenchmarkId::new("encode_iupac", size),
            &dna_seq,
            |b, seq| {
                let profiler = Iupac::encode_query(b"NY").0;
                let mut result = profiler.alloc_out();
                b.iter(|| {
                    for chunk in seq.array_chunks() {
                        profiler.encode_ref(chunk, &mut result);
                        black_box(&mut result);
                    }
                })
            },
        );

        group.bench_with_input(BenchmarkId::new("encode_dna", size), &dna_seq, |b, seq| {
            let profiler = Dna::encode_query(b"ACTG").0;
            let mut result = profiler.alloc_out();
            b.iter(|| {
                for chunk in seq.array_chunks() {
                    profiler.encode_ref(chunk, &mut result);
                    black_box(&mut result);
                }
            })
        });

        group.bench_with_input(
            BenchmarkId::new("encode_ascii", size),
            &ascii_seq,
            |b, seq| {
                let profiler = Ascii::<true>::encode_query(b"ABCDEFGHIJKLMNOPQRSTUVWXYZ").0;
                let mut result = profiler.alloc_out();
                b.iter(|| {
                    for chunk in seq.array_chunks() {
                        profiler.encode_ref(chunk, &mut result);
                        black_box(&mut result);
                    }
                })
            },
        );

        let query = b"AHOVHSHJDHFAPVHAJDJ";
        group.bench_with_input(
            BenchmarkId::new("ascii_search_20", size),
            &ascii_seq,
            |b, seq| {
                let mut deltas = vec![];
                b.iter(|| {
                    sassy::search_positions::<Ascii>(black_box(query), seq, &mut deltas);
                    black_box(&deltas);
                })
            },
        );

        let query = b"ACTGCAACTGCAACGACGTA";
        group.bench_with_input(
            BenchmarkId::new("dna_search_20", size),
            &dna_seq,
            |b, seq| {
                let mut deltas = vec![];
                b.iter(|| {
                    sassy::search_positions::<Dna>(black_box(query), seq, &mut deltas);
                    black_box(&deltas);
                })
            },
        );

        let query = b"ACTGCAACTGCAACGACGTA";
        group.bench_with_input(
            BenchmarkId::new("dna_search_20_k3", size),
            &dna_seq,
            |b, seq| {
                let mut deltas = vec![];
                b.iter(|| {
                    sassy::search_positions_bounded::<Dna>(black_box(query), seq, 3, &mut deltas);
                    black_box(&deltas);
                })
            },
        );

        let query = b"ACTGCAACTGCAACGACGTA";
        group.bench_with_input(
            BenchmarkId::new("dna_search_20_k1", size),
            &dna_seq,
            |b, seq| {
                let mut deltas = vec![];
                b.iter(|| {
                    sassy::search_positions_bounded::<Dna>(black_box(query), seq, 1, &mut deltas);
                    black_box(&deltas);
                })
            },
        );

        let query = b"ACTGCAACTGCAACGACGTA";
        group.bench_with_input(
            BenchmarkId::new("iupac_search_20", size),
            &dna_seq,
            |b, seq| {
                let mut deltas = vec![];
                b.iter(|| {
                    sassy::search_positions::<Iupac>(black_box(query), seq, &mut deltas);
                    black_box(&deltas);
                })
            },
        );

        let query = b"ACTGCAANTGCAACGACGTA";
        group.bench_with_input(
            BenchmarkId::new("iupac_search_20_N", size),
            &dna_seq,
            |b, seq| {
                let mut deltas = vec![];
                b.iter(|| {
                    sassy::search_positions::<Iupac>(black_box(query), seq, &mut deltas);
                    black_box(&deltas);
                })
            },
        );

        let query = b"ACTGCAACTGCAACGACGTAACACCTACTAAC";
        group.bench_with_input(
            BenchmarkId::new("iupac_search_32", size),
            &dna_seq,
            |b, seq| {
                let mut deltas = vec![];
                b.iter(|| {
                    sassy::search_positions::<Iupac>(black_box(query), seq, &mut deltas);
                    black_box(&deltas);
                })
            },
        );

        let query = b"ACTGCAANTGCAACGAYGTAACARCTACTAAC";
        group.bench_with_input(
            BenchmarkId::new("iupac_search_32_NRY", size),
            &dna_seq,
            |b, seq| {
                let mut deltas = vec![];
                b.iter(|| {
                    sassy::search_positions::<Iupac>(black_box(query), seq, &mut deltas);
                    black_box(&deltas);
                })
            },
        );

        let query = b"ACTGCAACTGCAACGACGTAACACCTACTAAC";
        group.bench_with_input(
            BenchmarkId::new("iupac_find_32", size),
            &dna_seq,
            |b, seq| {
                let mut deltas = vec![];
                let mut positions = vec![];
                b.iter(|| {
                    sassy::search_positions::<Iupac>(black_box(query), seq, &mut deltas);
                    positions.clear();
                    let mut costs = vec![];
                    sassy::find_below_threshold(
                        black_box(query),
                        8,
                        &deltas,
                        &mut positions,
                        &mut costs,
                    );
                    black_box(&positions);
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("dna_input_validation_valid", size),
            &dna_seq,
            |b, seq| {
                b.iter(|| {
                    let profiler = Dna::encode_query(b"").0;
                    let is_valid = profiler.valid_seq(seq);
                    black_box(is_valid);
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("iupac_input_validation_valid", size),
            &dna_seq,
            |b, seq| {
                b.iter(|| {
                    let profiler = Iupac::encode_query(b"").0;
                    let is_valid = profiler.valid_seq(seq);
                    black_box(is_valid);
                })
            },
        );
    }
    group.finish();
}

criterion_group!(benches, benchmark_base_lookup);
criterion_main!(benches);
