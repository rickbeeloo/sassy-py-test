#![feature(portable_simd, array_chunks)]
use criterion::{BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main};
use sassy::profiles::{ascii::*, dna::*, iupac::*};
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

        let query_bases_defaults = b"NY";
        let query_ascii_bases = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ";

        group.bench_with_input(BenchmarkId::new("iupac_u32", size), &dna_seq, |b, seq| {
            b.iter(|| {
                let mut result = vec![0; 6];
                for chunk in seq.array_chunks() {
                    iupac_u32_search(chunk, black_box(query_bases_defaults), &mut result);
                    black_box(&mut result);
                }
            })
        });

        group.bench_with_input(BenchmarkId::new("iupac_u64", size), &dna_seq, |b, seq| {
            b.iter(|| {
                let mut result = vec![0; 6];
                for chunk in seq.array_chunks() {
                    iupac_u64_search(chunk, black_box(query_bases_defaults), &mut result);
                    black_box(&mut result);
                }
            })
        });

        group.bench_with_input(BenchmarkId::new("dna_u64", size), &dna_seq, |b, seq| {
            b.iter(|| {
                let mut result = [0u64; 4];
                for chunk in seq.array_chunks() {
                    dna_u64_search(chunk, &mut result);
                    black_box(&mut result);
                }
            })
        });

        group.bench_with_input(BenchmarkId::new("ascii_u64", size), &ascii_seq, |b, seq| {
            b.iter(|| {
                let mut result = vec![0; query_ascii_bases.len()];
                for chunk in seq.array_chunks() {
                    ascii_u64_search(chunk, black_box(query_ascii_bases), &mut result);
                    black_box(&mut result);
                }
            })
        });

        let query = b"AHOVHSHJDHFAPVHAJDJ";
        group.bench_with_input(
            BenchmarkId::new("ascii_search_20", size),
            &ascii_seq,
            |b, seq| {
                let mut deltas = vec![];
                b.iter(|| {
                    sassy::search::<Ascii>(black_box(query), seq, &mut deltas);
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
                    sassy::search::<Dna>(black_box(query), seq, &mut deltas);
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
                    sassy::search::<Iupac>(black_box(query), seq, &mut deltas);
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
                    sassy::search::<Iupac>(black_box(query), seq, &mut deltas);
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
                    sassy::search::<Iupac>(black_box(query), seq, &mut deltas);
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
                    sassy::search::<Iupac>(black_box(query), seq, &mut deltas);
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
                    sassy::search::<Iupac>(black_box(query), seq, &mut deltas);
                    positions.clear();
                    sassy::find_below_threshold(black_box(query), 8, &deltas, &mut positions);
                    black_box(&positions);
                })
            },
        );
    }

    group.finish();
}

criterion_group!(benches, benchmark_base_lookup);
criterion_main!(benches);
