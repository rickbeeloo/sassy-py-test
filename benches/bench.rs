#![feature(portable_simd)]
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
    for size in [1024 * 1024] {
        let seq = generate_test_data(size);

        group.throughput(Throughput::Bytes(size as u64));

        let query_bases = b"ACGTN";

        // SIMD
        group.bench_with_input(BenchmarkId::new("SIMD", size), &seq, |b, seq| {
            b.iter(|| {
                let mut a_count = 0u32;
                let mut t_count = 0u32;
                let mut g_count = 0u32;
                let mut c_count = 0u32;
                let mut n_count = 0u32;
                let mut result = vec![];
                for chunk in seq.chunks(32) {
                    let chunk: [u8; 32] = chunk.try_into().unwrap();
                    unsafe { match_bases(&chunk, query_bases, &mut result) };
                    a_count += result[0].count_ones();
                    t_count += result[1].count_ones();
                    g_count += result[2].count_ones();
                    c_count += result[3].count_ones();
                    n_count += result[4].count_ones();
                }
                assert_eq!(a_count + t_count + g_count + c_count, seq.len() as u32);
                assert_eq!(n_count, seq.len() as u32);

                black_box((a_count, t_count, g_count, c_count, n_count))
            })
        });

        // Scalar, just to compare to counting nts
        group.bench_with_input(BenchmarkId::new("Scalar", size), &seq, |b, seq| {
            b.iter(|| {
                let mut a_count = 0u32;
                let mut t_count = 0u32;
                let mut g_count = 0u32;
                let mut c_count = 0u32;
                let mut n_count = 0u32;

                for &nt in seq {
                    match nt {
                        b'A' | b'a' => a_count += 1,
                        b'C' | b'c' => c_count += 1,
                        b'G' | b'g' => g_count += 1,
                        b'T' | b't' => t_count += 1,
                        b'N' | b'n' => n_count += 1,
                        _ => {}
                    }
                }

                black_box((a_count, t_count, g_count, c_count, n_count))
            })
        });
    }

    group.finish();
}

criterion_group!(benches, benchmark_base_lookup);
criterion_main!(benches);
