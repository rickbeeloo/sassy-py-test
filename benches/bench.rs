#![feature(portable_simd)]
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use sassy::*;
use std::time::Duration;

fn generate_test_data(size: usize) -> Vec<u8> {
    let mut rng = rand::thread_rng();
    let bases = vec![b'A', b'T', b'G', b'C'];
    let mut seq = vec![b'A'; size];
    for i in 0..size {
        seq[i] = bases[rand::Rng::gen_range(&mut rng, 0..4)];
    }
    seq
}

fn benchmark_base_lookup(c: &mut Criterion) {
    let mut group = c.benchmark_group("Base lookup");
    group.measurement_time(Duration::from_secs(10));

    // Test different sequence sizes, for now just one
    for size in [56 * 10000].iter() {
        let seq = generate_test_data(*size);

        group.throughput(Throughput::Bytes(*size as u64));

        // SIMD
        group.bench_with_input(BenchmarkId::new("SIMD", size), &seq, |b, seq| {
            b.iter(|| {
                let mut a_count = 0u32;
                let mut t_count = 0u32;
                let mut g_count = 0u32;
                let mut c_count = 0u32;
                let mut n_count = 0u32;
                for chunk in seq.chunks(256) {
                    if chunk.len() == 256 {
                        let result = unsafe { match_bases::<5>(chunk, b"atgcn") };
                        a_count += result[0]
                            .to_array()
                            .iter()
                            .map(|&m| m.count_ones())
                            .sum::<u32>();
                        t_count += result[1]
                            .to_array()
                            .iter()
                            .map(|&m| m.count_ones())
                            .sum::<u32>();
                        g_count += result[2]
                            .to_array()
                            .iter()
                            .map(|&m| m.count_ones())
                            .sum::<u32>();
                        c_count += result[3]
                            .to_array()
                            .iter()
                            .map(|&m| m.count_ones())
                            .sum::<u32>();
                        n_count += result[4]
                            .to_array()
                            .iter()
                            .map(|&m| m.count_ones())
                            .sum::<u32>();
                    }
                }

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
