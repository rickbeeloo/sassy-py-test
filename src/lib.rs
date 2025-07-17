#![feature(portable_simd)]

// INTERNAL MODS

mod bitpacking;
mod delta_encoding;
mod minima;
mod trace;

// (PARTIALLY) PUBLIC MODS

pub mod profiles;
pub mod rec_iter;
pub mod search;

pub use search::Searcher;

// BINDINGS

#[cfg(feature = "python")]
mod python;

#[cfg(feature = "c")]
mod c;

// TYPEDEFS

use std::simd::Simd;

#[cfg(feature = "avx512")]
const LANES: usize = 8;
#[cfg(not(feature = "avx512"))]
const LANES: usize = 4;

type S = Simd<u64, LANES>;

// TESTS

/// Print info on CPU features and speed of searching.
pub fn test_cpu_features() {
    eprintln!("CPU features during compilation and runtime:");
    #[cfg(target_arch = "x86_64")]
    {
        eprintln!("Target architecture: x86_64");

        let sse = if is_x86_feature_detected!("sse") {
            "+"
        } else {
            "-"
        };
        #[cfg(target_feature = "sse")]
        eprintln!("SSE  + {sse}");
        #[cfg(not(target_feature = "sse"))]
        eprintln!("SSE  - {sse}");

        let avx2 = if is_x86_feature_detected!("avx2") {
            "+"
        } else {
            "-"
        };
        #[cfg(target_feature = "avx")]
        eprintln!("AVX2 + {avx2}");
        #[cfg(not(target_feature = "avx"))]
        eprintln!("AVX2 - {avx2}");

        let bmi2 = if is_x86_feature_detected!("bmi2") {
            "+"
        } else {
            "-"
        };
        #[cfg(target_feature = "bmi2")]
        eprintln!("BMI2 + {bmi2}");
        #[cfg(not(target_feature = "bmi2"))]
        eprintln!("BMI2 - {bmi2}");
    }
    #[cfg(target_arch = "aarch64")]
    {
        use std::arch::is_aarch64_feature_detected;

        eprintln!("Target architecture: aarch64 currently unsupported");

        let neon = if is_aarch64_feature_detected!("neon") {
            "+"
        } else {
            "-"
        };
        #[cfg(target_feature = "neon")]
        eprintln!("NEON + {neon}");
        #[cfg(not(target_feature = "neon"))]
        eprintln!("NEON - {neon}");
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        eprintln!("Unsupported target architecture");
    }
}

pub fn test_throughput() {
    eprintln!("Running a little test: aligning a 23bp pattern against 100kb text, with 1 error.");

    use rand::Rng;
    let n = 100000;
    let m = 23;
    let k = 1;

    let mut rng = rand::rng();
    let text: Vec<u8> = (0..n).map(|_| b"ACGT"[rng.random_range(0..4)]).collect();
    let pattern: Vec<u8> = (0..m).map(|_| b"ACGT"[rng.random_range(0..4)]).collect();

    let mut searcher = Searcher::<profiles::Dna>::new(false, None);
    let start = std::time::Instant::now();
    let _matches = searcher.search(&pattern, &text, k);
    let duration = start.elapsed();
    eprintln!(
        "Search in GB/s: {}",
        text.len() as f32 / duration.as_secs_f32() / 1_000_000_000.0
    );
}

#[cfg(test)]
mod test {
    #[test]
    fn test_cpu_features() {
        super::test_cpu_features();
    }
    #[test]
    fn test_throughput() {
        super::test_throughput();
    }
}
