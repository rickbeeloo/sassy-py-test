#![feature(portable_simd, int_roundings, let_chains)]

use std::simd::Simd;

mod bitpacking;
mod delta_encoding;

pub mod profiles {
    mod ascii;
    mod dna;
    mod iupac;
    mod profile;

    pub use ascii::{Ascii, CaseInsensitiveAscii, CaseSensitiveAscii};
    pub use dna::Dna;
    pub use iupac::Iupac;
    pub use profile::Profile;
}

mod minima;
pub mod search;
mod trace;

// Python bindings module
#[cfg(feature = "python")]
mod python;

#[doc(hidden)]
pub mod private {
    pub use crate::minima::prefix_min;
}

#[cfg(feature = "avx512")]
const LANES: usize = 8;
#[cfg(not(feature = "avx512"))]
const LANES: usize = 4;
type S = Simd<u64, LANES>;
