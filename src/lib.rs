#![feature(portable_simd)]
#![feature(register_tool)]
#![register_tool(cbindgen)]

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
pub mod rec_iter;
pub mod search;
mod trace;

// Python bindings module
#[cfg(feature = "python")]
mod python;

#[cfg(feature = "c")]
#[allow(non_snake_case)]
mod c;

#[doc(hidden)]
pub mod private {
    pub use crate::minima::prefix_min;
}

#[cfg(feature = "avx512")]
const LANES: usize = 8;
#[cfg(not(feature = "avx512"))]
const LANES: usize = 4;
type S = Simd<u64, LANES>;
