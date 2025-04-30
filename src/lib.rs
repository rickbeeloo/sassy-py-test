#![feature(portable_simd, int_roundings)]
mod bitpacking;
mod delta_encoding;
pub mod profiles {

    mod ascii;
    mod dna;
    mod iupac;
    mod profile;

    pub use profile::Profile;

    pub use ascii::{Ascii, CaseInsensitiveAscii, CaseSensitiveAscii};
    pub use dna::Dna;
    pub use iupac::Iupac;
}

mod minima;
mod search;

pub use minima::find_below_threshold;
pub use search::{search, search_bounded};
