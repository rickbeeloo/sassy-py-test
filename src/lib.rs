#![feature(portable_simd, int_roundings)]
mod bitpacking;
mod delta_encoding;
pub mod profiles;
mod search;

pub use profiles::iupac::Iupac;
pub use search::find_below_threshold;
pub use search::search;
