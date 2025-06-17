use std::ops::{Index, IndexMut};

pub trait Profile: Clone + std::fmt::Debug {
    /// Encoding for a single character in `a`.
    type A;
    /// Encoding for 64 characters in `b`.
    type B: Index<usize, Output = u64> + IndexMut<usize, Output = u64> + Copy;
    fn encode_query(a: &[u8]) -> (Self, Vec<Self::A>);
    fn encode_ref(&self, b: &[u8; 64], out: &mut Self::B);
    /// Given the encoding of an `a` and the encoding for 64 `b`s,
    /// return a bitmask of which characters of `b` equal the corresponding character of `a`.
    fn eq(ca: &Self::A, cb: &Self::B) -> u64;
    /// Allocate a buffer of at most n_bases in search (and reuse)
    fn alloc_out() -> Self::B;
    fn n_bases(&self) -> usize;
    /// Verify whether a seqeunce matching the profile characters
    fn valid_seq(seq: &[u8]) -> bool;
    /// Return true if the two characters are a match accroding to profile
    fn is_match(char1: u8, char2: u8) -> bool;
    /// Reverse-complement the input string.
    fn reverse_complement(_query: &[u8]) -> Vec<u8> {
        unimplemented!(
            "Profile::reverse_complement not implemented for {:?}",
            std::any::type_name::<Self>()
        );
    }
    fn complement(_query: &[u8]) -> Vec<u8> {
        unimplemented!(
            "Profile::reverse_complement not implemented for {:?}",
            std::any::type_name::<Self>()
        );
    }
    fn supports_overhang() -> bool {
        unimplemented!("Profile does not support overhang");
    }
}
