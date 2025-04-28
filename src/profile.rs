use crate::packed_nibbles_portable_64;

/// Builds a 'profile' of `b` in `64`-bit blocks, and compressed `a` into a `[0,1,2,3]` alphabet.
///
/// Returns a bitpacked `B` indicating which chars of `b` equal a given char of `a`.
///
/// TODO: 32bit variant.
pub trait Profile: Clone + std::fmt::Debug {
    /// Encoding for a single character in `a`.
    type A;
    /// Encoding for 64 characters in `b`.
    type B: Default;

    fn encode_query(a: &[u8]) -> (Self, Vec<Self::A>);
    fn encode_ref(&self, b: &[u8; 64], out: &mut Self::B);
    /// Given the encoding of an `a` and the encoding for 64 `b`s,
    /// return a bitmask of which characters of `b` equal the corresponding character of `a`.
    fn eq(ca: &Self::A, cb: &Self::B) -> u64;
}

#[derive(Clone, Debug)]
pub struct Iupac {
    bases: Vec<u8>,
}

impl Profile for Iupac {
    type A = usize;
    type B = Vec<u64>;

    fn encode_query(a: &[u8]) -> (Self, Vec<Self::A>) {
        let mut bases = vec![b'A', b'C', b'T', b'G'];
        let mut query_profile = Vec::with_capacity(a.len());
        for &c in a {
            if !bases.contains(&c) {
                bases.push(c);
            }
            query_profile.push(bases.iter().position(|&x| x == c).unwrap());
        }
        (Iupac { bases }, query_profile)
    }

    #[inline(always)]
    fn encode_ref(&self, b: &[u8; 64], out: &mut Self::B) {
        out.resize(self.bases.len(), 0);
        packed_nibbles_portable_64(b, &self.bases[4..], out);
    }

    #[inline(always)]
    fn eq(ca: &usize, cb: &Vec<u64>) -> u64 {
        unsafe { *cb.get_unchecked(*ca) }
    }
}
