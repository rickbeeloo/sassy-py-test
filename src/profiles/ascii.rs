use crate::profiles::Profile;
use std::{
    mem::transmute,
    simd::{
        cmp::{SimdPartialEq, SimdPartialOrd},
        u8x32,
    },
};

#[derive(Clone, Debug)]
pub struct Ascii<const CASE_SENSITIVE: bool = true> {
    bases: Vec<u8>,
}

pub type CaseSensitiveAscii = Ascii<true>;
pub type CaseInsensitiveAscii = Ascii<false>;

impl<const CASE_SENSITIVE: bool> Profile for Ascii<CASE_SENSITIVE> {
    type A = usize;
    type B = [u64; 256]; // Maximum number of ASCII characters

    fn encode_query(a: &[u8]) -> (Self, Vec<Self::A>) {
        let mut bases = Vec::new();
        let mut query_profile = Vec::with_capacity(a.len());
        for &c in a {
            if !bases.contains(&c) {
                bases.push(c);
            }
            query_profile.push(bases.iter().position(|&x| x == c).unwrap());
        }
        (Ascii { bases }, query_profile)
    }

    #[inline(always)]
    fn encode_ref(&self, b: &[u8; 64], out: &mut Self::B) {
        if CASE_SENSITIVE {
            ascii_u64_search(b, &self.bases, out);
        } else {
            ascii_u64_search_case_insensitive(b, &self.bases, out);
        }
    }

    #[inline(always)]
    fn eq(ca: &usize, cb: &[u64; 256]) -> u64 {
        unsafe { *cb.get_unchecked(*ca) }
    }

    #[inline(always)]
    fn is_match(char1: u8, char2: u8) -> bool {
        if CASE_SENSITIVE {
            char1.eq(&char2)
        } else {
            // Safe rust version to handle cases only in Az range
            char1.eq_ignore_ascii_case(&char2)
        }
    }

    #[inline(always)]
    fn alloc_out() -> Self::B {
        [0; 256]
    }

    #[inline(always)]
    fn n_bases(&self) -> usize {
        self.bases.len()
    }

    #[inline(always)]
    fn valid_seq(_seq: &[u8]) -> bool {
        true // assuming every u8 is valid ascii
    }
}

#[inline(always)]
pub fn ascii_u64_search(seq: &[u8; 64], bases: &[u8], out: &mut [u64]) {
    unsafe {
        let chunk0 = u8x32::from_array(seq[0..32].try_into().unwrap());
        let chunk1 = u8x32::from_array(seq[32..64].try_into().unwrap());

        for (i, &base) in bases.iter().enumerate() {
            let m = u8x32::splat(base);
            let eq0 = chunk0.simd_eq(m);
            let eq1 = chunk1.simd_eq(m);
            let low = eq0.to_bitmask();
            let high = eq1.to_bitmask();
            *out.get_unchecked_mut(i) = (high << 32) | low;
        }
    }
}

// FIXME: Tests
#[inline(always)]
fn ascii_u64_search_case_insensitive(seq: &[u8; 64], bases: &[u8], out: &mut [u64]) {
    unsafe {
        let chunk0 = u8x32::from_array(seq[0..32].try_into().unwrap());
        let chunk1 = u8x32::from_array(seq[32..64].try_into().unwrap());

        const A: u8 = b'A';
        const Z: u8 = b'Z';
        let to_lowercase = b'a' - b'A';
        let is_char0 = chunk0.simd_ge(u8x32::splat(A)) & chunk0.simd_le(u8x32::splat(Z));
        let is_char1 = chunk1.simd_ge(u8x32::splat(A)) & chunk1.simd_le(u8x32::splat(Z));
        // Transmute from i8x32 to u8x32
        let lower0 =
            chunk0 | (u8x32::splat(to_lowercase) & transmute::<_, u8x32>(is_char0.to_int()));
        let lower1 =
            chunk1 | (u8x32::splat(to_lowercase) & transmute::<_, u8x32>(is_char1.to_int()));

        for (i, &base) in bases.iter().enumerate() {
            let m = u8x32::splat(base | 0x20);
            let eq0 = lower0.simd_eq(m);
            let eq1 = lower1.simd_eq(m);
            let low = eq0.to_bitmask();
            let high = eq1.to_bitmask();
            *out.get_unchecked_mut(i) = (high << 32) | low;
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn get_match_positions(out: &[u64]) -> Vec<Vec<usize>> {
        let mut positions = vec![vec![]; out.len()];
        for (i, _) in out.iter().enumerate() {
            let bits = out[i];
            for j in 0..64 {
                if (bits & (1u64 << j)) != 0 {
                    positions[i].push(j);
                }
            }
        }
        positions
    }

    const HELLO_TEST_SEQ: [u8; 64] = {
        let mut seq = [b'H'; 64];
        seq[0] = b'E';
        seq[1] = b'l';
        seq[2] = b'L';
        seq[3] = b'o';
        seq
    };

    const HELLO_TEST_BASES: [u8; 3] = [b'H', b'l', b'o'];

    #[test]
    fn test_ascii_is_match() {
        // Case sensitive
        assert!(Ascii::<true>::is_match(b'H', b'H'));
        assert!(!Ascii::<true>::is_match(b'l', b'L')); // Should not match
        // Case insensitive
        assert!(Ascii::<false>::is_match(b'H', b'H'));
        assert!(Ascii::<false>::is_match(b'l', b'L')); // Should match now
    }

    #[test]
    fn test_ascii_u64_search() {
        let mut out = vec![0u64; 3];
        ascii_u64_search(&HELLO_TEST_SEQ, &HELLO_TEST_BASES, &mut out);
        let positions = get_match_positions(&out);
        assert_eq!(positions[0], (4..64).collect::<Vec<_>>());
        assert_eq!(positions[1], vec![1]);
        assert_eq!(positions[2], vec![3]);
    }

    #[test]
    fn test_ascii_u64_search_case_insensitive() {
        let mut out = vec![0u64; 3];
        ascii_u64_search_case_insensitive(&HELLO_TEST_SEQ, &HELLO_TEST_BASES, &mut out);
        let positions = get_match_positions(&out);
        assert_eq!(positions[1], vec![1, 2]); // l and L
    }

    #[test]
    fn test_ascii_u64_search_case_sensitive() {
        let mut out = vec![0u64; 3];
        ascii_u64_search(&HELLO_TEST_SEQ, &HELLO_TEST_BASES, &mut out);
        let positions = get_match_positions(&out);
        assert_eq!(positions[1], vec![1]); // only l
    }
}
