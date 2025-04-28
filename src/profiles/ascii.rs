use crate::profiles::trai_def::Profile;
use std::simd::{cmp::SimdPartialEq, u8x32};

#[derive(Clone, Debug)]
pub struct Ascii {
    bases: Vec<u8>,
}

impl Profile for Ascii {
    type A = usize;
    type B = Vec<u64>;

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
        ascii_u64_search(b, &self.bases, out);
    }

    #[inline(always)]
    fn eq(ca: &usize, cb: &Vec<u64>) -> u64 {
        unsafe { *cb.get_unchecked(*ca) }
    }

    #[inline(always)]
    fn alloc_out(&self) -> Self::B {
        vec![0; self.bases.len()]
    }

    #[inline(always)]
    fn n_bases(&self) -> usize {
        self.bases.len()
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

#[inline(always)]
pub fn ascii_u64_search_case_insensitive(seq: &[u8; 64], bases: &[u8], out: &mut [u64]) {
    unsafe {
        let chunk0 = u8x32::from_array(seq[0..32].try_into().unwrap());
        let chunk1 = u8x32::from_array(seq[32..64].try_into().unwrap());

        let lower0 = chunk0 | u8x32::splat(0x20);
        let lower1 = chunk1 | u8x32::splat(0x20);

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
