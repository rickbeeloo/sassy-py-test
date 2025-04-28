use crate::profiles::Profile;
use std::simd::{cmp::SimdPartialEq, u8x32};

#[derive(Clone, Debug)]
pub struct Dna {
    bases: Vec<u8>,
}

impl Profile for Dna {
    type A = usize;
    type B = [u64; 4];

    fn encode_query(a: &[u8]) -> (Self, Vec<Self::A>) {
        let bases = vec![b'A', b'C', b'T', b'G'];
        let mut query_profile = Vec::with_capacity(a.len());
        for &c in a {
            query_profile.push(bases.iter().position(|&x| x == c).unwrap());
        }
        (Dna { bases }, query_profile)
    }

    #[inline(always)]
    fn encode_ref(&self, b: &[u8; 64], out: &mut Self::B) {
        unsafe {
            let chunk0 = u8x32::from_array(b[0..32].try_into().unwrap());
            let chunk1 = u8x32::from_array(b[32..64].try_into().unwrap());
            let chunk0_shifted = chunk0 >> 1;
            let chunk1_shifted = chunk1 >> 1;
            let masked0 = chunk0_shifted & u8x32::splat(0x03);
            let masked1 = chunk1_shifted & u8x32::splat(0x03);
            for (i, code) in CODES.iter().enumerate() {
                let eq0 = masked0.simd_eq(*code);
                let eq1 = masked1.simd_eq(*code);
                let low = eq0.to_bitmask();
                let high = eq1.to_bitmask();
                *out.get_unchecked_mut(i) = (high << 32) | low;
            }
        };
    }

    #[inline(always)]
    fn eq(ca: &usize, cb: &[u64; 4]) -> u64 {
        unsafe { *cb.get_unchecked(*ca) }
    }

    #[inline(always)] // bit ugly to not use length though, but better to have static here
    fn alloc_out(&self) -> Self::B {
        [0; 4]
    }

    #[inline(always)]
    fn n_bases(&self) -> usize {
        self.bases.len()
    }
}

// Same order as iupac
const CODES: [u8x32; 4] = [
    u8x32::splat(0u8), // A
    u8x32::splat(1u8), // C
    u8x32::splat(2u8), // T
    u8x32::splat(3u8), // G
];

#[cfg(test)]
mod test {
    use super::*;

    fn get_match_positions(out: &[u64; 4]) -> Vec<Vec<usize>> {
        let mut positions = vec![vec![]; 4];
        for (i, _) in CODES.iter().enumerate() {
            let bits = out[i];
            for j in 0..64 {
                if (bits & (1u64 << j)) != 0 {
                    positions[i].push(j);
                }
            }
        }
        positions
    }

    #[test]
    fn test_dna_u64_search() {
        let mut seq = [b'G'; 64];
        seq[0] = b'A';
        seq[1] = b'A';
        seq[63] = b'C';
        let mut out = [0u64; 4];
        {
            let seq: &[u8; 64] = &seq;
            let out: &mut [u64; 4] = &mut out;
            unsafe {
                let chunk0 = u8x32::from_array(seq[0..32].try_into().unwrap());
                let chunk1 = u8x32::from_array(seq[32..64].try_into().unwrap());
                let chunk0_shifted = chunk0 >> 1;
                let chunk1_shifted = chunk1 >> 1;
                let masked0 = chunk0_shifted & u8x32::splat(0x03);
                let masked1 = chunk1_shifted & u8x32::splat(0x03);
                for (i, code) in CODES.iter().enumerate() {
                    let eq0 = masked0.simd_eq(*code);
                    let eq1 = masked1.simd_eq(*code);
                    let low = eq0.to_bitmask();
                    let high = eq1.to_bitmask();
                    *out.get_unchecked_mut(i) = (high << 32) | low;
                }
            }
        }; // A, C, T, G
        let positions = get_match_positions(&out);
        assert_eq!(positions[0], vec![0, 1]);
        assert_eq!(positions[1], vec![63]);
        assert_eq!(positions[2], vec![]);
        assert_eq!(positions[3], (2..63).collect::<Vec<_>>());
    }

    #[test]
    fn test_dna_u64_case_insensitive() {
        let mut seq = [b'G'; 64];
        seq[0] = b'a';
        seq[1] = b'A';
        let mut out = [0u64; 4];
        {
            let seq: &[u8; 64] = &seq;
            let out: &mut [u64; 4] = &mut out;
            unsafe {
                let chunk0 = u8x32::from_array(seq[0..32].try_into().unwrap());
                let chunk1 = u8x32::from_array(seq[32..64].try_into().unwrap());
                let chunk0_shifted = chunk0 >> 1;
                let chunk1_shifted = chunk1 >> 1;
                let masked0 = chunk0_shifted & u8x32::splat(0x03);
                let masked1 = chunk1_shifted & u8x32::splat(0x03);
                for (i, code) in CODES.iter().enumerate() {
                    let eq0 = masked0.simd_eq(*code);
                    let eq1 = masked1.simd_eq(*code);
                    let low = eq0.to_bitmask();
                    let high = eq1.to_bitmask();
                    *out.get_unchecked_mut(i) = (high << 32) | low;
                }
            }
        };
        let positions = get_match_positions(&out);
        assert_eq!(positions[0], vec![0, 1]);
    }
}
