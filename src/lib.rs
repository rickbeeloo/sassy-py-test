#![feature(portable_simd, int_roundings)]
mod bitpacking;
mod delta_encoding;
mod profile;
mod search;

use std::{
    arch::x86_64::*,
    mem::transmute,
    simd::{cmp::SimdPartialOrd, u8x32},
};

#[rustfmt::skip]
pub const IUPAC_CODE: [u8; 32] = {
    let mut t = [255u8; 32];
    // Standard bases
    // Map ACGT -> [0,1,3,2], like packed_seq does.
    const A: u8 = 1<<0;
    const C: u8 = 1<<1;
    const G: u8 = 1<<3;
    const T: u8 = 1<<2;

    // Map common chars.
    // Lower case has the same last 5 bits as upper case.
    // (Thanks ASCII :)
    t[b'A' as usize & 0x1F] = A;
    t[b'C' as usize & 0x1F] = C;
    t[b'G' as usize & 0x1F] = G;
    t[b'T' as usize & 0x1F] = T;
    t[b'U' as usize & 0x1F] = T;
    t[b'N' as usize & 0x1F] = A|C|T|G;
    
    // IUPAC ambiguity codes
    // https://www.bioinformatics.org/sms/iupac.html
    t[b'R' as usize & 0x1F] = A|G;
    t[b'Y' as usize & 0x1F] = C|T;
    t[b'S' as usize & 0x1F] = G|C;
    t[b'W' as usize & 0x1F] = A|T;
    t[b'K' as usize & 0x1F] = G|T;
    t[b'M' as usize & 0x1F] = A|C;
    t[b'B' as usize & 0x1F] = C|G|T;
    t[b'D' as usize & 0x1F] = A|G|T;
    t[b'H' as usize & 0x1F] = A|C|T;
    t[b'V' as usize & 0x1F] = A|C|G;
    
    // Gap/unknown
    t[b'X' as usize & 0x1F] = 0;
    
    t
};

#[inline(always)]
pub fn get_encoded(c: u8) -> u8 {
    IUPAC_CODE[(c & 0x1F) as usize]
}

#[repr(align(16))]
pub struct AlignedPacked([u8; 16]);

const PACKED_NIBBLES: AlignedPacked = AlignedPacked({
    let mut p = [0u8; 16];
    let mut i = 0;
    while i < 16 {
        let lo = IUPAC_CODE[i] & 0x0F;
        let hi = IUPAC_CODE[i + 16] & 0x0F;
        // packed 8 bit of low nibbles(0-3) and high nibbles(4-7)
        p[i] = (hi << 4) | lo;
        i += 1;
    }
    p
});

/// NOTE: `out` must have the same length as `4 + extra_bases`.
#[inline(always)]
pub fn packed_nibbles_portable_32(seq: &[u8; 32], extra_bases: &[u8], out: &mut [u32]) {
    unsafe {
        let zero = u8x32::splat(0);
        let mask4 = u8x32::splat(0x0F);
        let tbl256 = u8x32::from_array(transmute([PACKED_NIBBLES.0, PACKED_NIBBLES.0]));

        let chunk = u8x32::from_array(*seq);

        let idx5 = chunk & u8x32::splat(0x1F);
        let low4 = idx5 & mask4;
        let is_hi = idx5.simd_gt(u8x32::splat(15));

        let shuffled: u8x32 = transmute(_mm256_shuffle_epi8(transmute(tbl256), transmute(low4)));
        // This &mask4 is redundant because we only ever look at the low 4 bits anyway,
        // but removing it worsens codegen.
        let lo_nib = shuffled & mask4;
        let hi_nib = shuffled >> 4;
        let nib = is_hi.select(hi_nib, lo_nib);

        for (i, base) in [b'A', b'T', b'G', b'C'].iter().enumerate() {
            let m = u8x32::splat(get_encoded(*base));
            let nz = (nib & m).simd_gt(zero);
            *out.get_unchecked_mut(i) = nz.to_bitmask() as u32;
        }

        for (i, &c) in extra_bases.iter().enumerate() {
            let m = u8x32::splat(get_encoded(c));
            let nz = (nib & m).simd_gt(zero);
            *out.get_unchecked_mut(i + 4) = nz.to_bitmask() as u32;
        }
    }
}

/// NOTE: `out` must have the same length as `4 + extra_bases`.
#[inline(always)]
pub fn packed_nibbles_portable_64(seq: &[u8; 64], extra_bases: &[u8], out: &mut [u64]) {
    unsafe {
        let zero = u8x32::splat(0);
        let mask4 = u8x32::splat(0x0F);
        let tbl256 = u8x32::from_array(transmute([PACKED_NIBBLES.0, PACKED_NIBBLES.0]));

        let chunk0 = u8x32::from_array(seq[0..32].try_into().unwrap());
        let chunk1 = u8x32::from_array(seq[32..64].try_into().unwrap());

        let idx5_0 = chunk0 & u8x32::splat(0x1F);
        let idx5_1 = chunk1 & u8x32::splat(0x1F);
        let low4_0 = idx5_0 & mask4;
        let low4_1 = idx5_1 & mask4;

        let is_hi_0 = idx5_0.simd_ge(u8x32::splat(15));
        let is_hi_1 = idx5_1.simd_ge(u8x32::splat(15));

        let shuffled0: u8x32 = transmute(_mm256_shuffle_epi8(transmute(tbl256), transmute(low4_0)));
        let shuffled1: u8x32 = transmute(_mm256_shuffle_epi8(transmute(tbl256), transmute(low4_1)));

        let lo_nib0 = shuffled0 & mask4;
        let lo_nib1 = shuffled1 & mask4;

        let hi_nib0 = shuffled0 >> 4;
        let hi_nib1 = shuffled1 >> 4;

        let nib0 = is_hi_0.select(hi_nib0, lo_nib0);
        let nib1 = is_hi_1.select(hi_nib1, lo_nib1);

        for (i, &base) in [b'A', b'T', b'G', b'C'].iter().enumerate() {
            let m = u8x32::splat(get_encoded(base));

            let match0 = (nib0 & m).simd_gt(zero);
            let match1 = (nib1 & m).simd_gt(zero);

            let low = match0.to_bitmask() as u64;
            let high = match1.to_bitmask() as u64;

            *out.get_unchecked_mut(i) = (high << 32) | low;
        }

        for (i, &base) in extra_bases.iter().enumerate() {
            let m = u8x32::splat(get_encoded(base));

            let match0 = (nib0 & m).simd_gt(zero);
            let match1 = (nib1 & m).simd_gt(zero);

            let low = match0.to_bitmask() as u64;
            let high = match1.to_bitmask() as u64;

            *out.get_unchecked_mut(i + 4) = (high << 32) | low;
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn get_match_positions_u32(result: &[u32]) -> Vec<Vec<usize>> {
        result
            .iter()
            .filter_map(|&base_result| {
                if base_result == 0 {
                    None
                } else {
                    let positions: Vec<usize> = (0..32)
                        .filter(|&pos| (base_result & (1 << pos)) != 0)
                        .collect();
                    Some(positions)
                }
            })
            .collect()
    }

    fn get_match_positions_u64(result: &[u64]) -> Vec<Vec<usize>> {
        result
            .iter()
            .filter_map(|&base_result| {
                if base_result == 0 {
                    None
                } else {
                    let positions: Vec<usize> = (0..64)
                        .filter(|&pos| (base_result & (1 << pos)) != 0)
                        .collect();
                    Some(positions)
                }
            })
            .collect()
    }

    #[test]
    fn test_just_atgc() {
        let mut seq = [b'g'; 32];
        seq[0] = b'a';
        seq[1] = b'y'; // C or T
        let mut result = [0u32; 32];
        packed_nibbles_portable_32(&seq, b"", &mut result);
        let positions = get_match_positions_u32(&result);
        let a_positions = positions[0].clone();
        let t_positions = positions[1].clone();
        let g_positions = positions[2].clone();
        let c_positions = positions[3].clone();
        assert_eq!(a_positions, vec![0]);
        assert_eq!(t_positions, vec![1]);
        assert_eq!(g_positions, (2..32).collect::<Vec<_>>());
        assert_eq!(c_positions, vec![1]);
    }

    #[test]
    fn test_extra_bases_ny() {
        let mut seq = [b'g'; 32];
        seq[0] = b'a'; // Does not match Y
        seq[1] = b'y'; // Matches Y
        seq[2] = b'C'; // Matches Y
        let mut result = [0u32; 32];
        packed_nibbles_portable_32(&seq, b"NY", &mut result);
        let positions = get_match_positions_u32(&result);
        let n_positions = positions[4].clone();
        let y_positions = positions[5].clone();
        // N matches all positions
        assert_eq!(n_positions, (0..32).collect::<Vec<_>>());
        // Y matches 1,2
        assert_eq!(y_positions, vec![1, 2]);
    }

    #[test]
    fn test_just_atgc_64() {
        let mut seq = [b'g'; 64];
        seq[0] = b'a';
        seq[1] = b'y'; // C or T
        seq[34] = b'y'; // C or T
        let mut result = [0u64; 64];
        packed_nibbles_portable_64(&seq, b"", &mut result);
        let positions = get_match_positions_u64(&result);
        let a_positions = positions[0].clone();
        let t_positions = positions[1].clone();
        let g_positions = positions[2].clone();
        let c_positions = positions[3].clone();
        assert_eq!(a_positions, vec![0]);
        assert_eq!(t_positions, vec![1, 34]);
        assert_eq!(
            g_positions,
            [
                &(2..34).collect::<Vec<_>>()[..], // 34 not inclusive
                &(35..64).collect::<Vec<_>>()[..]
            ]
            .concat()
        );
        assert_eq!(c_positions, vec![1, 34]);
    }

    #[test]
    fn test_extra_bases_ny_64() {
        let mut seq = [b'g'; 64];
        seq[0] = b'a'; // Does not match Y
        seq[1] = b'y'; // Matches Y
        seq[2] = b'C'; // Matches Y
        seq[50] = b'y'; // Matches Y
        seq[63] = b'y'; // Matches Y
        let mut result = [0u64; 64];
        packed_nibbles_portable_64(&seq, b"NY", &mut result);
        let positions = get_match_positions_u64(&result);
        let n_positions = positions[4].clone();
        let y_positions = positions[5].clone();
        // N matches all positions
        assert_eq!(n_positions, (0..64).collect::<Vec<_>>());
        assert_eq!(y_positions, vec![1, 2, 50, 63]);
    }
}

pub use profile::Iupac;
pub use search::search;
