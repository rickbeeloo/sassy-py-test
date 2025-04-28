#![feature(portable_simd)]
use std::arch::x86_64::*;

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

/// NOTE: `out` must have the same length as `query_bases`.
#[inline(always)]
pub fn packed_nibbles(seq: &[u8; 32], extra_bases: &[u8], out: &mut [u32]) {
    unsafe {
        let zero = _mm256_setzero_si256();
        let mask5 = _mm256_set1_epi8(0x1F);
        let mask4 = _mm256_set1_epi8(0x0F);
        let thr15: __m256i = _mm256_set1_epi8(15);
        let tbl128 = _mm_load_si128(PACKED_NIBBLES.0.as_ptr() as *const __m128i);
        let tbl256 = _mm256_broadcastsi128_si256(tbl128);

        let ptr = seq.as_ptr() as *const __m256i;
        let chunk = _mm256_loadu_si256(ptr);

        let idx5 = _mm256_and_si256(chunk, mask5);
        let low4 = _mm256_and_si256(idx5, mask4);
        let is_hi = _mm256_cmpgt_epi8(idx5, thr15);

        let shuffled = _mm256_shuffle_epi8(tbl256, low4);
        let lo_nib = _mm256_and_si256(shuffled, mask4);
        let hi_nib = _mm256_and_si256(_mm256_srli_epi16(shuffled, 4), mask4);
        let nib = _mm256_blendv_epi8(lo_nib, hi_nib, is_hi);

        let a_mask = _mm256_set1_epi8(get_encoded(b'A') as i8);
        let t_mask = _mm256_set1_epi8(get_encoded(b'T') as i8);
        let g_mask = _mm256_set1_epi8(get_encoded(b'G') as i8);
        let c_mask = _mm256_set1_epi8(get_encoded(b'C') as i8);

        let a_match = _mm256_cmpgt_epi8(_mm256_and_si256(nib, a_mask), zero);
        let t_match = _mm256_cmpgt_epi8(_mm256_and_si256(nib, t_mask), zero);
        let g_match = _mm256_cmpgt_epi8(_mm256_and_si256(nib, g_mask), zero);
        let c_match = _mm256_cmpgt_epi8(_mm256_and_si256(nib, c_mask), zero);

        *out.get_unchecked_mut(0) = _mm256_movemask_epi8(a_match) as u32;
        *out.get_unchecked_mut(1) = _mm256_movemask_epi8(t_match) as u32;
        *out.get_unchecked_mut(2) = _mm256_movemask_epi8(g_match) as u32;
        *out.get_unchecked_mut(3) = _mm256_movemask_epi8(c_match) as u32;

        for (i, &c) in extra_bases.iter().enumerate() {
            let m = _mm256_set1_epi8(get_encoded(c) as i8);
            let nz = _mm256_cmpgt_epi8(_mm256_and_si256(nib, m), zero);
            *out.get_unchecked_mut(i + 4) = _mm256_movemask_epi8(nz) as u32;
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn get_match_positions(result: &[u32]) -> Vec<Vec<usize>> {
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

    #[test]
    fn test_just_atgc() {
        let mut seq = [b'g'; 32];
        seq[0] = b'a';
        seq[1] = b'y'; // C or T
        let mut result = [0u32; 32];
        packed_nibbles(&seq, b"", &mut result);
        let positions = get_match_positions(&result);
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
        packed_nibbles(&seq, b"NY", &mut result);
        let positions = get_match_positions(&result);
        let n_positions = positions[4].clone();
        let y_positions = positions[5].clone();
        // N matches all positions
        assert_eq!(n_positions, (0..32).collect::<Vec<_>>());
        // Y matches 1,2
        assert_eq!(y_positions, vec![1, 2]);
    }
}
