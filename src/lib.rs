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
pub fn match_bases_packed_nibbles(seq: &[u8; 32], query_bases: &[u8], out: &mut [u32]) {
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

        for (i, &c) in query_bases.iter().enumerate() {
            let m = _mm256_set1_epi8(get_encoded(c) as i8);
            let nz = _mm256_cmpgt_epi8(_mm256_and_si256(nib, m), zero);
            *out.get_unchecked_mut(i) = _mm256_movemask_epi8(nz) as u32;
        }
    }
}

/// NOTE: `out` must have the same length as `query_bases`.
#[inline(always)]
pub fn match_bases_2_table(seq: &[u8; 32], query_bases: &[u8], out: &mut [u32]) {
    unsafe {
        // constants & tables
        let zero = _mm256_setzero_si256();
        let mask5 = _mm256_set1_epi8(0x1F);
        let mask4 = _mm256_set1_epi8(0x0F);
        let thr16 = _mm256_set1_epi8(16);
        let tbl256 = _mm256_loadu_si256(IUPAC_CODE.as_ptr() as *const __m256i);
        let lo_tbl = _mm256_permute2x128_si256(tbl256, tbl256, 0x00);
        let hi_tbl = _mm256_permute2x128_si256(tbl256, tbl256, 0x11);

        let ptr = seq.as_ptr() as *const __m256i;
        let chunk = _mm256_loadu_si256(ptr);

        let idx5 = _mm256_and_si256(chunk, mask5);
        let is_lo = _mm256_cmpgt_epi8(thr16, idx5);
        let low4 = _mm256_and_si256(idx5, mask4);

        for (i, &c) in query_bases.iter().enumerate() {
            let m = _mm256_set1_epi8(get_encoded(c) as i8);
            let lo = _mm256_and_si256(lo_tbl, m);
            let hi = _mm256_and_si256(hi_tbl, m);

            let lo_sh = _mm256_shuffle_epi8(lo, low4);
            let hi_sh = _mm256_shuffle_epi8(hi, low4);
            let v = _mm256_blendv_epi8(hi_sh, lo_sh, is_lo);
            let nz = _mm256_cmpgt_epi8(v, zero);
            *out.get_unchecked_mut(i) = _mm256_movemask_epi8(nz) as u32;
        }
    }
}

pub fn get_match_positions(result: &[u32]) -> Vec<Vec<usize>> {
    let mut positions = Vec::with_capacity(result.len());

    for &base_result in result {
        let mut base_positions = Vec::new();

        if base_result == 0 {
            continue; // Skip if no matches
        }

        // Find set bits
        for bit_pos in 0..32 {
            if (base_result & (1u32 << bit_pos)) != 0 {
                base_positions.push(bit_pos);
            }
        }

        positions.push(base_positions);
    }

    positions
}

// Print the match positions in a readable format
pub fn print_match_positions(result: &[u32], bases: &[u8]) {
    let positions = get_match_positions(result);

    println!("Match positions for each base:");
    for (i, pos_list) in positions.iter().enumerate() {
        let base = if i < bases.len() {
            char::from(bases[i])
        } else {
            '-'
        };

        if pos_list.is_empty() {
            println!("  Base '{}': No matches", base);
        } else {
            // Just  print first 10
            let first_ten = pos_list.iter().take(10).collect::<Vec<_>>();
            println!(
                "  Base '{}': {} matches at positions {:?}",
                base,
                pos_list.len(),
                first_ten
            );
        }
    }
}
