#![feature(portable_simd)]
use arrayvec::ArrayVec;
use core::simd::Simd;
use std::arch::x86_64::*;

#[rustfmt::skip]
pub const IUPAC_CODE: [u8; 32] = {
    let mut t = [0u8; 32];
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
    t[b'N' as usize & 0x1F] = A|C|T|G;
    
    // IUPAC ambiguity codes
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

pub fn get_encoded(c: u8) -> u8 {
    IUPAC_CODE[(c & 0x1F) as usize]
}

#[target_feature(enable = "avx2")]
pub unsafe fn match_bases<const M: usize>(
    seq: &[u8],
    bases: &[u8; M],
) -> ArrayVec<Simd<u64, 4>, 32> {
    assert!(seq.len() == 256, "Sequence must be exactly 256 bytes");
    assert!(M <= 32, "At most 32 bases supported");

    // constants & tables
    let zero = _mm256_setzero_si256();
    let mask5 = _mm256_set1_epi8(0x1F);
    let mask4 = _mm256_set1_epi8(0x0F);
    let thr16 = _mm256_set1_epi8(16);
    let tbl256 = _mm256_loadu_si256(IUPAC_CODE.as_ptr() as *const __m256i);
    let lo_tbl = _mm256_permute2x128_si256(tbl256, tbl256, 0x00);
    let hi_tbl = _mm256_permute2x128_si256(tbl256, tbl256, 0x11);

    // Precompute lo/hi masks
    let mut lo_masks = [_mm256_setzero_si256(); 32]; // Use fixed size 32
    let mut hi_masks = [_mm256_setzero_si256(); 32];

    for i in 0..M {
        let m = _mm256_set1_epi8(get_encoded(bases[i]) as i8);
        lo_masks[i] = _mm256_and_si256(lo_tbl, m);
        hi_masks[i] = _mm256_and_si256(hi_tbl, m);
    }

    // M × 8×u32 bitmasks
    let mut acc = [[0u32; 8]; 32]; // Use fixed size 32

    for chunk_idx in 0..8 {
        let ptr = seq.as_ptr().add(chunk_idx * 32) as *const __m256i;
        let chunk = _mm256_loadu_si256(ptr);

        let idx5 = _mm256_and_si256(chunk, mask5);
        let is_lo = _mm256_cmpgt_epi8(thr16, idx5);
        let low4 = _mm256_and_si256(idx5, mask4);

        for j in 0..M {
            let lo_sh = _mm256_shuffle_epi8(lo_masks[j], low4);
            let hi_sh = _mm256_shuffle_epi8(hi_masks[j], low4);
            let v = _mm256_blendv_epi8(hi_sh, lo_sh, is_lo);
            let nz = _mm256_cmpgt_epi8(v, zero);
            acc[j][chunk_idx] = _mm256_movemask_epi8(nz) as u32;
        }
    }

    // Pack into lanes and store in arrayvec output
    let mut out = ArrayVec::<Simd<u64, 4>, 32>::new();
    for i in 0..M {
        out.push(Simd::from_array([
            (acc[i][0] as u64) | ((acc[i][1] as u64) << 32),
            (acc[i][2] as u64) | ((acc[i][3] as u64) << 32),
            (acc[i][4] as u64) | ((acc[i][5] as u64) << 32),
            (acc[i][6] as u64) | ((acc[i][7] as u64) << 32),
        ]));
    }
    out
}

pub fn get_match_positions(result: &ArrayVec<Simd<u64, 4>, 32>) -> Vec<Vec<usize>> {
    let mut positions = Vec::with_capacity(result.len());

    for (base_idx, base_result) in result.iter().enumerate() {
        let mut base_positions = Vec::new();

        // Process each of the 4 u64 chunks
        for chunk_idx in 0..4 {
            let chunk = base_result.as_array()[chunk_idx];
            if chunk == 0 {
                continue; // Skip if no matches
            }

            let base_offset = chunk_idx * 64;
            // Find set bits
            for bit_pos in 0..64 {
                if (chunk & (1u64 << bit_pos)) != 0 {
                    base_positions.push(base_offset + bit_pos);
                }
            }
        }

        positions.push(base_positions);
    }

    positions
}

// Print the match positions in a readable format
pub fn print_match_positions(result: &ArrayVec<Simd<u64, 4>, 32>, bases: &[u8]) {
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

fn main() {
    // Normal > IUPAC
    let mut seq = vec![b'g'; 256];
    seq[0] = b'a';
    seq[1] = b'y'; // C or T
    let bases = b"ATGC";
    let result = unsafe { match_bases(&seq, bases) };
    print_match_positions(&result, bases);

    // IUPAC > normal
    let mut seq = vec![b'g'; 256];
    seq[0] = b'a';
    seq[1] = b'y'; // C or T
    seq[2] = b'C';
    let bases = b"Y";
    let result = unsafe { match_bases(&seq, bases) };
    print_match_positions(&result, bases);
}
