use crate::bitpacking::compute_block_simd;
use crate::delta_encoding::V;
use crate::profiles::Dna;
use crate::profiles::Profile;
use std::simd::Simd;

fn compute_traceback(query: &[u8], text: &[u8]) {
    type Base = u64;
    type VV = V<Base>;
    type S = Simd<Base, 4>;

    let (profiler, query_profile) = Dna::encode_query(query);

    let mut vp = S::splat(0);
    let mut vm = S::splat(0);
    let mut hp = vec![S::splat(1); query.len()];
    let mut hm = vec![S::splat(0); query.len()];

    let mut profile_out = profiler.alloc_out();
    for block in text.chunks(64) {
        for (j, q) in query.iter().enumerate() {
            let mut slice: [u8; 64] = [b'X'; 64];
            slice[..block.len()].copy_from_slice(block);
            profiler.encode_ref(&slice, &mut profile_out);
            let eq = Dna::eq(&query_profile[j], &profile_out);
            let s = Simd::from_array([eq, eq, eq, eq]);
            compute_block_simd(&mut hp[j], &mut hm[j], &mut vp, &mut vm, s);
        }
    }

    println!("hp: {:?}", hp);
    println!("hm: {:?}", hm);
    println!("vp: {:?}", vp);
    println!("vm: {:?}", vm);
}

fn get_trace(
    query_len: usize,
    text_len: usize,
    vp: Simd<u64, 4>,
    vm: Simd<u64, 4>,
    hp: &[Simd<u64, 4>],
    hm: &[Simd<u64, 4>],
) -> Vec<Vec<(usize, usize)>> {
    let mut traces = Vec::new();

    for lane in 0..4 {
        let mut trace = Vec::new();
        let mut i = query_len;
        let mut j = text_len;

        let vp_bits = vp[lane];
        let vm_bits = vm[lane];

        while i > 0 && j > 0 {
            let hp_bit = (hp[i - 1][lane] >> (j - 1)) & 1;
            let hm_bit = (hm[i - 1][lane] >> (j - 1)) & 1;
            let vp_bit = (vp_bits >> (j - 1)) & 1;
            let vm_bit = (vm_bits >> (j - 1)) & 1;

            if hp_bit == 1 && hm_bit == 0 {
                // H+ del
                trace.push((i, j - 1));
                j -= 1;
            } else if vp_bit == 1 && vm_bit == 0 {
                // V+ ins
                trace.push((i - 1, j));
                i -= 1;
            } else if hp_bit == 0 && hm_bit == 1 {
                // H- del
                trace.push((i, j - 1));
                j -= 1;
            } else if vp_bit == 0 && vm_bit == 1 {
                // V- ins
                trace.push((i - 1, j));
                i -= 1;
            } else {
                // Diag
                trace.push((i - 1, j - 1));
                i -= 1;
                j -= 1;
            }
        }

        // Rest
        while i > 0 {
            trace.push((i - 1, j));
            i -= 1;
        }
        while j > 0 {
            trace.push((i, j - 1));
            j -= 1;
        }

        trace.reverse();
        traces.push(trace);
    }

    traces
}

mod test {
    use super::*;

    #[test]
    fn test_traceback() {
        let query = b"ACGTGGA";
        let text = b"TTTTACGTGGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        compute_traceback(query, text);
    }
}
