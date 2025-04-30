use crate::bitpacking::compute_block_simd;
use crate::delta_encoding::V;
use crate::profiles::Dna;
use crate::profiles::Profile;
use std::simd::Simd;

#[derive(Debug)]
struct BlockState {
    vp: Simd<u64, 4>,
    vm: Simd<u64, 4>,
    hp: Vec<Simd<u64, 4>>,
    hm: Vec<Simd<u64, 4>>,
}

fn compute_traceback(query: &[u8], text: &[u8]) -> Vec<BlockState> {
    type Base = u64;
    type S = Simd<Base, 4>;

    let (profiler, query_profile) = Dna::encode_query(query);
    let mut block_states = Vec::new();

    // Process each block of 64 characters
    for (block_idx, block) in text.chunks(64).enumerate() {
        let mut vp = S::splat(0);
        let mut vm = S::splat(0);
        let mut hp = vec![S::splat(1); query.len()];
        let mut hm = vec![S::splat(0); query.len()];

        let mut profile_out = profiler.alloc_out();
        let mut slice: [u8; 64] = [b'X'; 64];
        slice[..block.len()].copy_from_slice(block);

        for (j, q) in query.iter().enumerate() {
            profiler.encode_ref(&slice, &mut profile_out);
            let eq = Dna::eq(&query_profile[j], &profile_out);
            let s = Simd::from_array([eq, eq, eq, eq]);
            compute_block_simd(&mut hp[j], &mut hm[j], &mut vp, &mut vm, s);
        }

        // Store the state for this block
        block_states.push(BlockState { vp, vm, hp, hm });
    }

    block_states
}

fn get_trace(
    lane: usize,
    position: usize,
    k: usize,
    block_states: &[BlockState],
    query_len: usize,
) -> Vec<(usize, usize)> {
    let mut trace = Vec::new();
    let mut block_idx = position / 64;
    let mut curr_j = position % 64;

    if block_idx >= block_states.len() {
        return trace;
    }

    let mut curr_i = query_len;
    let mut edits = k;

    while curr_i > 0 && (block_idx > 0 || curr_j > 0) {
        let state = &block_states[block_idx];
        let vp_bits = state.vp[lane];
        let vm_bits = state.vm[lane];

        // If we've reached the start of a block, move to previous block
        if curr_j == 0 {
            block_idx -= 1;
            curr_j = 64;
            continue;
        }

        let hp_bit = (state.hp[curr_i - 1][lane] >> (curr_j - 1)) & 1;
        let hm_bit = (state.hm[curr_i - 1][lane] >> (curr_j - 1)) & 1;
        let vp_bit = (vp_bits >> (curr_j - 1)) & 1;
        let vm_bit = (vm_bits >> (curr_j - 1)) & 1;

        let global_j = block_idx * 64 + curr_j;

        if hp_bit == 1 && hm_bit == 0 {
            trace.push((curr_i, global_j - 1));
            curr_j -= 1;
            edits += 1;
        } else if vp_bit == 1 && vm_bit == 0 {
            trace.push((curr_i - 1, global_j));
            curr_i -= 1;
            edits += 1;
        } else if hp_bit == 0 && hm_bit == 1 {
            trace.push((curr_i, global_j - 1));
            curr_j -= 1;
            edits += 1;
        } else if vp_bit == 0 && vm_bit == 1 {
            trace.push((curr_i - 1, global_j));
            curr_i -= 1;
            curr_j -= 1;
            edits += 1;
        } else {
            trace.push((curr_i - 1, global_j - 1));
            curr_i -= 1;
            curr_j -= 1;
        }
    }

    if edits == k {
        trace.reverse();
    }

    trace
}

#[test]
fn test_traceback() {
    let query = b"ACGTGGA";
    let text = b"TTTTACGTGGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACGTGGATTTTTTT";
    let block_states = compute_traceback(query, text);
    let trace = get_trace(0, 10, 0, &block_states, query.len());
    println!("Trace 1: {:?}", trace);
    //  [(6, 10), (5, 9), (4, 8), (3, 7), (2, 5), (1, 5), (0, 3)]
    // should be [(6, 10), (5, 9), (4, 8), (3, 7), (2, 6), (1, 5), (0, 4)]
    let trace = get_trace(0, text.len() - 14, 0, &block_states, query.len());
    println!("Trace 2: {:?}", trace); // FIXME: This is wrong
}

#[test]
fn test_and_block_boundary() {
    let query = b"ACGTGGA";
    let mut text = [b'G'; 128];
    text[64 - 3..64 + 4].copy_from_slice(query);
    let block_states = compute_traceback(query, &text);
    let trace = get_trace(0, 64 + 3, 0, &block_states, query.len());
    println!("Trace 1: {:?}", trace); // FIXME: This is wrong when crossing block boundary
}
