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
    let block_idx = position / 64;
    let local_pos = position % 64;

    if block_idx >= block_states.len() {
        return trace;
    }

    let state = &block_states[block_idx];
    let mut curr_i = query_len;
    let mut curr_j = local_pos;
    let mut edits = 0;

    let vp_bits = state.vp[lane];
    let vm_bits = state.vm[lane];

    while curr_i > 0 && curr_j > 0 && edits < k {
        let hp_bit = (state.hp[curr_i - 1][lane] >> (curr_j - 1)) & 1;
        let hm_bit = (state.hm[curr_i - 1][lane] >> (curr_j - 1)) & 1;
        let vp_bit = (vp_bits >> (curr_j - 1)) & 1;
        let vm_bit = (vm_bits >> (curr_j - 1)) & 1;

        let global_j = block_idx * 64 + curr_j; // Convert to global position

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
            edits += 1;
        } else {
            trace.push((curr_i - 1, global_j - 1));
            curr_i -= 1;
            curr_j -= 1;
        }

        // If we reach the start of this block but still have edits to find
        if curr_j == 0 && block_idx > 0 {
            panic!("Reached end already??")
        }
    }

    if edits == k {
        trace.reverse();
    } else {
        trace.clear();
    }

    trace
}

#[test]
fn test_traceback() {
    let query = b"ACGTGGA";
    let text = b"TTTTACGTGGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACGTGGATTTTTTT";
    let block_states = compute_traceback(query, text);
    let trace = get_trace(0, 10, 1, &block_states, query.len());
    println!("Trace 1: {:?}", trace);
    let trace = get_trace(0, text.len() - 8, 1, &block_states, query.len());
    println!("Trace 2: {:?}", trace);
}
