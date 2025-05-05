use crate::minima::prefix_min;
use pa_types::Cost;
use std::{
    array::from_fn,
    simd::{Simd, cmp::SimdPartialOrd},
};

use crate::{
    bitpacking::compute_block_simd,
    delta_encoding::{V, VEncoding},
    profiles::Profile,
};

/// Search for query in text.
/// Do something like Farrar's striped SIMD, but on a higher level:
/// - Split the text into 4 equally long 64byte-aligned chunks.
/// - Process those chunks in parallel.
/// - Do at least 2*query.len() overlap between the chunks.
///
/// Returns the results in a reused output vector.
// FIXME: We really need some proper tests for this.
pub fn search_positions<P: Profile>(query: &[u8], text: &[u8], deltas: &mut Vec<V<u64>>) {
    search_positions_maybe_bounded::<P, false>(query, text, 0, deltas);
}

/// Like [`search`], but early-abort upon reaching a threshold of `k`.
// FIXME: We really need some proper tests for this.
pub fn search_positions_bounded<P: Profile>(
    query: &[u8],
    text: &[u8],
    k: Cost,
    deltas: &mut Vec<V<u64>>,
) {
    search_positions_maybe_bounded::<P, true>(query, text, k, deltas);
}

fn search_positions_maybe_bounded<P: Profile, const BOUNDED: bool>(
    query: &[u8],
    text: &[u8],
    k: Cost,
    deltas: &mut Vec<V<u64>>,
) {
    let (profiler, query_profile) = P::encode_query(query);

    // Terminology:
    // - chunk: roughly 1/4th of the input text, with small overlaps.
    // - block: 64 bytes of text.
    // - lane: a u64 of a SIMD vec.

    // The query will match a pattern of length at most query.len() + k.
    // We round that up to a multiple of 64 to find the number of blocks overlap between chunks.
    let overlap_blocks = (query.len() + k as usize).next_multiple_of(64) / 64;

    // Total number of blocks to be processed, including overlaps.
    let text_blocks = text.len().div_ceil(64);
    let total_blocks = text_blocks + 3 * overlap_blocks;
    let blocks_per_chunk = total_blocks.div_ceil(4);
    // Length of each of the four chunks.
    let chunk_offset = blocks_per_chunk.saturating_sub(overlap_blocks);
    deltas.resize(text_blocks, V::zero());

    type Base = u64;
    type VV = V<Base>;
    type S = Simd<Base, 4>;

    let mut hp = vec![S::splat(1); query.len()];
    let mut hm = vec![S::splat(0); query.len()];

    let mut text_profile: [P::B; 4] = [
        profiler.alloc_out(),
        profiler.alloc_out(),
        profiler.alloc_out(),
        profiler.alloc_out(),
    ];

    let mut text_chunks: [[u8; 64]; 4] = [[0; 64]; 4];

    let mut prev_max_j = 0;

    for i in 0..blocks_per_chunk {
        // The alignment can start anywhere, so start with deltas of 0.
        let mut vp = S::splat(0);
        let mut vm = S::splat(0);

        let mut x = S::splat(0);

        // Collect the 4 slices of input text.
        // Out-of-bounds characters are replaced by 'X', which doesn't match anything.
        for lane in 0..4 {
            let start = lane * chunk_offset * 64 + 64 * i;
            if start + 64 <= text.len() {
                text_chunks[lane] = text[start..start + 64].try_into().unwrap();
            } else {
                text_chunks[lane] = [b'X'; 64];
                if start <= text.len() {
                    let slice = &text[start..];
                    text_chunks[lane][..slice.len()].copy_from_slice(slice);
                }
            }
            profiler.encode_ref(&text_chunks[lane], &mut text_profile[lane])
        }

        let mut dist_to_start_of_lane = S::splat(0);
        let mut dist_to_end_of_lane = S::splat(0);

        let mut cur_max_j = 0;

        // Iterate over query chars.
        'query: for j in 0..query.len() {
            dist_to_start_of_lane += hp[j];
            dist_to_start_of_lane -= hm[j];

            let eq = from_fn(|lane| P::eq(&query_profile[j], &text_profile[lane])).into();
            compute_block_simd(&mut hp[j], &mut hm[j], &mut vp, &mut vm, eq);

            if BOUNDED {
                // For DNA, the distance between random/unrelated sequences is around q.len()/2.
                // Thus, for threshold k, we can expect random matches between seqs of length ~2*k.
                // To have some buffer, we start filtering at length 3*k.
                //
                // TODO: Currently this filtering is much too slow to be useable for small k.
                if j >= 16 && j > prev_max_j && j.is_power_of_two() {
                    // Check for each lane
                    for lane in 0..4 {
                        let v = V(vp.as_array()[lane], vm.as_array()[lane]);
                        let min_in_lane = dist_to_start_of_lane.as_array()[lane] as Cost
                            + prefix_min(v).0 as Cost;
                        if min_in_lane <= k {
                            continue 'query;
                        }
                    }
                    // eprintln!("Break after col {j}");
                    // All lanes only have values > k. We set remaining horizontal deltas to +1.
                    for j2 in j + 1..query.len() {
                        hp[j2] = S::splat(1);
                        hm[j2] = S::splat(0);
                    }
                    break;
                }

                dist_to_end_of_lane += hp[j];
                dist_to_end_of_lane -= hm[j];
                // FIXME: Somehwo the deltas are ocasionally off by 1. We need
                // this +1 here to fix things.
                // This feels like a bug in bitpacking.rs.
                cur_max_j = if (dist_to_end_of_lane.simd_le(S::splat(k as u64 + 1))).any() {
                    j
                } else {
                    cur_max_j
                };
            }
        }
        for lane in 0..4 {
            let idx = lane * chunk_offset + i;
            if idx < deltas.len() {
                deltas[idx] = <VV as VEncoding<Base>>::from(vp[lane], vm[lane]);
            }
        }
        prev_max_j = cur_max_j;
    }
}
