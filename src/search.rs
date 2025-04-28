use std::{arch::x86_64::_pext_u64, array::from_fn, iter::zip, simd::Simd};

use pa_types::Cost;

use crate::{
    bitpacking::compute_block_simd,
    delta_encoding::{V, VEncoding},
    profiles::trai_def::Profile,
};

/// Search for query in text.
/// Do something like Farrar's striped SIMD, but on a higher level:
/// - Split the text into 4 equally long 64byte-aligned chunks.
/// - Process those chunks in parallel.
/// - Do at least 2*query.len() overlap between the chunks.
///
/// Returns the results in a reused output vector.
pub fn search<P: Profile>(query: &[u8], text: &[u8], deltas: &mut Vec<V<u64>>) {
    let (profiler, query_profile) = P::encode_query(query);

    assert!(
        query.len() <= 32,
        "For longer queries, we may need more than 64 overlap?"
    );
    // Number of 64char lanes, including 3 overlap between chunks.
    let num_lanes = text.len().div_ceil(64) + 3;
    // Number of simd units to cover everything.
    // TODO: div_ceil
    let num_simds = num_lanes.div_ceil(4);
    // Length of each of the four chunks.
    let chunk_len = num_simds - 1;

    deltas.resize(num_lanes, V::zero());

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

    for i in 0..num_simds {
        // The alignment can start anywhere, so start with deltas of 0.
        let mut vp = S::splat(0);
        let mut vm = S::splat(0);

        for lane in 0..4 {
            let start = lane * chunk_len * 64 + 64 * i;
            if start <= text.len() - 64 {
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

        for (q_char, (hp, hm)) in zip(&query_profile, zip(&mut hp, &mut hm)) {
            let eq = from_fn(|lane| P::eq(&q_char, &text_profile[lane])).into();
            compute_block_simd(hp, hm, &mut vp, &mut vm, eq);
        }
        for lane in 0..4 {
            deltas[lane * chunk_len + i] = <VV as VEncoding<Base>>::from(vp[lane], vm[lane]);
        }
    }
}

pub fn find_below_threshold(
    query: &[u8],
    threshold: Cost,
    deltas: &[V<u64>],
    positions: &mut Vec<usize>,
) {
    let mut cur_cost = query.len() as Cost;
    for (i, v) in deltas.iter().enumerate() {
        let (min, delta) = prefix_min(*v);
        if cur_cost + (min as Cost) <= threshold {
            positions.push(i * 64);
        }
        cur_cost += delta as Cost;
    }
}

/// For each byte: (min_cost, end_cost)
/// Each 1 in a byte indicates -1.
/// Each 0 in a byte indicates +1.
const TABLE: [(i8, i8); 256] = {
    let mut table = [(0, 0); 256];

    let mut i = 0;
    while i < 256 {
        let mut min = 0;
        let mut cur = 0;
        let mut j = 0;
        while j < 8 {
            let bit = (i >> j) & 1;
            let delta = if bit == 1 { -1 } else { 1 };
            cur += delta;
            if cur < min {
                min = cur;
            }
            j += 1;
        }
        table[i] = (min, cur);
        i += 1;
    }

    table
};

/// Return the min_cost in the word, and total cost across the word.
fn prefix_min(v: V<u64>) -> (i8, i8) {
    // extract only the relevant chars
    let (p, m) = v.pm();
    let delta = p | m;
    let num_p = p.count_ones();
    let num_m = m.count_ones();
    let deltas = unsafe { _pext_u64(m, delta) };
    let mut min = 0;
    let mut cur = 0;
    for i in 0..8 {
        let byte = (deltas >> (i * 8)) as u8 as usize;
        let (min_cost, end_cost) = TABLE[byte];
        min = min.min(cur + min_cost);
        cur += end_cost;
    }

    (min, num_p as i8 - num_m as i8)
}

#[test]
fn test_prefix_min() {
    let v = V(0, 0);
    let (min, delta) = prefix_min(v);
    assert_eq!(min, 0);
    assert_eq!(delta, 0);

    let v = V(u64::MAX, 0);
    let (min, delta) = prefix_min(v);
    assert_eq!(min, 0);
    assert_eq!(delta, 64);

    let v = V(0, u64::MAX);
    let (min, delta) = prefix_min(v);
    assert_eq!(min, -64);
    assert_eq!(delta, -64);

    let v = V(
        0b11111111000000001111111100000000,
        0b00000000111111110000000011111111,
    );
    let (min, delta) = prefix_min(v);
    assert_eq!(min, -8);
    assert_eq!(delta, 0);
}

#[test]
fn test_search() {
    let query = b"ACTGNA";

    for len in (256 + 64..512).step_by(21) {
        let mut text = vec![b'A'; len];
        let offset = 10;
        text[offset + 1] = b'C';
        text[offset + 2] = b'T';
        text[offset + 3] = b'G';
        text[offset + 4] = b'N';
        text[offset + 5] = b'A';

        let mut deltas = vec![];
        search::<crate::Iupac>(query, &text, &mut deltas);
        let mut positions = vec![];
        find_below_threshold(query, 2, &deltas, &mut positions);
        // Note that this only returns the start index of the lane; not the exact position.
        assert_eq!(positions, vec![0], "Failure for len {len}");
    }
}
