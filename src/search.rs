use std::{arch::x86_64::_pext_u64, array::from_fn, iter::zip, simd::Simd};

use pa_types::Cost;

use crate::{
    bitpacking::compute_block_simd,
    delta_encoding::{VEncoding, V},
    profile::Profile,
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

    // Number of 64char lanes.
    let num_lanes = text.len().div_ceil(64);
    // Number of simd units to cover everything.
    let num_simds = text.len().div_ceil(256);
    // Length of each of the four chunks.
    // FIXME: overlap.
    let chunk_len = num_lanes.div_ceil(4);

    deltas.resize(num_lanes, V::zero());

    type Base = u64;
    type VV = V<Base>;
    type S = Simd<Base, 4>;

    let mut hp = vec![S::splat(1); query.len()];
    let mut hm = vec![S::splat(0); query.len()];

    let mut text_profile: [_; 4] = Default::default();

    for i in 0..num_simds {
        // The alignment can start anywhere, so start with deltas of 0.
        let mut vp = S::splat(0);
        let mut vm = S::splat(0);

        for lane in 0..4 {
            profiler.encode_ref(
                &text[lane * chunk_len * 64 + 64 * i..][..64]
                    .try_into()
                    .unwrap(),
                &mut text_profile[lane],
            )
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
    deltas: &Vec<V<u64>>,
    positions: &mut Vec<usize>,
) {
    let mut cur_cost = query.len() as Cost;
    for (i, v) in deltas.iter().enumerate() {
        let V(p, m) = v;
        eprintln!("{p:>064b} {m:>064b}");
        let (min, delta) = prefix_min(*v);
        if cur_cost + (min as Cost) <= threshold {
            positions.push(i * 64);
        }
        eprintln!("{}: {} {} {}", i, cur_cost, delta, min);
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

    let mut text = [b'A'; 512];
    text[11] = b'C';
    text[12] = b'T';
    text[13] = b'G';
    text[14] = b'N';
    text[15] = b'A';

    let mut deltas = vec![];
    search::<crate::Iupac>(query, &text, &mut deltas);
    eprintln!("deltas: {deltas:?}");
    let mut positions = vec![];
    find_below_threshold(query, 0, &deltas, &mut positions);
    assert_eq!(positions, vec![5]);
}
