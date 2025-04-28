use std::{array::from_fn, iter::zip, simd::Simd};

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
pub fn search<P: Profile>(query: &[u8], text: &[u8], output: &mut Vec<V<u64>>) {
    let (profiler, query_profile) = P::encode_query(query);

    // Number of 64char lanes.
    let num_lanes = text.len().div_ceil(64);
    // Number of simd units to cover everything.
    let num_simds = text.len().div_ceil(256);
    // Length of each of the four chunks.
    // FIXME: overlap.
    let chunk_len = num_lanes.div_ceil(4);

    output.resize(num_lanes, V::zero());

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
            output[lane * chunk_len + i] = <VV as VEncoding<Base>>::from(vp[lane], vm[lane]);
        }
    }
}

#[test]
fn test_search() {
    let query = b"ACTGNA";

    let text = [b'A'; 512];
    search::<crate::Iupac>(query, &text, &mut vec![]);
}
