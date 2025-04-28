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
pub fn search<P: Profile>(query: &[u8], text: &[u8]) -> Vec<V<u64>> {
    let query_profile = P::encode_a(query);

    // Number of 64char lanes.
    let num_lanes = text.len().div_ceil(64);
    // Number of simd units to cover everything.
    let num_simds = text.len().div_ceil(256);
    // Length of each of the four chunks.
    // TODO: overlap.
    let chunk_len = num_lanes.div_ceil(4);

    let mut output = vec![V::zero(); num_lanes];

    type Base = u64;
    type VV = V<Base>;
    type S = Simd<Base, 4>;

    let mut hp = vec![S::splat(1); query.len()];
    let mut hm = vec![S::splat(0); query.len()];

    for i in 0..num_simds {
        // The alignment can start anywhere, so start with deltas of 0.
        let mut vp = S::splat(0);
        let mut vm = S::splat(0);

        let text_profile: [_; 4] = from_fn(|lane| {
            P::encode_b(
                &text[lane * chunk_len * 64 + 64 * i..][..64]
                    .try_into()
                    .unwrap(),
            )
        });

        for (q_char, (hp, hm)) in zip(&query_profile, zip(&mut hp, &mut hm)) {
            let eq = from_fn(|lane| P::eq(&q_char, &text_profile[lane])).into();
            compute_block_simd(hp, hm, &mut vp, &mut vm, eq);
        }
        for lane in 0..4 {
            output[lane * chunk_len + i] = <VV as VEncoding<Base>>::from(vp[lane], vm[lane]);
        }
    }
    output
}
