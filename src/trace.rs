use pa_types::Cost;
use pa_types::I;

use crate::bitpacking::compute_block;
use crate::delta_encoding::V;
use crate::delta_encoding::VEncoding;
use crate::profiles::Dna;
use crate::profiles::Profile;

use crate::bitpacking::compute_block_simd;
use std::array::from_fn;
use std::simd::Simd;

/// Costs for states in a single column of an alignment (corresponding to 1 query char vs all the text).
#[derive(Debug, Clone)]
struct ColCosts {
    /// Cost to the top of the col.
    /// Typically simple its index.
    offset: Cost,
    /// Deltas between adjacent rows.
    deltas: Vec<V<u64>>,
}

/// Compute the full n*m matrix corresponding to the query * text alignment.
/// TODO: SIMD variant that takes 1 query, and 4 text slices of the same length.
fn fill(query: &[u8], text: &[u8]) -> Vec<ColCosts> {
    let (profiler, query_profile) = Dna::encode_query(query);
    let mut col_costs: Vec<_> = (0..=query.len())
        .map(|i| ColCosts {
            offset: i as Cost,
            deltas: vec![V::<u64>::zero(); text.len().div_ceil(64)],
        })
        .collect();

    let mut h = vec![(1, 0); query.len()];

    let mut text_profile = profiler.alloc_out();

    // Process chunks of 64 chars, that end exactly at the end of the text.
    for (i, block) in text.rchunks(64).rev().enumerate() {
        let mut slice: [u8; 64] = [b'X'; 64];
        slice[64 - block.len()..].copy_from_slice(block);
        profiler.encode_ref(&slice, &mut text_profile);

        println!("Chunk:: {}", String::from_utf8_lossy(&slice));

        let mut v = V::<u64>::zero();

        for j in 0..query.len() {
            compute_block::<Dna, _, _>(&mut h[j], &mut v, &query_profile[j], &text_profile);
            col_costs[j + 1].deltas[i] = v;
        }
    }

    col_costs
}

fn calc_pad(text_lengths: &[usize]) -> (usize, Vec<usize>) {
    let max_len = text_lengths.iter().max().unwrap().next_multiple_of(64);

    // Calculate how many chars of padding each text needs on the left
    let paddings = text_lengths
        .iter()
        .map(|len| max_len - len) // Simple: difference between max length and this text's length
        .collect();

    (max_len, paddings)
}

fn fill_chunk(
    chunk: &mut [u8; 64],
    text: &[u8],
    text_offset: &mut usize,
    left_padding: &mut usize,
) {
    if *left_padding >= 64 {
        // Full block of padding, leave as all X's
        *left_padding -= 64;
    } else if *left_padding > 0 {
        // Partial padding block: some X's followed by text
        let text_len = 64 - *left_padding;
        let text_end = *text_offset + text_len;
        chunk[*left_padding..].copy_from_slice(&text[*text_offset..text_end]);
        *text_offset = text_end;
        *left_padding = 0;
    } else {
        // Full text block
        let text_end = (*text_offset + 64).min(text.len());
        let text_len = text_end - *text_offset;
        chunk[..text_len].copy_from_slice(&text[*text_offset..text_end]);
        *text_offset = text_end;
    }
}

fn simd_fill<P: Profile>(query: &[u8], texts: [&[u8]; 4]) -> [Vec<ColCosts>; 4] {
    let (profiler, query_profile) = P::encode_query(query);
    let text_lengths = texts.map(|t| t.len());
    let (max_len, mut paddings) = calc_pad(&text_lengths);

    type Base = u64;
    type VV = V<Base>;
    type S = Simd<Base, 4>;

    let mut hp = vec![S::splat(1); query.len()];
    let mut hm = vec![S::splat(0); query.len()];
    let mut text_profile = [
        profiler.alloc_out(),
        profiler.alloc_out(),
        profiler.alloc_out(),
        profiler.alloc_out(),
    ];

    let mut lane_col_costs: [Vec<ColCosts>; 4] = (0..4)
        .map(|_lane| {
            (0..=query.len())
                .map(|i| ColCosts {
                    offset: i as Cost,
                    deltas: vec![V::<u64>::zero(); max_len.div_ceil(64)],
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();

    // Total blocks to calculatep adding
    let total_blocks = max_len.div_ceil(64);
    let mut current_text_offsets = [0; 4];

    // Process blocks left-to-right
    for block_idx in 0..total_blocks {
        for lane in 0..4 {
            let mut chunk = [b'X'; 64];
            fill_chunk(
                &mut chunk,
                texts[lane],
                &mut current_text_offsets[lane],
                &mut paddings[lane],
            );
            println!("Chunk:: {}", String::from_utf8_lossy(&chunk));
            profiler.encode_ref(&chunk, &mut text_profile[lane]);
        }
        let mut vp = S::splat(0);
        let mut vm = S::splat(0);
        for j in 0..query.len() {
            let eq = from_fn(|lane| P::eq(&query_profile[j], &text_profile[lane])).into();
            compute_block_simd(&mut hp[j], &mut hm[j], &mut vp, &mut vm, eq);
            for lane in 0..4 {
                let v = <VV as VEncoding<Base>>::from(vp[lane], vm[lane]);
                lane_col_costs[lane][j + 1].deltas[block_idx] = v;
            }
        }
    }

    lane_col_costs
}

fn get_trace(col_costs: &[ColCosts]) -> Vec<(usize, usize)> {
    let mut trace = Vec::new();
    let mut i = col_costs.len() - 1;
    let mut j = 64 * col_costs[0].deltas.len();

    let cost = |i: usize, j: usize| -> Cost {
        col_costs[i].offset + V::<u64>::value_to(&col_costs[i].deltas, j as I)
    };

    // remaining dist to (i,j)
    let mut g = cost(i, j);

    while i > 0 {
        assert!(j > 0, "Traceback reached top of filled region.");

        trace.push((i, j));

        // Match
        if cost(i - 1, j - 1) == g {
            i -= 1;
            j -= 1;
            println!("Match");
            continue;
        }
        // We make some kind of mutation.
        g -= 1;

        // Insert text char.
        if cost(i, j - 1) == g {
            j -= 1;
            g -= 1;
            println!("Insert");
            continue;
        }
        // Mismatch.
        if cost(i - 1, j - 1) == g {
            i -= 1;
            j -= 1;
            println!("Mismatch");
            continue;
        }
        // Delete query char.
        if cost(i - 1, j) == g {
            i -= 1;
            println!("Delete");
            continue;
        }
        panic!("Trace failed!");
    }

    trace.reverse();

    trace
}

#[test]
fn test_traceback() {
    let query = b"ATTTTCCCGGGGATTTT".as_slice();
    let text2: &[u8] = b"ATTTTGGGGATTTT".as_slice();
    let text1 = b"ATTTTCCCGGGGATTTT".as_slice();
    let text3 = b"TGGGGATTTT".as_slice();
    let text4 = b"TTTTTTTTTTATTTTGGGGATTTT".as_slice();

    let col_costs = simd_fill::<Dna>(&query, [&text1, &text2, &text3, &text4]);
    let trace = get_trace(&col_costs[0]);
    println!("Trace: {:?}", trace);
}

// let text1 = b"ATCGACTAGC".as_slice();

// let text3 = b"CTAGC".as_slice();
// let text4 = b"TGGC".as_slice();

// let col_costs = fill(query, text2);
// for c in col_costs {
//     println!("col: {:?}", c);
//     let (p, m) = c.deltas[0].pm();
//     println!("p: {:064b} \nm: {:064b}\n\n", p, m);
// }

// let col_costs = simd_fill::<Iupac>(&query, [&text2, &text2, &text2, &text2]);
// let c = &col_costs[0];
// for col in c {
//     let (p, m) = col.deltas[0].pm();
//     // println!("(p, m): {:?}", (p, m));
//     // print the binary 0100101
//     println!("p: {:064b} \nm: {:064b}\n\n", p, m);
// }

// // Simd
// let col_costs = simd_fill::<Dna>(&query, [&text2, &text2, &text2, &text2]);

// for lane in 0..4 {
//     //println!("\nCol costs for lane {}\n{:?}", lane, col_costs[lane]);
//     let trace = get_trace(&col_costs[lane]);
//     println!("Trace {}: {:?}", lane, trace);
// }

// #[test]
// fn test_and_block_boundary() {
//     let query = b"ACGTGGA";
//     let mut text = [b'G'; 128];
//     text[64 - 3..64 + 4].copy_from_slice(query);
//     let col_costs = fill(query, &text[..64 + 4]);
//     let trace = get_trace(&col_costs);
//     println!("Trace 1: {:?}", trace); // FIXME: This is wrong when crossing block boundary
// }
/*
query:   ATTTTCCCGGGGATTTT
text2: ...GGATTTTCCGGATTTT
 */
