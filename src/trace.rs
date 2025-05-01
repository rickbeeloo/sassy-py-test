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
    for (i, block) in text.chunks(64).enumerate() {
        let mut slice: [u8; 64] = [b'X'; 64];
        slice[..block.len()].copy_from_slice(block);
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

fn simd_fill<P: Profile>(query: &[u8], texts: [&[u8]; 4]) -> [Vec<ColCosts>; 4] {
    let (profiler, query_profile) = P::encode_query(query);
    let max_len = texts.iter().map(|t| t.len()).max().unwrap();
    let num_chunks = max_len.div_ceil(64);

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

    let mut lane_col_costs: [Vec<ColCosts>; 4] = from_fn(|_lane| {
        (0..=query.len())
            .map(|i| ColCosts {
                offset: i as Cost,
                deltas: vec![V::<u64>::zero(); num_chunks],
            })
            .collect::<Vec<_>>()
    });

    for i in 0..num_chunks {
        for lane in 0..4 {
            let mut slice = [b'X'; 64];
            let block = texts[lane].get(64 * i..).unwrap_or_default();
            let block = block.get(..64).unwrap_or(block);
            slice[..block.len()].copy_from_slice(block);
            profiler.encode_ref(&slice, &mut text_profile[lane]);
        }
        let mut vp = S::splat(0);
        let mut vm = S::splat(0);
        for j in 0..query.len() {
            let eq = from_fn(|lane| P::eq(&query_profile[j], &text_profile[lane])).into();
            compute_block_simd(&mut hp[j], &mut hm[j], &mut vp, &mut vm, eq);
            for lane in 0..4 {
                let v = <VV as VEncoding<Base>>::from(vp[lane], vm[lane]);
                lane_col_costs[lane][j + 1].deltas[i] = v;
            }
        }
    }

    lane_col_costs
}

fn get_trace<P: Profile>(query: &[u8], text: &[u8], col_costs: &[ColCosts]) -> Vec<(usize, usize)> {
    let mut trace = Vec::new();
    let mut i = query.len();
    let mut j = text.len();

    let cost = |i: usize, j: usize| -> Cost {
        col_costs[i].offset + V::<u64>::value_to(&col_costs[i].deltas, j as I)
    };

    // remaining dist to (i,j)
    let mut g = cost(i, j);

    while i > 0 {
        eprintln!("({i}, {j}) {g}");
        trace.push((i, j));

        // Match
        if j > 0 && cost(i - 1, j - 1) == g && P::is_match(query[i - 1], text[j - 1]) {
            eprintln!("match");
            i -= 1;
            j -= 1;
            continue;
        }
        // We make some kind of mutation.
        g -= 1;

        // Insert text char.
        if j > 0 && cost(i, j - 1) == g {
            eprintln!("insert");
            j -= 1;
            continue;
        }
        // Mismatch.
        if j > 0 && cost(i - 1, j - 1) == g {
            eprintln!("mismatch");
            i -= 1;
            j -= 1;
            continue;
        }
        // Delete query char.
        if cost(i - 1, j) == g {
            eprintln!("delete");
            i -= 1;
            continue;
        }
        panic!(
            "Trace failed! No ancestor found of {i} {j} at distance {}",
            g + 1
        );
    }

    assert_eq!(g, 0, "Remaining cost after the trace must be 0.");

    trace.reverse();

    trace
}

#[test]
fn test_traceback() {
    let query = b"ATTTTCCCGGGGATTTT".as_slice();
    let text2: &[u8] = b"ATTTTGGGGATTTT".as_slice();

    let col_costs = fill(query, text2);
    for c in &col_costs {
        let (p, m) = c.deltas[0].pm();
        println!("{}\np: {:064b} \nm: {:064b}\n", c.offset, p, m);
    }

    let trace = get_trace::<Dna>(query, text2, &col_costs);
    println!("Trace: {:?}", trace);
}

#[test]
fn test_traceback_simd() {
    let query = b"ATTTTCCCGGGGATTTT".as_slice();
    let text1 = b"ATTTTCCCGGGGATTTT".as_slice();
    let text2 = b"ATTTTGGGGATTTT".as_slice();
    let text3 = b"TGGGGATTTT".as_slice();
    let text4 = b"TTTTTTTTTTATTTTGGGGATTTT".as_slice();

    let col_costs = simd_fill::<Dna>(&query, [&text1, &text2, &text3, &text4]);
    let _trace = get_trace::<Dna>(query, text1, &col_costs[0]);
    let _trace = get_trace::<Dna>(query, text2, &col_costs[1]);
    let _trace = get_trace::<Dna>(query, text3, &col_costs[2]);
    let trace = get_trace::<Dna>(query, text4, &col_costs[3]);
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
