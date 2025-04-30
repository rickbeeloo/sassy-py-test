use pa_types::Cost;
use pa_types::I;

use crate::bitpacking::compute_block;
use crate::delta_encoding::VEncoding;
use crate::delta_encoding::V;
use crate::profiles::Dna;
use crate::profiles::Profile;

struct ColCosts {
    offset: Cost,
    deltas: Vec<V<u64>>,
}

fn compute_traceback(query: &[u8], text: &[u8]) -> Vec<ColCosts> {
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

        let mut v = V::<u64>::zero();

        for j in 0..query.len() {
            compute_block::<Dna, _, _>(&mut h[j], &mut v, &query_profile[j], &text_profile);
            col_costs[j + 1].deltas[i] = v.clone();
        }
    }

    col_costs
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
            continue;
        }
        // We make some kind of mutation.
        g -= 1;

        // Insert text char.
        if cost(i, j - 1) == g {
            j -= 1;
            g -= 1;
            continue;
        }
        // Mismatch.
        if cost(i - 1, j - 1) == g {
            i -= 1;
            j -= 1;
            continue;
        }
        // Delete query char.
        if cost(i - 1, j) == g {
            i -= 1;
            continue;
        }
        panic!("Trace failed!");
    }

    trace.reverse();

    trace
}

#[test]
fn test_traceback() {
    let query = b"ACGTGGA";
    let text = b"TTTTACGTGGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACGTGGATTTTTTT";
    let end = 11;
    let col_costs = compute_traceback(query, &text[..end]);
    let trace = get_trace(&col_costs);
    println!("Trace 1: {:?}", trace);
    //  [(6, 10), (5, 9), (4, 8), (3, 7), (2, 5), (1, 5), (0, 3)]
    // should be [(6, 10), (5, 9), (4, 8), (3, 7), (2, 6), (1, 5), (0, 4)]
    let col_costs = compute_traceback(query, &text[..text.len() - 7]);
    let trace = get_trace(&col_costs);
    println!("Trace 2: {:?}", trace); // FIXME: This is wrong
}

#[test]
fn test_and_block_boundary() {
    let query = b"ACGTGGA";
    let mut text = [b'G'; 128];
    text[64 - 3..64 + 4].copy_from_slice(query);
    let col_costs = compute_traceback(query, &text[..64 + 4]);
    let trace = get_trace(&col_costs);
    println!("Trace 1: {:?}", trace); // FIXME: This is wrong when crossing block boundary
}
