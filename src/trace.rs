use pa_types::Cigar;
use pa_types::Cost;
use pa_types::I;
use pa_types::Pos;

use crate::bitpacking::compute_block;
use crate::delta_encoding::V;
use crate::delta_encoding::VEncoding;
use crate::profiles::Dna;
use crate::profiles::Profile;

use crate::Match;
use crate::bitpacking::compute_block_simd;
use std::array::from_fn;
use std::simd::Simd;

#[derive(Debug, Clone, Default)]
pub struct CostMatrix {
    /// Query lenght.
    q: usize,
    deltas: Vec<V<u64>>,
}

impl CostMatrix {
    /// i: text idx
    /// j: query idx
    pub fn get(&self, i: usize, j: usize) -> Cost {
        let mut s = j as Cost;
        for idx in (j..j + i / 64 * (self.q + 1)).step_by(self.q + 1) {
            s += self.deltas[idx].value();
        }
        if i % 64 != 0 {
            s += self.deltas[j + i / 64 * (self.q + 1)].value_of_prefix(i as I % 64);
        }
        s
    }
}

/// Compute the full n*m matrix corresponding to the query * text alignment.
/// TODO: SIMD variant that takes 1 query, and 4 text slices of the same length.
#[allow(unused)] // FIXME
fn fill(query: &[u8], text: &[u8], m: &mut CostMatrix) {
    m.q = query.len();
    m.deltas.clear();
    m.deltas.reserve((m.q + 1) * text.len().div_ceil(64));
    let (profiler, query_profile) = Dna::encode_query(query);
    let mut h = vec![(1, 0); query.len()];

    let mut text_profile = profiler.alloc_out();

    // Process chunks of 64 chars, that end exactly at the end of the text.
    for (i, block) in text.chunks(64).enumerate() {
        let mut slice: [u8; 64] = [b'X'; 64];
        slice[..block.len()].copy_from_slice(block);
        profiler.encode_ref(&slice, &mut text_profile);

        println!("Chunk:: {}", String::from_utf8_lossy(&slice));

        let mut v = V::<u64>::zero();

        m.deltas.push(v);
        for j in 0..query.len() {
            compute_block::<Dna, _, _>(&mut h[j], &mut v, &query_profile[j], &text_profile);
            m.deltas.push(v);
        }
    }
}

pub fn simd_fill<P: Profile>(query: &[u8], texts: &[&[u8]], m: &mut [CostMatrix; 4]) {
    assert!(texts.len() <= 4);
    let lanes = texts.len();

    let (profiler, query_profile) = P::encode_query(query);
    let max_len = texts.iter().map(|t| t.len()).max().unwrap();
    let num_chunks = max_len.div_ceil(64);

    for m in &mut *m {
        m.q = query.len();
        m.deltas.clear();
        m.deltas.reserve((m.q + 1) * num_chunks);
    }

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

    for i in 0..num_chunks {
        for lane in 0..lanes {
            let mut slice = [b'X'; 64];
            let block = texts[lane].get(64 * i..).unwrap_or_default();
            let block = block.get(..64).unwrap_or(block);
            slice[..block.len()].copy_from_slice(block);
            profiler.encode_ref(&slice, &mut text_profile[lane]);
        }
        let mut vp = S::splat(0);
        let mut vm = S::splat(0);
        for lane in 0..lanes {
            let v = <VV as VEncoding<Base>>::from(vp[lane], vm[lane]);
            m[lane].deltas.push(v);
        }
        // FIXME: for large queries, use the SIMD within this single block, rather than spreading it thin over 4 'matches' when there is only a single candidate match.
        for j in 0..query.len() {
            let eq = from_fn(|lane| P::eq(&query_profile[j], &text_profile[lane])).into();
            compute_block_simd(&mut hp[j], &mut hm[j], &mut vp, &mut vm, eq);
            for lane in 0..lanes {
                let v = <VV as VEncoding<Base>>::from(vp[lane], vm[lane]);
                m[lane].deltas.push(v);
            }
        }
    }

    for lane in 0..lanes {
        assert_eq!(m[lane].deltas.len(), num_chunks * (m[lane].q + 1));
    }
}

pub fn get_trace<P: Profile>(
    query: &[u8],
    text_offset: usize,
    text: &[u8],
    m: &CostMatrix,
) -> Match {
    let mut trace = Vec::new();
    let mut i = query.len();
    let mut j = text.len();

    let cost = |i: usize, j: usize| -> Cost { m.get(j, i) };

    // remaining dist to (i,j)
    let mut g = cost(i, j);
    let total_cost = g;

    let mut cigar = Cigar::default();

    loop {
        // eprintln!("({i}, {j}) {g}");
        trace.push((i, text_offset + j));

        if i == 0 {
            break;
        }

        // Match
        if j > 0 && cost(i - 1, j - 1) == g && P::is_match(query[i - 1], text[j - 1]) {
            cigar.push(pa_types::CigarOp::Match);
            i -= 1;
            j -= 1;
            continue;
        }
        // We make some kind of mutation.
        g -= 1;

        // Insert text char.
        if j > 0 && cost(i, j - 1) == g {
            cigar.push(pa_types::CigarOp::Ins);
            j -= 1;
            continue;
        }
        // Mismatch.
        if j > 0 && cost(i - 1, j - 1) == g {
            cigar.push(pa_types::CigarOp::Sub);
            i -= 1;
            j -= 1;
            continue;
        }
        // Delete query char.
        if cost(i - 1, j) == g {
            cigar.push(pa_types::CigarOp::Del);
            i -= 1;
            continue;
        }
        panic!(
            "Trace failed! No ancestor found of {i} {j} at distance {}",
            g + 1
        );
    }

    assert_eq!(g, 0, "Remaining cost after the trace must be 0.");

    cigar.reverse();

    Match {
        cost: total_cost,
        start: Pos(0, (text_offset + j) as I),
        end: Pos(query.len() as I, (text_offset + text.len()) as I),
        strand: crate::Strand::Fwd,
        cigar,
    }
}

#[test]
fn test_traceback() {
    let query = b"ATTTTCCCGGGGATTTT".as_slice();
    let text2: &[u8] = b"ATTTTGGGGATTTT".as_slice();

    let mut cost_matrix = Default::default();
    fill(query, text2, &mut cost_matrix);

    let trace = get_trace::<Dna>(query, 0, text2, &cost_matrix);
    println!("Trace: {:?}", trace);
}

#[test]
fn test_traceback_simd() {
    let query = b"ATTTTCCCGGGGATTTT".as_slice();
    let text1 = b"ATTTTCCCGGGGATTTT".as_slice();
    let text2 = b"ATTTTGGGGATTTT".as_slice();
    let text3 = b"TGGGGATTTT".as_slice();
    let text4 = b"TTTTTTTTTTATTTTGGGGATTTT".as_slice();

    let mut cost_matrix = Default::default();
    simd_fill::<Dna>(&query, &[&text1, &text2, &text3, &text4], &mut cost_matrix);
    let _trace = get_trace::<Dna>(query, 0, text1, &cost_matrix[0]);
    let _trace = get_trace::<Dna>(query, 0, text2, &cost_matrix[1]);
    let _trace = get_trace::<Dna>(query, 0, text3, &cost_matrix[2]);
    let trace = get_trace::<Dna>(query, 0, text4, &cost_matrix[3]);
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
