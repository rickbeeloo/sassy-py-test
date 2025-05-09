#![feature(portable_simd, int_roundings, let_chains)]
mod bitpacking;
mod delta_encoding;
pub mod profiles {

    mod ascii;
    mod dna;
    mod iupac;
    mod profile;

    pub use profile::Profile;

    pub use crate::LocalMinK;
    pub use ascii::{Ascii, CaseInsensitiveAscii, CaseSensitiveAscii};
    pub use dna::Dna;
    pub use iupac::Iupac;
}

pub mod implementations {
    pub mod crispr;
    pub mod search;
}

mod minima;
mod search;
mod trace;

pub use minima::{find_all_minima, find_below_threshold, find_local_minima};
use pa_types::{Cigar, Cost, Pos};
use search::Deltas;
pub use search::{search_positions, search_positions_bounded};
use std::{array::from_fn, cell::RefCell, simd::Simd};

use profiles::Profile;
use trace::{CostMatrix, get_trace, simd_fill};

#[cfg(feature = "avx512")]
const LANES: usize = 8;
#[cfg(not(feature = "avx512"))]
const LANES: usize = 4;
type S = Simd<u64, LANES>;

#[derive(Debug, Clone)]
pub struct Match {
    pub start: Pos,
    pub end: Pos,
    pub cost: Cost,
    pub strand: Strand,
    pub cigar: Cigar,
}

#[derive(Debug, Clone, Copy)]
pub enum Strand {
    Fwd,
    Rc,
}

pub fn search<P: Profile>(query: &[u8], text: &[u8], k: usize) -> Vec<Match> {
    let mut deltas = vec![];
    search_positions::<P>(query, text, &mut deltas);

    let matches = find_local_minima(query, &mut deltas, k as Cost, text.len());

    let mut traces = Vec::with_capacity(matches.len());

    thread_local! {
        static M: RefCell<[CostMatrix;LANES]> = RefCell::new(from_fn(|_| CostMatrix::default()));
    }

    let fill_len = query.len() + k;
    for matches in matches.chunks(LANES) {
        let mut text_slices = [[].as_slice(); LANES];
        let mut offsets = [0; LANES];
        for i in 0..matches.len() {
            let end_pos = matches[i].0;
            let offset = end_pos.saturating_sub(fill_len);
            offsets[i] = offset;
            text_slices[i] = &text[offset..end_pos];
        }
        let text_slices = &text_slices[..matches.len()];
        M.with(|m| {
            let mut m = m.borrow_mut();
            simd_fill::<P>(query, text_slices, &mut m);

            for lane in 0..matches.len() {
                let m = get_trace::<P>(query, offsets[lane], text_slices[lane], &m[lane]);
                assert!(
                    m.cost <= k as Cost,
                    "Match has cost {} > {}: {m:?}",
                    m.cost,
                    k
                );
                traces.push(m);
            }
        });
    }

    traces
}

pub fn search_bounded<P: Profile>(query: &[u8], text: &[u8], k: usize) -> Vec<Match> {
    let mut deltas = vec![];
    search_positions_bounded::<P>(query, text, k as Cost, &mut deltas);

    let matches = find_local_minima(query, &mut deltas, k as Cost, text.len());

    let mut traces = Vec::with_capacity(matches.len());

    thread_local! {
        static M: RefCell<[CostMatrix;LANES]> = RefCell::new(from_fn(|_| CostMatrix::default()));
    }

    let fill_len = query.len() + k;
    for matches in matches.chunks(LANES) {
        let mut text_slices = [[].as_slice(); LANES];
        let mut offsets = [0; LANES];
        for i in 0..matches.len() {
            let end_pos = matches[i].0;
            let offset = end_pos.saturating_sub(fill_len);
            offsets[i] = offset;
            text_slices[i] = &text[offset..end_pos];
        }
        let text_slices = &text_slices[..matches.len()];

        M.with(|m| {
            let mut m = m.borrow_mut();
            simd_fill::<P>(query, text_slices, &mut m);

            for lane in 0..matches.len() {
                let m = get_trace::<P>(query, offsets[lane], text_slices[lane], &m[lane]);
                assert!(
                    m.cost <= k as Cost,
                    "Match has cost {} > {k}: {m:?}",
                    m.cost,
                );
                traces.push(m);
            }
        });
    }
    traces
}

pub fn search_with_rc<P: Profile>(query: &[u8], text: &[u8], k: usize) -> Vec<Match> {
    let mut matches = search::<P>(query, text, k);
    let query_rc = &P::reverse_complement(query);
    let mut rc_matches = search::<P>(query_rc, text, k);
    // Patch up the rc matches.
    for m in &mut rc_matches {
        m.strand = Strand::Rc;
    }
    matches.extend(rc_matches);
    matches
}

pub fn search_maybe_rc<P: Profile>(query: &[u8], text: &[u8], k: usize, rc: bool) -> Vec<Match> {
    if rc {
        search_with_rc::<P>(query, text, k)
    } else {
        search::<P>(query, text, k)
    }
}

/*

    Below mostly to benchmark the different behaviors, see examples/bench.rs
    mmaking sure all options are commpiled
*/

pub struct WithRc;
pub struct WithoutRc;
pub struct Bounded;
pub struct Unbounded;

pub trait BoundBehavior<P: Profile> {
    fn search(query: &[u8], text: &[u8], k: usize) -> Vec<Match>;
    fn search_positions(query: &[u8], text: &[u8], k: Cost, deltas: &mut Deltas);
}

pub trait RcBehavior<P: Profile> {
    fn search<F: Fn(&[u8], &[u8], usize) -> Vec<Match>>(
        query: &[u8],
        text: &[u8],
        k: usize,
        search_fn: F,
    ) -> Vec<Match>;
}

pub struct AllK;
pub struct LocalMinK;

pub trait MinimaBehavior {
    fn find_minima(
        query: &[u8],
        deltas: &mut Deltas,
        k: Cost,
        text_len: usize,
    ) -> Vec<(usize, Cost)>;
}

impl MinimaBehavior for AllK {
    fn find_minima(
        query: &[u8],
        deltas: &mut Deltas,
        k: Cost,
        text_len: usize,
    ) -> Vec<(usize, Cost)> {
        find_all_minima(query, deltas, k, text_len)
    }
}

impl MinimaBehavior for LocalMinK {
    fn find_minima(
        query: &[u8],
        deltas: &mut Deltas,
        k: Cost,
        text_len: usize,
    ) -> Vec<(usize, Cost)> {
        find_local_minima(query, deltas, k, text_len)
    }
}

pub fn search_generic<P: Profile, RcOpt, BoundOpt, MinOpt>(
    query: &[u8],
    text: &[u8],
    k: usize,
) -> Vec<Match>
where
    RcOpt: RcBehavior<P>,
    BoundOpt: BoundBehavior<P>,
    MinOpt: MinimaBehavior,
{
    RcOpt::search(query, text, k, |q, t, k| {
        let mut deltas = vec![];
        BoundOpt::search_positions(q, t, k as Cost, &mut deltas);
        let matches = MinOpt::find_minima(q, &mut deltas, k as Cost, t.len());

        let mut traces = Vec::with_capacity(matches.len());

        thread_local! {
            static M: RefCell<[CostMatrix;LANES]> = RefCell::new(from_fn(|_| CostMatrix::default()));
        }

        let fill_len = q.len() + k;

        for matches in matches.chunks(LANES) {
            let mut text_slices = [[].as_slice(); LANES];
            let mut offsets = [0; LANES];
            for i in 0..matches.len() {
                let end_pos = matches[i].0;
                let offset = end_pos.saturating_sub(fill_len);
                offsets[i] = offset;
                text_slices[i] = &t[offset..end_pos];
            }
            let text_slices = &text_slices[..matches.len()];
            M.with(|m| {
                let mut m = m.borrow_mut();
                simd_fill::<P>(q, text_slices, &mut m);

                for lane in 0..matches.len() {
                    let m = get_trace::<P>(q, offsets[lane], text_slices[lane], &m[lane]);
                    assert!(
                        m.cost <= k as Cost,
                        "Match has cost {} > {}: {m:?}",
                        m.cost,
                        k
                    );
                    traces.push(m);
                }
            });
        }

        traces
    })
}

// Compiled version without reverse complement search
impl<P: Profile> RcBehavior<P> for WithoutRc {
    fn search<F: Fn(&[u8], &[u8], usize) -> Vec<Match>>(
        query: &[u8],
        text: &[u8],
        k: usize,
        search_fn: F,
    ) -> Vec<Match> {
        search_fn(query, text, k)
    }
}

// Compiled version with reverse complement search
impl<P: Profile> RcBehavior<P> for WithRc {
    fn search<F: Fn(&[u8], &[u8], usize) -> Vec<Match>>(
        query: &[u8],
        text: &[u8],
        k: usize,
        search_fn: F,
    ) -> Vec<Match> {
        let mut matches = search_fn(query, text, k);
        let query_rc = &P::reverse_complement(query);
        let mut rc_matches = search_fn(query_rc, text, k);
        for m in &mut rc_matches {
            m.strand = Strand::Rc;
        }
        matches.extend(rc_matches);
        matches
    }
}

// Compiled version without bounds
impl<P: Profile> BoundBehavior<P> for Unbounded {
    fn search(query: &[u8], text: &[u8], k: usize) -> Vec<Match> {
        search::<P>(query, text, k)
    }
    fn search_positions(query: &[u8], text: &[u8], k: Cost, deltas: &mut Deltas) {
        search_positions::<P>(query, text, deltas)
    }
}

// Compiled version with bounds
impl<P: Profile> BoundBehavior<P> for Bounded {
    fn search(query: &[u8], text: &[u8], k: usize) -> Vec<Match> {
        search_bounded::<P>(query, text, k)
    }
    fn search_positions(query: &[u8], text: &[u8], k: Cost, deltas: &mut Deltas) {
        search_positions_bounded::<P>(query, text, k, deltas)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::random_range;

    use crate::profiles::Dna;

    #[test]
    fn test_case1() {
        let query = b"AGATGTGTCC";
        let text = b"GAGAGATAACCGTGCGCTCACTGTTACAGTTTATGTGTCGAATTCTTTTAGGGCTGGTCACTGCCCATGCGTAGGAATGAATAGCGGTTGCGGATAACTAAGCAGTGCTTGTCGCTATTAAAGTTAGACCCCGCTCCCACCTTGCCACTCCTAGATGTCGACGAGATTCTCACCGCCAAGGGATCAGGCATATCATAGTAGCTGGCGCAGCCCGTCGTTTAAGGAGGCCCATATACTAATTCAAACGATGGGGTGGGCACATCCCCTAGAAGCACTTGTCCCTTGAGTCACGACTGATGCGTGGTCTCTCGCTAAATGTTCCGGCCTCTCGGACATTTAAAGGGTGGCATGTGACCATGGAGGATTAGTGAATGAGAGGTGTCCGCCTTTGTTCGCCAGAACTCTTATAGCGTAGGGGAGTGTACTCACCGCGAACCCGTATCAGCAATCTTGTCAGTGGCTCCTGACTCAAACATCGATGCGCTGCACATGGCCTTAGAATGAAGCAGCCCCTTTCTATTGTGGCCGGGCTGATTCTTTGTTTTGTTGAATGGTCGGGCCTGTCTGCCTTTTCCTAGTGTTGAAACTCCGAACCGCATGAACTGCGTTGCTAGCGAGCTATCACTGGACTGGCCGGGGGACGAAAGTTCGCGGGACCCACTACCCGCGCCCAGAAGACCACACTAGGGAGAAGGATTCTATCGGCATAGCCGTC";
        let matches = search::<Dna>(query, text, 2);
        println!("Matches: {:?}", matches);
    }

    #[test]
    fn test_case2() {
        let query = b"CGCAGAAGAA";
        let text = b"CGTCCGAGCTTACGACGTAGGTCAATTTTCCGTCCATTTCCAAACAGCATTTATGAGATCCATTATGGATGTACATCAGAATCATAGGCAAATGAGCCTGCGGCGACTTTAGTCGCAAGCGAATGCGCCGATGTAATCGGGACATAGAGAAGCGGGCGGCTACTAAGCTCCTCTTAGTGACGTATCGGCGACGTGTGGCGTTCATACTCTGGGCCAGGGTTGAAACGTTTATGACTCTAACGAGCCGGATGCCGCTTGCGAAAGCGGAGTTGATAGAACATGGGAACGTGGTATGGTCATCTATCAGTTAGGAGGGTTTAACCTTGGTTTAATCAGGGTTGTCCCTGGCTACTTGAGAGACGCACCGCTCTAACTATAAAGGGGAATCTCGCTTCTCCGGACTTTCCCAACACTGTACGATAAGGGAGCCAAGTGAAGTTCGGAGCTATTGTACTCTTCTATAGTTATCACCCTCACGGGGGGACTAAAAATCAGAATCTATTGAGATCCTATACGCTCGCAAGTGCTTCCTGAGTCCGATTCGGGTGGGAGAGAGATGAAGGCACCTCGCCGCCGTTCTGCCGAAGACGAGAACCCCCCGCGCTTTGTCTTAATCTAAGTGCGGCCATTACGACCGAAGACGGAGTGTATGGTACTGCGACCTACAGATAGGGCTAGCCGCATGAGAAAATGCCATCTAGACAAGGCTACGTGTGCTACTTAGCAATTCTGTTTCACCTCGAATTGTATACACGGCCTGCAAAACCAAGCGGAGCTCTGTAATAGGAGTTCAGAAGACCCCTTCTTAAACTCCAATATTTTGTAAATTGGTATGTTATCTCCTAGCGAGTGGTTAAAGGTGTTCATTTCCCAAACCATTGTACATCGCGCTAGACTCACTTATGATAATCGAGAATACGGTACGCCTGCGAGAGCAAATTTAGAAACACTCGACAG";
        let matches = search::<Dna>(query, text, 2);
        println!("Matches: {:?}", matches);
    }

    #[test]
    fn test_case3() {
        let query = b"AAAACCCAGT";
        let text = b"AAAACCAAGT";
        let matches = search::<Dna>(query, text, 2);
        assert!(matches.len() > 0);
    }

    #[test]
    fn no_extra_matches() {
        let edits = 6;
        let expected_idx = 277;
        let query = b"TAAGCAGAAGGGAGGTATAAAGTCTGTCAGCGGTGCTTAAG";
        let text = b"ACCGTAACCGCTTGGTACCATCCGGCCAGTCGCTCGTTGCGCCCCACTATCGGGATCGACGCGCAGTAATTAAACACCACCCACGCCACGAGGTAGAACGAGAGCGGGGGGCTAGCAAATAATAGTGAGAGTGCGTTCAAAGGGTCTTTCGTAACCTCAGCGGGCGGGTACGGGGGAAATATCGCACCAATTTTGGAGATGCGATTAGCTCAGCGTAACGCGAATTCCCTATAACTTGCCTAGTGTGTGTGAATGGACAATTCGTTTTACAGTTTCAAGGTAGCAGAAGGGCAGGATAAGTCTGTCGCGGTGCTTAAGGCTTTCCATCCATGTTGCCCCCTACATGAATCGGATCGCCAGCCAGAATATCACATGGTTCCAAAAGTTGCAAGCTTCCCCGTACCGCTACTTCACCTCACGCCAGAGGCCTATCGCCGCTCGGCCGTTCCGTTTTGGGGAAGAATCTGCCTGTTCTCGTCACAAGCTTCTTAGTCCTTCCACCATGGTGCTGTTACTCATGCCATCAAATATTCGAGCTCTTGCCTAGGGGGGTTATACCTGTGCGATAGATACACCCCCTATGACCGTAGGTAGAGAGCCTATTTTCAACGTGTCGATCGTTTAATGACACCAACTCCCGGTGTCGAGGTCCCCAAGTTTCGTAGATCTACTGAGCGGGGGAATATTTGACGGTAAGGCATCGCTTGTAGGATCGTATCGCGACGGTAGATACCCATAAGCGTTGCTAACCTGCCAATAACTGTCTCGCGATCCCAATTTAGCACAAGTCGGTGGCCTTGATAAGGCTAACCAGTTTCGCACCGCTTCCGTTCCATTTTACGATCTACCGCTCGGATGGATCCGAAATACCGAGGTAGTAATATCAACACGTACCCAATGTCC";
        let matches = super::search_bounded::<Dna>(query, text, edits);
        let m = matches
            .iter()
            .find(|m| (m.start.1 as usize).abs_diff(expected_idx) <= edits);
        assert!(m.is_some());
    }

    #[test]
    fn search_fuzz() {
        let mut query_lens = (10..20)
            .chain((0..10).map(|_| random_range(10..100)))
            .chain((0..10).map(|_| random_range(100..1000)))
            .collect::<Vec<_>>();

        let mut text_lens = (10..20)
            .chain((0..10).map(|_| random_range(10..100)))
            .chain((0..10).map(|_| random_range(100..1000)))
            .chain((0..10).map(|_| random_range(1000..10000)))
            .collect::<Vec<_>>();

        // let mut query_lens = [1000, 2000, 5000, 10_000].repeat(1);
        // let mut text_lens = [10_000, 100_000, 1_000_000].repeat(1);

        query_lens.sort();
        text_lens.sort();

        for q in query_lens {
            for t in text_lens.clone() {
                println!("q {q} t {t}");
                let query = (0..q)
                    .map(|_| b"ACGT"[random_range(0..4)])
                    .collect::<Vec<_>>();
                let mut text = (0..t)
                    .map(|_| b"ACGT"[random_range(0..4)])
                    .collect::<Vec<_>>();

                let edits = random_range(0..q / 3);
                let mut p = query.clone();
                for _ in 0..edits {
                    let tp = random_range(0..3);
                    match tp {
                        0 => {
                            // insert
                            let idx = random_range(0..=p.len());
                            p.insert(idx, b"ACGT"[random_range(0..4)]);
                        }
                        1 => {
                            // del
                            let idx = random_range(0..p.len());
                            p.remove(idx);
                        }
                        2 => {
                            // replace
                            let idx = random_range(0..p.len());
                            p[idx] = b"ACGT"[random_range(0..4)];
                        }
                        _ => panic!(),
                    }
                }

                fn show(x: &[u8]) -> &str {
                    str::from_utf8(x).unwrap()
                }
                eprintln!("");
                eprintln!("edits {edits}");
                eprintln!("query q={q} {}", show(&query));
                eprintln!("pattern {}", show(&p));

                if p.len() > text.len() {
                    continue;
                }

                let idx = random_range(0..=text.len() - p.len());
                eprintln!("text len {}", text.len());
                eprintln!("planted idx {idx}");
                let expected_idx = (idx + p.len()).saturating_sub(q);
                eprintln!("expected idx {expected_idx}");

                text.splice(idx..idx + p.len(), p);
                eprintln!("text {}", show(&text));

                // Plain
                let matches = super::search::<Dna>(&query, &text, edits);
                eprintln!("matches {matches:?}");
                let m = matches
                    .iter()
                    .find(|m| (m.start.1 as usize).abs_diff(expected_idx) <= edits);
                assert!(m.is_some());

                // Bounded
                let matches = super::search_bounded::<Dna>(&query, &text, edits);
                eprintln!("matches {matches:?}");
                let m = matches
                    .iter()
                    .find(|m| (m.start.1 as usize).abs_diff(expected_idx) <= edits);
                assert!(m.is_some());
            }
        }
    }
}
