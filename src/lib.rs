#![feature(portable_simd, int_roundings)]
mod bitpacking;
mod delta_encoding;
pub mod profiles {

    mod ascii;
    mod dna;
    mod iupac;
    mod profile;

    pub use profile::Profile;

    pub use ascii::{Ascii, CaseInsensitiveAscii, CaseSensitiveAscii};
    pub use dna::Dna;
    pub use iupac::Iupac;
}

mod minima;
mod search;
mod trace;

pub use minima::{find_below_threshold, find_local_minima, find_local_minima_slow};
use pa_types::{Cigar, Cost, Pos};
pub use search::{search_positions, search_positions_bounded};

use profiles::Profile;
use trace::{get_trace, simd_fill};

#[derive(Debug, Clone)]
pub struct Match {
    pub start: Pos,
    pub end: Pos,
    pub cost: Cost,
    pub strand: Strand,
    pub cigar: Cigar,
}

#[derive(Debug, Clone)]
pub enum Strand {
    Fwd,
    Rc,
}

pub fn search<P: Profile>(query: &[u8], text: &[u8], k: usize) -> Vec<Match> {
    let mut deltas = vec![];
    search_positions::<P>(query, text, &mut deltas);

    let matches = find_local_minima_slow(query, &deltas, k as Cost);

    let mut traces = Vec::with_capacity(matches.len());

    let fill_len = query.len() + k;
    for matches in matches.chunks(4) {
        let mut text_slices = [[].as_slice(); 4];
        let mut offsets = [0; 4];
        for i in 0..matches.len() {
            let end_pos = matches[i].0;
            let offset = end_pos.saturating_sub(fill_len);
            offsets[i] = offset;
            text_slices[i] = &text[offset..end_pos];
        }
        // TODO: Reuse allocated costs.
        let costs = simd_fill::<P>(query, text_slices);

        for lane in 0..matches.len() {
            traces.push(get_trace::<P>(
                query,
                offsets[lane],
                text_slices[lane],
                &costs[lane],
            ));
        }
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
    fn search_fuzz() {
        let mut query_lens = (10..20)
            .chain((0..10).map(|_| random_range(10..100)))
            .collect::<Vec<_>>();
        query_lens.sort();
        let mut text_lens = (0..100).map(|_| random_range(0..1000)).collect::<Vec<_>>();
        text_lens.sort();
        for q in query_lens {
            for t in text_lens.clone() {
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
                eprintln!("query {}", show(&query));
                eprintln!("pattern {}", show(&p));

                if p.len() > text.len() {
                    continue;
                }

                let idx = random_range(0..=text.len() - p.len());
                eprintln!("text len {}", text.len());
                eprintln!("idx {idx}");
                text.splice(idx..idx + p.len(), p);
                eprintln!("text {}", show(&text));

                let matches = super::search::<Dna>(&query, &text, edits);
                eprintln!("matches {matches:?}");
                let m = matches
                    .iter()
                    .find(|m| (m.start.1 as usize).abs_diff(idx) <= edits);

                assert!(m.is_some());
            }
        }
    }
}
