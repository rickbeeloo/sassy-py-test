pub use crate::minima::{find_all_minima, find_local_minima};
use crate::profiles::Profile;
use crate::trace::{CostMatrix, get_trace, simd_fill};
use pa_types::{Cigar, Cost, Pos};
use std::borrow::Cow;
use std::simd::cmp::SimdPartialOrd;
use std::{array::from_fn, cell::RefCell};

use crate::minima::prefix_min;
use crate::{LANES, S};

use crate::{
    bitpacking::compute_block_simd,
    delta_encoding::{V, VEncoding},
};

pub type Deltas = Vec<(Cost, V<u64>)>;

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

#[derive(Clone, Copy)]
pub enum SearchMode {
    Local,
    All,
}

#[derive(Clone, Copy)]
pub enum RcMode {
    No,
    Yes,
}

pub trait SearchAble {
    /// The forward text
    fn text(&self) -> &[u8];
    /// Produce the reverse‐complement (or reverse) when requested
    fn rc_text(&self) -> Cow<[u8]>;
}

impl<T> SearchAble for T
where
    T: AsRef<[u8]>,
{
    fn text(&self) -> &[u8] {
        self.as_ref()
    }

    fn rc_text(&self) -> Cow<[u8]> {
        Cow::Owned(self.as_ref().iter().rev().copied().collect())
    }
}

pub struct StaticText<'a> {
    pub text: &'a [u8],
    pub rc: Vec<u8>,
}

impl<'a> StaticText<'a> {
    pub fn new(text: &'a [u8]) -> Self {
        let rc = text.iter().rev().copied().collect();
        StaticText { text, rc }
    }
}

impl<'a> SearchAble for StaticText<'a> {
    fn text(&self) -> &[u8] {
        self.text
    }
    fn rc_text(&self) -> Cow<[u8]> {
        // borrow stored, is free
        Cow::Borrowed(&self.rc)
    }
}

#[derive(Clone)]
pub struct Searcher<P: Profile, const RC: bool, const ALL_MINIMA: bool> {
    _phantom: std::marker::PhantomData<P>,
    deltas: Deltas,
}

impl<P: Profile, const RC: bool, const ALL_MINIMA: bool> Searcher<P, RC, ALL_MINIMA> {
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
            deltas: vec![],
        }
    }

    pub fn search<I: SearchAble>(&mut self, query: &[u8], input: &I, k: usize) -> Vec<Match> {
        let mut matches = self.search_internal(query, input.text(), k);
        if RC {
            let rc_matches = self.search_internal(&P::complement(query), &input.rc_text(), k);
            matches.extend(rc_matches.into_iter().map(|mut m| {
                m.strand = Strand::Rc;
                m
            }));
        }
        matches
    }

    fn search_internal(&mut self, query: &[u8], text: &[u8], k: usize) -> Vec<Match> {
        let matches = self.find_matches(query, text, k);
        if !matches.is_empty() {
            self.process_matches(matches, query, text, k)
        } else {
            vec![]
        }
    }

    fn find_matches(&mut self, query: &[u8], text: &[u8], k: usize) -> Vec<(usize, Cost)> {
        self.search_positions_bounded(query, text, k as Cost);

        if ALL_MINIMA {
            find_all_minima(query, &mut self.deltas, k as Cost, text.len())
        } else {
            find_local_minima(query, &mut self.deltas, k as Cost, text.len())
        }
    }

    fn search_positions_bounded(&mut self, query: &[u8], text: &[u8], k: Cost) {
        let (profiler, query_profile) = P::encode_query(query);

        // Terminology:
        // - chunk: roughly 1/4th of the input text, with small overlaps.
        // - block: 64 bytes of text.
        // - lane: a u64 of a SIMD vec.

        // The query will match a pattern of length at most query.len() + k.
        // We round that up to a multiple of 64 to find the number of blocks overlap between chunks.
        let overlap_blocks = (query.len() + k as usize).next_multiple_of(64) / 64;

        // Total number of blocks to be processed, including overlaps.
        let text_blocks = text.len().div_ceil(64);
        let total_blocks = text_blocks + (LANES - 1) * overlap_blocks;
        let blocks_per_chunk = total_blocks.div_ceil(LANES);
        // Length of each of the four chunks.
        let chunk_offset = blocks_per_chunk.saturating_sub(overlap_blocks);
        self.deltas.resize(text_blocks, Default::default());

        // We should also clear?
        self.deltas.fill(Default::default());

        type Base = u64;
        type VV = V<Base>;

        let mut hp = vec![S::splat(1); query.len()];
        let mut hm = vec![S::splat(0); query.len()];

        let mut text_profile: [P::B; LANES] = from_fn(|_| profiler.alloc_out());

        let mut text_chunks: [[u8; 64]; LANES] = [[0; 64]; LANES];

        // Up to where the previous column was computed.
        let mut prev_max_j = 0;
        // The max row where the right of the previous column was <=k
        let mut prev_end_last_below = usize::MAX;

        let first_check = (3 * k as usize).max(8);

        'text_chunk: for i in 0..blocks_per_chunk {
            // The alignment can start anywhere, so start with deltas of 0.
            let mut vp = S::splat(0);
            let mut vm = S::splat(0);

            // Collect the LANES slices of input text.
            // Out-of-bounds characters are replaced by 'X', which doesn't match anything.
            for lane in 0..LANES {
                let start = lane * chunk_offset * 64 + 64 * i;
                if start + 64 <= text.len() {
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

            let mut dist_to_start_of_lane = S::splat(0);
            let mut dist_to_end_of_lane = S::splat(0);

            let mut cur_end_last_below = 0;

            let mut next_check = first_check;

            // Iterate over query chars.
            for j in 0..query.len() {
                dist_to_start_of_lane += hp[j];
                dist_to_start_of_lane -= hm[j];

                let eq = from_fn(|lane| {
                    P::eq(
                        unsafe { query_profile.get_unchecked(j) },
                        &text_profile[lane],
                    )
                })
                .into();
                compute_block_simd(&mut hp[j], &mut hm[j], &mut vp, &mut vm, eq);

                // For DNA, the distance between random/unrelated sequences is around q.len()/2.
                // Thus, for threshold k, we can expect random matches between seqs of length ~2*k.
                // To have some buffer, we start filtering at length 3*k.
                'check: {
                    dist_to_end_of_lane += hp[j];
                    dist_to_end_of_lane -= hm[j];

                    if j >= first_check {
                        cur_end_last_below =
                            if (dist_to_end_of_lane.simd_le(S::splat(k as u64))).any() {
                                j
                            } else {
                                cur_end_last_below
                            };
                        if j == next_check {
                            next_check *= 2;

                            if j > prev_end_last_below {
                                // Check for each lane
                                for lane in 0..LANES {
                                    let v = V(vp.as_array()[lane], vm.as_array()[lane]);
                                    let min_in_lane = dist_to_start_of_lane.as_array()[lane]
                                        as Cost
                                        + prefix_min(v.0, v.1).0 as Cost;
                                    if min_in_lane <= k {
                                        // Go to the post-processing below.
                                        break 'check;
                                    }
                                }
                                // All lanes only have values > k. We set remaining horizontal deltas to +1.
                                for j2 in j + 1..=prev_max_j {
                                    hp[j2] = S::splat(1);
                                    hm[j2] = S::splat(0);
                                }
                                prev_end_last_below = cur_end_last_below;
                                for lane in 0..LANES {
                                    let idx = lane * chunk_offset + i;
                                    if idx < self.deltas.len() {
                                        self.deltas[idx] =
                                            (Cost::MAX, <VV as VEncoding<Base>>::from(u64::MAX, 0));
                                    }
                                }
                                prev_max_j = j;
                                continue 'text_chunk;
                            }
                        }
                    }
                }
            }
            for lane in 0..LANES {
                let idx = lane * chunk_offset + i;
                if idx < self.deltas.len() {
                    self.deltas[idx] = (
                        dist_to_start_of_lane.as_array()[lane] as _,
                        <VV as VEncoding<Base>>::from(vp[lane], vm[lane]),
                    );
                }
            }
            prev_end_last_below = cur_end_last_below;
            prev_max_j = query.len() - 1;
        }
    }

    fn process_matches(
        &self,
        matches: Vec<(usize, Cost)>,
        query: &[u8],
        text: &[u8],
        k: usize,
    ) -> Vec<Match> {
        let mut traces = Vec::with_capacity(matches.len());
        let fill_len = query.len() + k;

        thread_local! {
            static M: RefCell<[CostMatrix;LANES]> = RefCell::new(from_fn(|_| CostMatrix::default()));
        }

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
}

impl<P: Profile, const RC: bool, const ALL_MINIMA: bool> Default for Searcher<P, RC, ALL_MINIMA> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profiles::Dna;
    use rand::random_range;

    #[test]
    fn test_case1() {
        let query = b"AGATGTGTCC";
        let text = b"GAGAGATAACCGTGCGCTCACTGTTACAGTTTATGTGTCGAATTCTTTTAGGGCTGGTCACTGCCCATGCGTAGGAATGAATAGCGGTTGCGGATAACTAAGCAGTGCTTGTCGCTATTAAAGTTAGACCCCGCTCCCACCTTGCCACTCCTAGATGTCGACGAGATTCTCACCGCCAAGGGATCAGGCATATCATAGTAGCTGGCGCAGCCCGTCGTTTAAGGAGGCCCATATACTAATTCAAACGATGGGGTGGGCACATCCCCTAGAAGCACTTGTCCCTTGAGTCACGACTGATGCGTGGTCTCTCGCTAAATGTTCCGGCCTCTCGGACATTTAAAGGGTGGCATGTGACCATGGAGGATTAGTGAATGAGAGGTGTCCGCCTTTGTTCGCCAGAACTCTTATAGCGTAGGGGAGTGTACTCACCGCGAACCCGTATCAGCAATCTTGTCAGTGGCTCCTGACTCAAACATCGATGCGCTGCACATGGCCTTAGAATGAAGCAGCCCCTTTCTATTGTGGCCGGGCTGATTCTTTGTTTTGTTGAATGGTCGGGCCTGTCTGCCTTTTCCTAGTGTTGAAACTCCGAACCGCATGAACTGCGTTGCTAGCGAGCTATCACTGGACTGGCCGGGGGACGAAAGTTCGCGGGACCCACTACCCGCGCCCAGAAGACCACACTAGGGAGAAGGATTCTATCGGCATAGCCGTC";
        let matches = Searcher::<Dna, true, false>::new().search(query, &text, 2);
        println!("Matches: {:?}", matches);
    }

    #[test]
    fn no_extra_matches() {
        let edits = 6;
        let expected_idx = 277;
        let query = b"TAAGCAGAAGGGAGGTATAAAGTCTGTCAGCGGTGCTTAAG";
        let text = b"ACCGTAACCGCTTGGTACCATCCGGCCAGTCGCTCGTTGCGCCCCACTATCGGGATCGACGCGCAGTAATTAAACACCACCCACGCCACGAGGTAGAACGAGAGCGGGGGGCTAGCAAATAATAGTGAGAGTGCGTTCAAAGGGTCTTTCGTAACCTCAGCGGGCGGGTACGGGGGAAATATCGCACCAATTTTGGAGATGCGATTAGCTCAGCGTAACGCGAATTCCCTATAACTTGCCTAGTGTGTGTGAATGGACAATTCGTTTTACAGTTTCAAGGTAGCAGAAGGGCAGGATAAGTCTGTCGCGGTGCTTAAGGCTTTCCATCCATGTTGCCCCCTACATGAATCGGATCGCCAGCCAGAATATCACATGGTTCCAAAAGTTGCAAGCTTCCCCGTACCGCTACTTCACCTCACGCCAGAGGCCTATCGCCGCTCGGCCGTTCCGTTTTGGGGAAGAATCTGCCTGTTCTCGTCACAAGCTTCTTAGTCCTTCCACCATGGTGCTGTTACTCATGCCATCAAATATTCGAGCTCTTGCCTAGGGGGGTTATACCTGTGCGATAGATACACCCCCTATGACCGTAGGTAGAGAGCCTATTTTCAACGTGTCGATCGTTTAATGACACCAACTCCCGGTGTCGAGGTCCCCAAGTTTCGTAGATCTACTGAGCGGGGGAATATTTGACGGTAAGGCATCGCTTGTAGGATCGTATCGCGACGGTAGATACCCATAAGCGTTGCTAACCTGCCAATAACTGTCTCGCGATCCCAATTTAGCACAAGTCGGTGGCCTTGATAAGGCTAACCAGTTTCGCACCGCTTCCGTTCCATTTTACGATCTACCGCTCGGATGGATCCGAAATACCGAGGTAGTAATATCAACACGTACCCAATGTCC";
        let matches = Searcher::<Dna, false, false>::new().search(query, &text, edits);
        let m = matches
            .iter()
            .find(|m| (m.start.1 as usize).abs_diff(expected_idx) <= edits);
        assert!(m.is_some());
    }

    fn random_dna_string(len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[random_range(0..4)]).collect()
    }

    // Just for profiling
    // use std::hint::black_box;
    // #[test]
    // fn random_big_search() {
    //     let mut total_matches = 0;
    //     for i in 0..1000 {
    //         let query = random_dna_string(random_range(10..100));
    //         let text = random_dna_string(1_000_000);
    //         let matches = Search::<Dna, false, false>::new(&query, &text, 5).search();
    //         total_matches += matches.len();
    //     }
    //     println!("total matches: {total_matches}");
    // }

    #[test]
    fn test_fwd_rc_search() {
        let query = b"ATCGATCA";
        let rc = Dna::reverse_complement(query);
        let text = [b"GGGGGGGG".as_ref(), &rc, b"GGGGGGGG"].concat();
        let matches = Searcher::<Dna, true, false>::new().search(query, &text, 1);
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].start.1, 8);
        assert_eq!(matches[0].end.1, 8 + query.len() as i32);
        // Now disableing rc search should yield no matches
        let matches = Searcher::<Dna, false, false>::new().search(query, &text, 1);
        assert_eq!(matches.len(), 0);
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

        // Create single searcher for all tests to check proper resetting of internal states
        let mut searcher = Searcher::<Dna, false, false>::new();
        let mut rc_searcher = Searcher::<Dna, true, false>::new();

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

                // Just fwd
                let matches = searcher.search(&query, &text, edits);
                eprintln!("matches {matches:?}");
                let m = matches
                    .iter()
                    .find(|m| (m.start.1 as usize).abs_diff(expected_idx) <= edits);
                assert!(m.is_some());

                // Also rc search, should still find the same match
                let matches = rc_searcher.search(&query, &text, edits);

                eprintln!("matches {matches:?}");
                let m = matches
                    .iter()
                    .find(|m| (m.start.1 as usize).abs_diff(expected_idx) <= edits);
                assert!(m.is_some());
            }
        }
    }
}
