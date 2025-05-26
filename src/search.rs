use crate::minima::prefix_min;
use crate::minima::prefix_min_k;
pub use crate::minima::{find_all_minima, find_local_minima};
use crate::profiles::Profile;
use crate::trace::{CostMatrix, fill, get_trace, simd_fill};
use crate::{LANES, S};
use crate::{
    bitpacking::compute_block_simd,
    delta_encoding::{V, VEncoding},
};
use pa_types::{Cigar, Cost, Pos};
use std::borrow::Cow;
use std::ops::Range;
use std::simd::cmp::SimdPartialOrd;

pub type Deltas = Vec<(Cost, V<u64>)>;
type Base = u64;
type VV = V<Base>;

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
struct LaneState<P: Profile> {
    decreasing: bool,
    text_slice: [u8; 64],
    text_profile: P::B,
    matches: Vec<(usize, Cost)>,
    chunk_offset: usize,
}

impl<P: Profile> LaneState<P> {
    fn new(text_profile: P::B, chunk_offset: usize) -> Self {
        Self {
            decreasing: false,
            text_slice: [0; 64],
            text_profile,
            matches: Vec::new(),
            chunk_offset,
        }
    }

    fn update_and_encode(&mut self, text: &[u8], i: usize, profiler: &P) {
        let start = self.chunk_offset * 64 + 64 * i;
        if start + 64 <= text.len() {
            self.text_slice.copy_from_slice(&text[start..start + 64]);
        } else {
            self.text_slice.fill(b'X');
            if start <= text.len() {
                let slice = &text[start..];
                self.text_slice[..slice.len()].copy_from_slice(slice);
            }
        }
        profiler.encode_ref(&self.text_slice, &mut self.text_profile);
    }
}

#[derive(Clone)]
pub struct Searcher<P: Profile, const RC: bool, const ALL_MINIMA: bool> {
    _phantom: std::marker::PhantomData<P>,
    cost_matrices: [CostMatrix; LANES],
    hp: Vec<S>,
    hm: Vec<S>,
    lanes: [LaneState<P>; LANES],
}

impl<P: Profile, const RC: bool, const ALL_MINIMA: bool> Searcher<P, RC, ALL_MINIMA> {
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
            cost_matrices: std::array::from_fn(|_| CostMatrix::default()),
            hp: Vec::new(),
            hm: Vec::new(),
            lanes: std::array::from_fn(|_| LaneState::new(P::alloc_out(), 0)),
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
        self.find_matches(query, text, k);
        self.process_matches(query, text, k)
    }

    fn find_matches(&mut self, query: &[u8], text: &[u8], k: usize) {
        self.search_positions_bounded(query, text, k as Cost)
    }

    #[inline(always)]
    fn check_lanes(
        &self,
        vp: &S,
        vm: &S,
        dist_to_start_of_lane: &S,
        k: Cost,
        j: usize,
    ) -> Option<usize> {
        for lane in 0..LANES {
            let v = V(vp.as_array()[lane], vm.as_array()[lane]);

            //Fixme: to go back to old impl. we could use prefix_min here again. Check speed difference
            let min_in_lane =
                prefix_min(v.0, v.1).0 as Cost + dist_to_start_of_lane.as_array()[lane] as Cost;
            if min_in_lane <= k {
                return Some(j + 4.max((k - min_in_lane) as usize));
            }
        }
        None
    }

    fn search_positions_bounded(&mut self, query: &[u8], text: &[u8], k: Cost) {
        let (profiler, query_profile) = P::encode_query(query);

        // Terminology:
        // - chunk: roughly 1/4th of the input text, with small overlaps.
        // - block: 64 bytes of text.
        // - lane: a u64 of a SIMD vec.

        // The query will match a pattern of length at most query.len() + k.
        // We round that up to a multiple of 64 to find the number of blocks overlap between chunks.
        let max_overlap_blocks = (query.len() + k as usize).next_multiple_of(64) / 64;
        let overlap_blocks = 0;

        // Total number of blocks to be processed, including overlaps.
        let text_blocks = text.len().div_ceil(64);
        let total_blocks = text_blocks + (LANES - 1) * overlap_blocks;
        let blocks_per_chunk = total_blocks.div_ceil(LANES);
        // Length of each of the four chunks.
        let chunk_offset = blocks_per_chunk.saturating_sub(overlap_blocks);

        // Update chunk offsets
        for lane in 0..LANES {
            self.lanes[lane].chunk_offset = lane * chunk_offset;
        }

        // Clear matches in each lane
        for lane in 0..LANES {
            self.lanes[lane].matches.clear();
        }

        // Up to where the previous column was computed.
        let mut prev_max_j = 0;
        // The max row where the right of the previous column was <=k
        let mut prev_end_last_below = 0;

        // SmallVec will now stack alloc 128 elements, and heap allocate beyond that.
        // Only little difference in practice, perhaps switch back to regular vecs
        self.hp.clear();
        self.hm.clear();
        self.hp.resize(query.len(), S::splat(1));
        self.hm.resize(query.len(), S::splat(0));

        'text_chunk: for i in 0..blocks_per_chunk + max_overlap_blocks {
            let mut vp = S::splat(0);
            let mut vm = S::splat(0);

            // Update text slices and profiles
            for lane in 0..LANES {
                self.lanes[lane].update_and_encode(text, i, &profiler);
            }

            let mut dist_to_start_of_lane = S::splat(0);
            let mut dist_to_end_of_lane = S::splat(0);
            let mut cur_end_last_below = 0;

            // Iterate over query chars
            for j in 0..query.len() {
                dist_to_start_of_lane += self.hp[j];
                dist_to_start_of_lane -= self.hm[j];

                let query_char = unsafe { query_profile.get_unchecked(j) };
                let eq: std::simd::Simd<u64, 4> = S::from(std::array::from_fn(|lane| {
                    P::eq(query_char, &self.lanes[lane].text_profile)
                }));

                compute_block_simd(&mut self.hp[j], &mut self.hm[j], &mut vp, &mut vm, eq);

                // For DNA, the distance between random/unrelated sequences is around q.len()/2.
                // Thus, for threshold k, we can expect random matches between seqs of length ~2*k.
                // To have some buffer, we start filtering at length 3*k.
                'check: {
                    dist_to_end_of_lane += self.hp[j];
                    dist_to_end_of_lane -= self.hm[j];

                    let cmp = dist_to_end_of_lane.simd_le(S::splat(k as u64));
                    let bitmask = cmp.to_bitmask(); // less assembly than .any
                    let end_leq_k = bitmask != 0;

                    cur_end_last_below = if end_leq_k { j } else { cur_end_last_below };

                    if j > prev_end_last_below {
                        if let Some(new_end) =
                            self.check_lanes(&vp, &vm, &dist_to_start_of_lane, k, j)
                        {
                            prev_end_last_below = new_end;
                            break 'check;
                        }

                        // No lanes had values <= k
                        for j2 in j + 1..=prev_max_j {
                            self.hp[j2] = S::splat(1);
                            self.hm[j2] = S::splat(0);
                        }
                        prev_end_last_below = cur_end_last_below.max(8);
                        prev_max_j = j;

                        if i >= blocks_per_chunk
                            && (64 * (i - blocks_per_chunk)).saturating_sub(j) > k as usize
                        {
                            break 'text_chunk;
                        }
                        continue 'text_chunk;
                    }
                }
            }

            // Save positions with cost <= k directly after processing each row
            for lane in 0..LANES {
                let v = <VV as VEncoding<Base>>::from(vp[lane], vm[lane]);

                let base_pos = self.lanes[lane].chunk_offset * 64 + 64 * i;
                let cost = dist_to_start_of_lane.as_array()[lane] as Cost;

                if ALL_MINIMA {
                    let matches_before = self.lanes[lane].matches.clone();
                    self.find_all_minima(v, cost, k, text.len(), base_pos, lane);
                    let matches_after = self.lanes[lane].matches.clone();
                } else {
                    self.find_local_minima(v, cost, k, text.len(), base_pos, lane);
                }
            }

            prev_end_last_below = cur_end_last_below.max(8);
            prev_max_j = query.len() - 1;
        }

        for j in 0..=prev_max_j {
            self.hp[j] = S::splat(1);
            self.hm[j] = S::splat(0);
        }
    }

    #[inline(always)]
    pub fn find_local_minima(
        &mut self,
        v: V<u64>,
        cur_cost: Cost,
        k: Cost,
        text_len: usize,
        base_pos: usize,
        lane: usize,
    ) {
        let (p, m) = v.pm();
        let changes = p | m;
        if changes == 0 {
            return;
        }

        let max_pos = (text_len.saturating_sub(base_pos)).min(64);
        let mut prev_cost = cur_cost;
        let mut cur_cost = cur_cost;
        let mut changes = changes;

        while changes != 0 {
            let pos = changes.trailing_zeros() as usize;

            if pos >= max_pos {
                break;
            }

            let delta = ((p >> pos) & 1) as Cost - ((m >> pos) & 1) as Cost;
            cur_cost += delta;

            // Check for local minimum
            let was_decreasing = self.lanes[lane].decreasing;
            let is_increasing = cur_cost > prev_cost;
            let is_now_decreasing = cur_cost < prev_cost;

            // Add match if we were decreasing and now increasing (local minimum)
            if was_decreasing && is_increasing && prev_cost <= k {
                self.lanes[lane].matches.push((base_pos + pos, prev_cost));
            }

            self.lanes[lane].decreasing = is_now_decreasing || (was_decreasing && !is_increasing);
            prev_cost = cur_cost;
            changes &= changes - 1;
        }

        // Check final position
        if self.lanes[lane].decreasing && cur_cost <= k && base_pos + 64 >= text_len {
            self.lanes[lane]
                .matches
                .push(((base_pos + 64).min(text_len), cur_cost));
        }
    }

    #[inline(always)]
    fn find_all_minima(
        &mut self,
        v: V<u64>,
        cur_cost: Cost,
        k: Cost,
        text_len: usize,
        base_pos: usize,
        lane: usize,
    ) {
        let (p, m) = v.pm();
        let mut cost = cur_cost;
        // All <=k end points
        for bit in 0..64 {
            let pos = base_pos + bit;
            if pos >= text_len {
                break;
            }
            // Check if this position is a match (cost <= k)
            if cost <= k {
                self.lanes[lane].matches.push((pos, cost));
            }

            // Update cost based on the P/M bit patterns
            let p_bit = ((p >> bit) & 1) as Cost;
            let m_bit = ((m >> bit) & 1) as Cost;
            cost += p_bit;
            cost -= m_bit;
        }
    }

    fn process_matches(&mut self, query: &[u8], text: &[u8], k: usize) -> Vec<Match> {
        let mut traces = Vec::new();
        let fill_len = query.len() + k;
        let mut last_processed_pos = 0;
        let mut text_slices = [[].as_slice(); LANES];
        let mut offsets = [0; LANES];
        let mut num_slices = 0;

        // We "pull" matches from each lane (left>right). As soon as we collect LANES slices
        // we use SIMD fill to compute the costs and traceback the path
        for lane in 0..LANES {
            for &(end_pos, _) in &self.lanes[lane].matches {
                // FIXME: We could be slightly more precise, and instead discard
                // up to the true end col of the overlap.
                if end_pos <= last_processed_pos {
                    continue;
                }

                let offset = end_pos.saturating_sub(fill_len);
                offsets[num_slices] = offset;
                text_slices[num_slices] = &text[offset..end_pos];
                num_slices += 1;

                // Process when we have a full chunk
                if num_slices == LANES {
                    // Fill cost matrices for all lanes
                    simd_fill::<P>(query, &text_slices[..num_slices], &mut self.cost_matrices);

                    // Process matches
                    for i in 0..num_slices {
                        let m = get_trace::<P>(
                            query,
                            offsets[i],
                            text_slices[i],
                            &self.cost_matrices[i],
                        );

                        assert!(
                            m.cost <= k as Cost,
                            "Match has cost {} > {}: {m:?}",
                            m.cost,
                            k
                        );
                        traces.push(m);
                    }
                    num_slices = 0;
                }
                last_processed_pos = end_pos;
            }
        }

        // We now have less < LANES text slices to process, we use SIMD fill if number of slices > 1
        // else we fall back to scalar fill to avoid SIMD overhead with "empty" lanes
        if num_slices > 0 {
            if num_slices > 1 {
                println!("num_slices: {num_slices} left");
                simd_fill::<P>(query, &text_slices[..num_slices], &mut self.cost_matrices);
                for i in 0..num_slices {
                    let m =
                        get_trace::<P>(query, offsets[i], text_slices[i], &self.cost_matrices[i]);
                    assert!(
                        m.cost <= k as Cost,
                        "Match has cost {} > {}: {m:?}",
                        m.cost,
                        k
                    );
                    traces.push(m);
                }
            } else {
                let mut cost_matrix = CostMatrix::default();
                fill::<P>(query, text_slices[0], &mut cost_matrix);
                let m = get_trace::<P>(query, offsets[0], text_slices[0], &cost_matrix);
                assert!(
                    m.cost <= k as Cost,
                    "Match has cost {} > {}: {m:?}",
                    m.cost,
                    k
                );
                traces.push(m);
            }
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
        let matches = Searcher::<Dna, true, false>::new().search(query, &text, 0);
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].start.1, 8);
        assert_eq!(matches[0].end.1, 8 + query.len() as i32);
        // Now disableing rc search should yield no matches
        let matches = Searcher::<Dna, false, false>::new().search(query, &text, 0);
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

                let idx = random_range(0..=text.len().saturating_sub(p.len()));
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

    #[test]
    fn test_fixed_matches() {
        let query = b"ATCGATCA";
        let mut text = vec![b'G'; 1000]; // Create a text of 1000 G's

        // Insert 5 matches at fixed positions
        let positions = [50, 150, 250, 350, 450, 800];
        let mut expected_matches = Vec::new();

        for &pos in &positions {
            // Insert the query at each position
            text.splice(pos..pos + query.len(), query.iter().copied());

            // Record expected match position
            let expected_idx = (pos + query.len()).saturating_sub(query.len());
            expected_matches.push(expected_idx);
        }

        // Test forward search
        let mut searcher = Searcher::<Dna, false, true>::new();
        let matches = searcher.search(query, &text, 1);

        // Verify all matches are found
        for expected_idx in expected_matches {
            let found = matches.iter().any(|m| m.start.1 as usize == expected_idx);
            assert!(found, "Expected match at {} not found", expected_idx);
        }

        for m in matches {
            println!("match: {:?}", m);
        }
    }
}
