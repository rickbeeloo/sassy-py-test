use crate::minima::prefix_min;
use crate::profiles::Profile;
use crate::trace::{CostMatrix, fill, get_trace, simd_fill};
use crate::{LANES, S};
use crate::{
    bitpacking::compute_block_simd,
    delta_encoding::{V, VEncoding},
};
use pa_types::{Cigar, CigarOp, Cost, Pos};
use std::borrow::Cow;
use std::simd::cmp::SimdPartialOrd;

pub type Deltas = Vec<(Cost, V<u64>)>;

#[derive(Debug, Clone)]
pub struct Match {
    pub start: Pos,
    pub end: Pos,
    pub cost: Cost,
    pub strand: Strand,
    pub cigar: Cigar,
}

impl Match {
    pub fn to_path(&self) -> Vec<Pos> {
        let mut pos = self.start;
        let mut path = vec![pos];
        for el in &self.cigar.ops {
            for _ in 0..el.cnt {
                pos += match el.op {
                    CigarOp::Match => Pos(1, 1),
                    CigarOp::Sub => Pos(1, 1),
                    CigarOp::Del => Pos(1, 0),
                    CigarOp::Ins => Pos(0, 1),
                };
                path.push(pos);
            }
        }
        path.pop();
        path
    }

    pub fn without_cigar(&self) -> Match {
        Match {
            start: self.start,
            end: self.end,
            cost: self.cost,
            strand: self.strand,
            cigar: Cigar::default(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Fwd,
    Rc,
}

pub trait SearchAble {
    /// The forward text
    fn text(&self) -> &[u8];
    /// Produce the reverse‐complement (or reverse) when requested
    fn rev_text(&self) -> Cow<[u8]>;
}

impl<T> SearchAble for T
where
    T: AsRef<[u8]>,
{
    fn text(&self) -> &[u8] {
        self.as_ref()
    }

    fn rev_text(&'_ self) -> Cow<'_, [u8]> {
        Cow::Owned(self.as_ref().iter().rev().copied().collect())
    }
}

pub struct StaticText<'a> {
    pub text: &'a [u8],
    pub rev: Vec<u8>,
}

impl<'a> StaticText<'a> {
    pub fn new(text: &'a [u8]) -> Self {
        let rev = text.iter().rev().copied().collect();
        StaticText { text, rev }
    }
}

impl<'a> SearchAble for StaticText<'a> {
    fn text(&self) -> &[u8] {
        self.text
    }
    fn rev_text(&'_ self) -> Cow<'_, [u8]> {
        // borrow stored, is free
        Cow::Borrowed(&self.rev)
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
            // Pad with N, so that costs at the end are diagonally preserved.
            self.text_slice.fill(b'N');
            if start <= text.len() {
                let slice = &text[start..];
                self.text_slice[..slice.len()].copy_from_slice(slice);
            }
        }
        profiler.encode_ref(&self.text_slice, &mut self.text_profile);
    }
}

#[derive(Clone)]
pub struct Searcher<P: Profile> {
    // Config
    rc: bool,
    /// overhang cost
    /// If set, must satisfy `0.0 <= alpha <= 1.0`
    alpha: Option<f32>,

    // Internal caches
    cost_matrices: [CostMatrix; LANES],
    hp: Vec<S>,
    hm: Vec<S>,
    lanes: [LaneState<P>; LANES],

    _phantom: std::marker::PhantomData<P>,
}

impl<P: Profile> Searcher<P> {
    // The number of rows (query chars) we *at least*
    // mainly to avoid branching
    pub const CHECK_AT_LEAST_ROWS: usize = 8;

    pub fn new_fwd() -> Self {
        Self::new(false, None)
    }

    pub fn new_rc() -> Self {
        Self::new(true, None)
    }

    fn _overhang_check(alpha: f32) {
        if !P::supports_overhang() {
            panic!(
                "Overhang is not supported for {:?}",
                std::any::type_name::<P>()
            );
        }
        if !(0.0..=1.0).contains(&alpha) {
            panic!("Alpha must be in range 0.0 <= alpha <= 1.0");
        }
    }

    pub fn new_fwd_with_overhang(alpha: f32) -> Self {
        Self::_overhang_check(alpha);
        Self::new(false, Some(alpha))
    }

    pub fn new_rc_with_overhang(alpha: f32) -> Self {
        Self::_overhang_check(alpha);
        Self::new(true, Some(alpha))
    }

    pub fn new(rc: bool, alpha: Option<f32>) -> Self {
        Self {
            rc,
            _phantom: std::marker::PhantomData,
            cost_matrices: std::array::from_fn(|_| CostMatrix::default()),
            hp: Vec::new(),
            hm: Vec::new(),
            lanes: std::array::from_fn(|_| LaneState::new(P::alloc_out(), 0)),
            alpha,
        }
    }

    /// Returns matches for *only local minima* end positions with score <=k.
    pub fn search<I: SearchAble>(&mut self, query: &[u8], input: &I, k: usize) -> Vec<Match> {
        self.search_handle_rc(
            query,
            input,
            k,
            false,
            None::<fn(&[u8], &[u8], Strand) -> bool>,
        )
    }

    /// Returns matches for *all* end positions with score <=k.
    pub fn search_all<I: SearchAble>(&mut self, query: &[u8], input: &I, k: usize) -> Vec<Match> {
        self.search_handle_rc(
            query,
            input,
            k,
            true,
            None::<fn(&[u8], &[u8], Strand) -> bool>,
        )
    }

    /// Returns matches for *all* end positions where end_filter_fn returns true
    pub fn search_with_fn<I: SearchAble>(
        &mut self,
        query: &[u8],
        input: &I,
        k: usize,
        all_minima: bool,
        filter_fn: impl Fn(&[u8], &[u8], Strand) -> bool,
    ) -> Vec<Match> {
        self.search_handle_rc(query, input, k, all_minima, Some(filter_fn))
    }

    fn search_handle_rc<I: SearchAble>(
        &mut self,
        query: &[u8],
        input: &I,
        k: usize,
        all_minima: bool,
        filter_fn: Option<impl Fn(&[u8], &[u8], Strand) -> bool>,
    ) -> Vec<Match> {
        let mut matches =
            self.search_one_strand(query, input.text(), k, all_minima, &filter_fn, Strand::Fwd);
        if self.rc {
            let rc_matches = self.search_one_strand(
                &P::complement(query),
                &input.rev_text(),
                k,
                all_minima,
                &filter_fn,
                Strand::Rc,
            );
            matches.extend(rc_matches.into_iter().map(|mut m| {
                m.strand = Strand::Rc;
                // Also adjust start and end positions to original text orientation
                let org_start = m.start.1;
                let org_end = m.end.1;
                m.start.1 = input.text().len() as i32 - org_end;
                m.end.1 = input.text().len() as i32 - org_start;
                m
            }));
        }
        matches
    }

    fn search_one_strand(
        &mut self,
        query: &[u8],
        text: &[u8],
        k: usize,
        all_minima: bool,
        filter_fn: &Option<impl Fn(&[u8], &[u8], Strand) -> bool>,
        strand: Strand,
    ) -> Vec<Match> {
        self.search_positions_bounded(query, text, k as Cost, all_minima);
        // If there is a filter fn, filter end positions based on function before processing matches
        if let Some(filter_fn) = filter_fn {
            self.lanes.iter_mut().for_each(|lane| {
                lane.matches.retain(|(end_pos, _)| {
                    let text_till_end = &text[..*end_pos];
                    filter_fn(query, text_till_end, strand)
                });
            });
        }
        self.process_matches(query, text, k as Cost)
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
            // Get the current cost state for this lane
            let v = V(vp.as_array()[lane], vm.as_array()[lane]);

            // Calculate the minimum possible cost in this lane
            // This is the best case scenario - if even this minimum is > k,
            // then no matches are possible in this lane
            let min_in_lane =
                prefix_min(v.0, v.1).0 as Cost + dist_to_start_of_lane.as_array()[lane] as Cost;

            if min_in_lane <= k {
                // Promising lane, we "estimate" how many rows more we need to check
                // as the difference between the minimum in the lane and the maximum edits
                // we can afford.
                let rows_needed = (k - min_in_lane) as usize;
                let new_end = j + Self::CHECK_AT_LEAST_ROWS.max(rows_needed);
                return Some(new_end);
            }
        }

        // No lanes are promising
        None
    }

    fn search_positions_bounded(&mut self, query: &[u8], text: &[u8], k: Cost, all_minima: bool) {
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
        let text_blocks = if self.alpha.is_some() {
            // When allowing overlaps, for simplicity we 'extend' the text a bit more with N.
            (text.len() + query.len()).div_ceil(64)
        } else {
            text.len().div_ceil(64)
        };
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
            self.lanes[lane].decreasing = true;
        }

        // State tracking for early termination optimization:
        // - prev_max_j: tracks the highest query row we've computed so far
        // - prev_end_last_below: tracks the highest row where any lane had cost <= k
        let mut prev_max_j = 0;
        let mut prev_end_last_below = 0;

        self.hp.clear();
        self.hm.clear();
        self.hp.resize(query.len(), S::splat(1));
        self.hm.resize(query.len(), S::splat(0));

        init_deltas_for_overshoot(&mut self.hp, self.alpha);

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

            // Iterate over query chars (rows in the DP matrix)
            for j in 0..query.len() {
                dist_to_start_of_lane += self.hp[j];
                dist_to_start_of_lane -= self.hm[j];

                let query_char = unsafe { query_profile.get_unchecked(j) };
                let eq: std::simd::Simd<u64, 4> = S::from(std::array::from_fn(|lane| {
                    P::eq(query_char, &self.lanes[lane].text_profile)
                }));

                compute_block_simd(&mut self.hp[j], &mut self.hm[j], &mut vp, &mut vm, eq);

                // Early termination check: If we've moved past the last row that had any
                // promising matches (cost <= k), we can potentially skip ahead or terminate
                'check: {
                    dist_to_end_of_lane += self.hp[j];
                    dist_to_end_of_lane -= self.hm[j];

                    // Check if any lane has cost <= k at the current row
                    let cmp = dist_to_end_of_lane.simd_le(S::splat(k as u64));
                    let bitmask = cmp.to_bitmask();
                    let end_leq_k = bitmask != 0;

                    // Track the highest row where we found any promising matches
                    cur_end_last_below = if end_leq_k { j } else { cur_end_last_below };

                    // Only do early termination checks if we've moved past the last promising row
                    if j > prev_end_last_below {
                        // Check if any lane has a minimum cost that could lead to matches <= k
                        if let Some(new_end) =
                            self.check_lanes(&vp, &vm, &dist_to_start_of_lane, k, j)
                        {
                            // Found a promising lane - update our tracking and continue
                            prev_end_last_below = new_end;
                            break 'check;
                        }

                        // No lanes have promising matches - we can skip ahead
                        self.reset_rows(j + 1, prev_max_j);
                        prev_end_last_below = cur_end_last_below.max(Self::CHECK_AT_LEAST_ROWS);
                        prev_max_j = j;

                        // Early termination: if we're in overlap region and too far from text end
                        if self.should_terminate_early(i, blocks_per_chunk, j, k) {
                            break 'text_chunk;
                        }
                        continue 'text_chunk;
                    }
                }
            }

            // Save positions with cost <= k directly after processing each row
            for lane in 0..LANES {
                let v = <V<u64> as VEncoding<u64>>::from(vp[lane], vm[lane]);
                let base_pos = self.lanes[lane].chunk_offset * 64 + 64 * i;
                let cost = dist_to_start_of_lane.as_array()[lane] as Cost;

                self.find_minima_with_overhang(
                    v,
                    cost,
                    k,
                    text.len(),
                    query.len(),
                    base_pos,
                    lane,
                    all_minima,
                );
            }

            prev_end_last_below = cur_end_last_below.max(Self::CHECK_AT_LEAST_ROWS);
            prev_max_j = query.len() - 1;
        }

        // Clean up any remaining rows that weren't reset
        self.reset_rows(0, prev_max_j);
    }

    /// Reset rows that are no longer needed for future computations
    #[inline(always)]
    fn reset_rows(&mut self, from_row: usize, to_row: usize) {
        for j2 in from_row..=to_row {
            self.hp[j2] = S::splat(1);
            self.hm[j2] = S::splat(0);
        }
    }

    /// Check if we should terminate early based on position and distance from text end
    #[inline(always)]
    fn should_terminate_early(
        &self,
        current_block: usize,
        blocks_per_chunk: usize,
        current_row: usize,
        k: Cost,
    ) -> bool {
        // Only consider early termination in the overlap region (after main chunks)
        if current_block < blocks_per_chunk {
            return false;
        }

        // Calculate how far we are from the end of the text
        let distance_from_end =
            (64 * (current_block - blocks_per_chunk)).saturating_sub(current_row);

        // If we're too far from the end, no matches are possible
        distance_from_end > k as usize
    }

    #[inline(always)]
    fn add_overshoot_cost(&self, cost: Cost, pos: usize, text_len: usize) -> Cost {
        let overshoot = pos.saturating_sub(text_len);
        let overshoot_cost = (self.alpha.unwrap_or(0.0) * overshoot as f32).floor() as Cost;
        cost + overshoot_cost
    }

    #[inline(always)]
    fn find_minima_with_overhang(
        &mut self,
        v: V<u64>,
        cur_cost: Cost,
        k: Cost,
        text_len: usize,
        query_len: usize,
        base_pos: usize,
        lane: usize,
        all_minima: bool,
    ) {
        let (p, m) = v.pm();
        let max_pos = if self.alpha.is_some() {
            text_len + query_len
        } else {
            text_len
        };

        if base_pos >= max_pos {
            return;
        }

        let mut cost = cur_cost;
        let mut prev_cost = self.add_overshoot_cost(cur_cost, base_pos, text_len);
        let mut prev_pos = base_pos;

        for bit in 1..=64 {
            cost += ((p >> (bit - 1)) & 1) as Cost;
            cost -= ((m >> (bit - 1)) & 1) as Cost;

            let pos: usize = base_pos + bit;
            if pos > max_pos {
                if !all_minima && self.lanes[lane].decreasing && prev_cost <= k {
                    self.lanes[lane].matches.push((prev_pos, prev_cost));
                }
                break;
            }

            let total_cost = self.add_overshoot_cost(cost, pos, text_len);

            if all_minima {
                if total_cost <= k {
                    self.lanes[lane].matches.push((pos, total_cost));
                }
            } else {
                // Local minima
                // Check how costs are changing
                let costs_are_equal = total_cost == prev_cost;
                let costs_are_increasing = total_cost > prev_cost;
                let costs_are_decreasing = total_cost < prev_cost;

                // Found a local minimum if we were decreasing and now costs are increasing
                if self.lanes[lane].decreasing && costs_are_increasing && prev_cost <= k {
                    self.lanes[lane].matches.push((prev_pos, prev_cost));
                }

                // Update decreasing state:
                // - If costs are decreasing, we're in a decreasing sequence
                // - If costs are equal, keep the previous state
                // - If costs are increasing, we're not decreasing
                self.lanes[lane].decreasing =
                    costs_are_decreasing || (self.lanes[lane].decreasing && costs_are_equal);
            }

            prev_cost = total_cost;
            prev_pos = pos;
        }
    }

    fn process_matches(&mut self, query: &[u8], text: &[u8], k: Cost) -> Vec<Match> {
        let mut traces = Vec::new();
        let fill_len = query.len() + k as usize;
        let mut last_processed_pos = 0;

        for m in &mut self.cost_matrices {
            m.alpha = self.alpha;
        }

        // Collect slices to process in batches
        let mut batch = MatchBatch::new();

        for lane in 0..LANES {
            for &(end_pos, cost) in &self.lanes[lane].matches {
                if end_pos <= last_processed_pos {
                    continue;
                }

                let offset = end_pos.saturating_sub(fill_len);
                let slice = &text[offset..end_pos.min(text.len())];

                batch.add(slice, offset, end_pos, cost);

                if batch.is_full() {
                    traces.extend(batch.process::<P>(
                        query,
                        fill_len,
                        &mut self.cost_matrices,
                        self.alpha,
                        k,
                    ));
                    batch.clear();
                }

                last_processed_pos = end_pos;
            }
        }

        if !batch.is_empty() {
            traces.extend(batch.process::<P>(
                query,
                fill_len,
                &mut self.cost_matrices,
                self.alpha,
                k,
            ));
        }

        traces
    }
}

struct MatchBatch<'a> {
    slices: [&'a [u8]; LANES],
    offsets: [usize; LANES],
    ends: [usize; LANES],
    expected_costs: [Cost; LANES],
    count: usize,
}

impl<'a> MatchBatch<'a> {
    fn new() -> Self {
        Self {
            slices: [b""; LANES],
            offsets: [0; LANES],
            ends: [0; LANES],
            expected_costs: [0; LANES],
            count: 0,
        }
    }

    fn add(&mut self, slice: &'a [u8], offset: usize, end: usize, cost: Cost) {
        self.slices[self.count] = slice;
        self.offsets[self.count] = offset;
        self.ends[self.count] = end;
        self.expected_costs[self.count] = cost;
        self.count += 1;
    }

    fn is_full(&self) -> bool {
        self.count == LANES
    }

    fn is_empty(&self) -> bool {
        self.count == 0
    }

    fn clear(&mut self) {
        self.count = 0;
        // Don't have to clear all data as add will keep track of count
        // which process uses to make sure it only uses filled data
    }

    fn process<P: Profile>(
        &self,
        query: &[u8],
        fill_len: usize,
        cost_matrices: &mut [CostMatrix; LANES],
        alpha: Option<f32>,
        k: Cost,
    ) -> Vec<Match> {
        if self.count > 1 {
            simd_fill::<P>(
                query,
                &self.slices[..self.count],
                fill_len,
                cost_matrices,
                alpha,
            );
        } else {
            fill::<P>(
                query,
                self.slices[0],
                fill_len,
                &mut cost_matrices[0],
                alpha,
            );
        }

        let mut results = Vec::with_capacity(self.count);

        for i in 0..self.count {
            let m = get_trace::<P>(
                query,
                self.offsets[i],
                self.ends[i],
                self.slices[i],
                &cost_matrices[i],
                alpha,
            );

            // Check if get_trace cost is same as expected end position cost
            assert!(
                m.cost <= self.expected_costs[i],
                "Match has unexpected cost {} > {}: {m:?}",
                m.cost,
                self.expected_costs[i],
            );

            // Make sure it's also <=k
            assert!(
                m.cost <= k,
                "Match exceeds k after traceback: m.cost={}, k={}",
                m.cost,
                k,
            );

            results.push(m);
        }

        results
    }
}

/// Assumes hp and hm are already the right size, hm=0 and hp=1.
/// Then sets hp according to the given alpha, if needed.
pub(crate) fn init_deltas_for_overshoot_scalar(h: &mut [(u64, u64)], alpha: Option<f32>) {
    if let Some(alpha) = alpha {
        for i in 0..h.len() {
            // Alternate 0 and 1 costs at very left of the matrix.
            // (Note: not at start of later chunks.)
            // FIXME: floor, round, or ceil?
            h[i].0 =
                (((i + 1) as f32) * alpha).floor() as u64 - ((i as f32) * alpha).floor() as u64;
        }
    }
}

/// Assumes hp and hm are already the right size, hm=0 and hp=1.
/// Then sets hp according to the given alpha, if needed.
pub(crate) fn init_deltas_for_overshoot(hp: &mut [S], alpha: Option<f32>) {
    if let Some(alpha) = alpha {
        for i in 0..hp.len() {
            // Alternate 0 and 1 costs at very left of the matrix.
            // (Note: not at start of later chunks.)
            // FIXME: floor, round, or ceil?
            hp[i].as_mut_array()[0] =
                (((i + 1) as f32) * alpha).floor() as u64 - ((i as f32) * alpha).floor() as u64;
        }
    }
}

/// Assumes hp and hm are already the right size, hm=0 and hp=1.
/// Then sets hp according to the given alpha, if needed.
pub(crate) fn init_deltas_for_overshoot_all_lanes(hp: &mut [S], alpha: Option<f32>) {
    if let Some(alpha) = alpha {
        for i in 0..hp.len() {
            // Alternate 0 and 1 costs at very left of the matrix.
            // (Note: not at start of later chunks.)
            // FIXME: floor, round, or ceil?
            let bit =
                (((i + 1) as f32) * alpha).floor() as u64 - ((i as f32) * alpha).floor() as u64;
            hp[i].as_mut_array().fill(bit);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profiles::{Dna, Iupac};
    use rand::random_range;

    #[test]
    fn overhang_test() {
        let query = b"CTTAAGCACTACCGGCTAAT";
        let text = b"AGTCGTCCTTTGCGAGCTCGGACATCTCCAGGCGAACCTGCAAGTTTTAATGTTCCCACAGTCCCTCATATGTTCTGAATTTCGTGATGTTTGTTTACCG";
        let mut s = Searcher::<Iupac>::new_fwd_with_overhang(0.0);
        let _matches = s.search_all(query, text, 100);
    }

    #[test]
    #[should_panic()]
    fn overhang_test_panic_for_dna() {
        let mut searcher = Searcher::<Dna>::new_fwd_with_overhang(0.0);
    }

    #[test]
    fn overshoot() {
        let query = b"CCCTTTCCCGGG";
        let text = b"AAAAAAAAACCCTTT";
        let mut s = Searcher::<Iupac>::new_fwd();
        s.alpha = Some(0.5);
        s.search_positions_bounded(query, text, 10, true);
        for l in s.lanes {
            println!("Matches: {:?}", l.matches);
        }
    }

    #[test]
    fn overshoot_test_prefix_trace() {
        let query = b"CCCTTTCCCGGG";
        let text = b"AAAAAAAAACCCTTT";
        let mut s = Searcher::<Iupac>::new_fwd();
        s.alpha = Some(0.5);
        s.search_all(query, text, 10);
        // First not error
    }

    #[test]
    fn overshoot_simple_prefix() {
        /*
        AAAAGGGG
            ||||
            GGGGTTTTTTTTTTTTTTTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
        01234567---
            0123---
        */
        let prefix = "AAAAGGGG";
        let text = "GGGGTTTTTTTTTTTTTTTT";
        let mut s = Searcher::<Iupac>::new_fwd();
        s.alpha = Some(0.5);
        s.search_positions_bounded(prefix.as_bytes(), text.as_bytes(), 2, true);
        let expected_idx = 3;
        let expected_edits = 2 as Cost;
        let m = s.lanes[0]
            .matches
            .iter()
            .find(|m| m.0 == expected_idx && m.1 <= expected_edits);
        assert!(m.is_some());
    }

    #[test]
    fn overshoot_simple_suffix() {
        /*
                            GGGGAAAA
                            ||||
            TTTTTTTTTTTTTTTTGGGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
            0123456789-123456789-12345
                                   ^23,24,25, etc
        */
        let prefix = "GGGGAAAA";
        let text = "TTTTTTTTTTTTTTTTGGGG";
        let mut s = Searcher::<Iupac>::new_fwd();
        s.alpha = Some(0.5);
        s.search_positions_bounded(prefix.as_bytes(), text.as_bytes(), 2, true);
        let expected_idx = 24;
        let expected_edits = 2 as Cost;
        let m = s.lanes[0]
            .matches
            .iter()
            .find(|m| m.0 == expected_idx && m.1 <= expected_edits);
        assert!(m.is_some());
    }

    #[test]
    fn overshoot_simple_suffix_local_minima() {
        /*
                            GGGGAAAA
                            ||||
            TTTTTTTTTTTTTTTTGGGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
            0123456789-123456789-12345
                                   ^23,24,25, etc
        */
        let prefix = b"GGGGAAAA";
        let text = b"TTTTTTTTTTTTTTTTGGGG";
        let mut s = Searcher::<Iupac>::new_fwd();
        s.alpha = Some(0.5);
        let matches = s.search(prefix, text, 4);
        let expected_end_pos = Pos(3, 20);
        let expected_edits = 2 as Cost;
        for m in matches.iter() {
            println!("Match: {:?}", m.without_cigar());
        }
        let m = matches
            .iter()
            .find(|m| m.end == expected_end_pos && m.cost == expected_edits);
        assert!(m.is_some());
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn overshoot_test_prefix_and_suffix() {
        /*
              AAAAGGGG                 AAAAGGGG
                  ||||                 ||||
                  GGGGGAAAAA     GGGGGAAAAANNNN
                  0123456789     0123456789-123
                     ^ 3                      ^ 13
        */
        let contained = "AAAAGGGG";
        let text = "GGGGGAAAAA";
        let mut s = Searcher::<Iupac>::new_fwd();
        s.alpha = Some(0.5);
        s.search_positions_bounded(contained.as_bytes(), text.as_bytes(), 2, true);
        let expected_indices = [3, 13];
        let expected_edits = [2, 2];
        let mut found = [false, false];
        for m in s.lanes[0].matches.iter() {
            for j in 0..expected_indices.len() {
                if m.0 == expected_indices[j] && m.1 == expected_edits[j] {
                    found[j] = true;
                }
            }
        }
        assert!(found[0]);
        assert!(found[1]);
    }

    #[test]
    fn test_case1() {
        let query = b"AGATGTGTCC";
        let text = b"GAGAGATAACCGTGCGCTCACTGTTACAGTTTATGTGTCGAATTCTTTTAGGGCTGGTCACTGCCCATGCGTAGGAATGAATAGCGGTTGCGGATAACTAAGCAGTGCTTGTCGCTATTAAAGTTAGACCCCGCTCCCACCTTGCCACTCCTAGATGTCGACGAGATTCTCACCGCCAAGGGATCAGGCATATCATAGTAGCTGGCGCAGCCCGTCGTTTAAGGAGGCCCATATACTAATTCAAACGATGGGGTGGGCACATCCCCTAGAAGCACTTGTCCCTTGAGTCACGACTGATGCGTGGTCTCTCGCTAAATGTTCCGGCCTCTCGGACATTTAAAGGGTGGCATGTGACCATGGAGGATTAGTGAATGAGAGGTGTCCGCCTTTGTTCGCCAGAACTCTTATAGCGTAGGGGAGTGTACTCACCGCGAACCCGTATCAGCAATCTTGTCAGTGGCTCCTGACTCAAACATCGATGCGCTGCACATGGCCTTAGAATGAAGCAGCCCCTTTCTATTGTGGCCGGGCTGATTCTTTGTTTTGTTGAATGGTCGGGCCTGTCTGCCTTTTCCTAGTGTTGAAACTCCGAACCGCATGAACTGCGTTGCTAGCGAGCTATCACTGGACTGGCCGGGGGACGAAAGTTCGCGGGACCCACTACCCGCGCCCAGAAGACCACACTAGGGAGAAGGATTCTATCGGCATAGCCGTC";
        let matches = Searcher::<Dna>::new_rc().search(query, &text, 2);
        println!("Matches: {:?}", matches);
    }

    #[test]
    fn no_extra_matches() {
        let edits = 6;
        let expected_idx = 277;
        let query = b"TAAGCAGAAGGGAGGTATAAAGTCTGTCAGCGGTGCTTAAG";
        let text = b"ACCGTAACCGCTTGGTACCATCCGGCCAGTCGCTCGTTGCGCCCCACTATCGGGATCGACGCGCAGTAATTAAACACCACCCACGCCACGAGGTAGAACGAGAGCGGGGGGCTAGCAAATAATAGTGAGAGTGCGTTCAAAGGGTCTTTCGTAACCTCAGCGGGCGGGTACGGGGGAAATATCGCACCAATTTTGGAGATGCGATTAGCTCAGCGTAACGCGAATTCCCTATAACTTGCCTAGTGTGTGTGAATGGACAATTCGTTTTACAGTTTCAAGGTAGCAGAAGGGCAGGATAAGTCTGTCGCGGTGCTTAAGGCTTTCCATCCATGTTGCCCCCTACATGAATCGGATCGCCAGCCAGAATATCACATGGTTCCAAAAGTTGCAAGCTTCCCCGTACCGCTACTTCACCTCACGCCAGAGGCCTATCGCCGCTCGGCCGTTCCGTTTTGGGGAAGAATCTGCCTGTTCTCGTCACAAGCTTCTTAGTCCTTCCACCATGGTGCTGTTACTCATGCCATCAAATATTCGAGCTCTTGCCTAGGGGGGTTATACCTGTGCGATAGATACACCCCCTATGACCGTAGGTAGAGAGCCTATTTTCAACGTGTCGATCGTTTAATGACACCAACTCCCGGTGTCGAGGTCCCCAAGTTTCGTAGATCTACTGAGCGGGGGAATATTTGACGGTAAGGCATCGCTTGTAGGATCGTATCGCGACGGTAGATACCCATAAGCGTTGCTAACCTGCCAATAACTGTCTCGCGATCCCAATTTAGCACAAGTCGGTGGCCTTGATAAGGCTAACCAGTTTCGCACCGCTTCCGTTCCATTTTACGATCTACCGCTCGGATGGATCCGAAATACCGAGGTAGTAATATCAACACGTACCCAATGTCC";
        let matches = Searcher::<Dna>::new_fwd().search(query, &text, edits);
        let m = matches
            .iter()
            .find(|m| (m.start.1 as usize).abs_diff(expected_idx) <= edits);
        assert!(m.is_some());
    }

    fn random_dna_string(len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[random_range(0..4)]).collect()
    }

    // Just for profiling
    #[test]
    #[ignore = "for profiling only"]
    fn random_big_search() {
        let mut total_matches = 0;
        for _ in 0..1000 {
            let query = random_dna_string(random_range(10..100));
            let text = random_dna_string(1_000_000);
            let matches = Searcher::<Dna>::new_fwd().search(&query, &text, 5);
            total_matches += matches.len();
        }
        println!("total matches: {total_matches}");
    }

    #[test]
    fn test_fwd_rc_search() {
        let query = b"ATCGATCA";
        let rc = Dna::reverse_complement(query);
        let text = [b"GGGGGGGG".as_ref(), &rc, b"GGGGGGGG"].concat();
        let matches = Searcher::<Dna>::new_rc().search(query, &text, 0);
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].start.1, 8);
        assert_eq!(matches[0].end.1, 8 + query.len() as i32);
        // Now disableing rc search should yield no matches
        let matches = Searcher::<Dna>::new_fwd().search(query, &text, 0);
        assert_eq!(matches.len(), 0);
    }

    #[test]
    fn test_filter_fn_simple() {
        let query = b"ATCGATCA";
        let mut text = vec![b'G'; 100];

        // Insert match once before 10 and once after 10
        text.splice(10..10, query.iter().copied());
        text.splice(50..50, query.iter().copied());
        let end_filter = |q: &[u8], text: &[u8], _strand: Strand| text.len() > 10 + q.len();
        let matches = Searcher::<Dna>::new_fwd().search_with_fn(query, &text, 0, false, end_filter);
        assert_eq!(matches.len(), 1); // First match *ending* at 10 should be discarded
        assert_eq!(matches[0].start.1, 50);

        // Sanity check, run the same without filter
        let matches = Searcher::<Dna>::new_fwd().search(query, &text, 0);
        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0].start.1, 10);
        assert_eq!(matches[1].start.1, 50);
    }

    fn rc(text: &[u8]) -> Vec<u8> {
        let mut rc = text.to_vec();
        rc.reverse();
        rc.iter_mut().for_each(|c| {
            *c = match *c {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => *c,
            }
        });
        rc
    }

    fn complement(text: &[u8]) -> Vec<u8> {
        let mut complement = text.to_vec();
        complement.iter_mut().for_each(|c| {
            *c = match *c {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => *c,
            }
        });
        complement
    }

    #[test]
    fn test_filter_fn_rc() {
        let query_fwd = b"ATCGATCA";
        let query_rc = rc(query_fwd);
        let mut text = vec![b'G'; 100];

        // Insert match once before 10 and once after 10
        text.splice(10..10, query_fwd.iter().copied()); // FWD
        text.splice(50..50, query_rc.iter().copied()); // RC

        let end_filter = |q: &[u8], text: &[u8], strand: Strand| match strand {
            Strand::Fwd => text[text.len() - q.len()..] == *query_fwd,
            Strand::Rc => {
                complement(&text[text.len() - q.len()..]) == *query_fwd // NOTE complement call 
            }
        };

        let matches =
            Searcher::<Dna>::new_rc().search_with_fn(query_fwd, &text, 0, false, end_filter);
        assert_eq!(matches.len(), 2); // Both matches should be found
        assert_eq!(matches[0].start.1, 10);
        assert_eq!(matches[1].start.1, 50);
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
        let mut searcher = Searcher::<Dna>::new_fwd();
        let mut rc_searcher = Searcher::<Dna>::new_rc();

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

    use rand::Rng;
    use serde::{Deserialize, Serialize};

    #[derive(Copy, Clone, Debug, PartialEq, Deserialize, Serialize)]
    #[serde(rename_all = "lowercase")]
    pub enum Alphabet {
        Dna,
        Iupac,
        Ascii,
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
        let mut searcher = Searcher::<Dna>::new_fwd();
        let matches = searcher.search_all(query, &text, 1);

        // Verify all matches are found
        for expected_idx in expected_matches {
            let found = matches.iter().any(|m| m.start.1 as usize == expected_idx);
            assert!(found, "Expected match at {} not found", expected_idx);
        }

        for m in matches {
            println!("match: {:?}", m);
        }
    }

    #[test]
    fn overhang_trace_fuzz() {
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};
        use std::iter::repeat_with;

        let mut rng = StdRng::seed_from_u64(42);
        let mut searcher = Searcher::<Iupac>::new_fwd();
        searcher.alpha = Some(0.5);

        fn rand_dna_w_seed(len: usize, rng: &mut StdRng) -> Vec<u8> {
            repeat_with(|| {
                let n = rng.random_range(0..4);
                match n {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                }
            })
            .take(len)
            .collect()
        }

        let mut skipped = 0;
        let iter = 1000;

        for _ in 0..iter {
            // Random query (short for testing)
            let query_len = rng.random_range(1..=100);
            let query = rand_dna_w_seed(query_len, &mut rng);

            // Random text (short for testing)
            let text_len = rng.random_range(1..=1000);
            let mut text = rand_dna_w_seed(text_len, &mut rng);

            // generate overlap at the prefix and suffix of the text
            let prefix_overlap = rng.random_range(1..=query_len.min(text_len));
            let suffix_overlap = rng.random_range(1..=query_len.min(text_len));

            // Ensure there's at least one character spacing between prefix and suffix
            if prefix_overlap + suffix_overlap >= text_len {
                skipped += 1;
                continue;
            }

            text.splice(
                0..prefix_overlap,
                query[query_len - prefix_overlap..].iter().copied(),
            );

            let expected_prefix_cost = ((query.len() as f32 - prefix_overlap as f32) * 0.5).floor();
            let expected_prefix_end_pos = prefix_overlap;

            // suffix overlap means we insert "suffix_overlap" start of the query at the end of the text
            text.splice(
                text_len - suffix_overlap..text_len,
                query[..suffix_overlap].iter().copied(),
            );
            let expected_suffix_cost = ((query.len() as f32 - suffix_overlap as f32) * 0.5).floor();
            let expected_suffix_end_pos = text_len;

            eprintln!("Q: {}", String::from_utf8_lossy(&query));
            eprintln!("T: {}", String::from_utf8_lossy(&text));
            eprintln!("Query len {query_len}");
            eprintln!("Text len {text_len}");
            eprintln!("[prefix] overlap {prefix_overlap}");
            eprintln!("[suffix] overlap {suffix_overlap}");
            eprintln!("[prefix] expected_cost {expected_prefix_cost}");
            eprintln!("[prefix] expected_end_pos {expected_prefix_end_pos}");
            eprintln!("[suffix] expected_cost {expected_suffix_cost}");
            eprintln!("[suffix] expected_end_pos {expected_suffix_end_pos}");
            eprintln!("--------------------------------");

            // Allow all k for now but later should be k
            let matches = searcher.search_all(&query, &text, query_len);
            // Check if matches are found with expected cost at expected positions
            let mut found = [false, false];
            let expected_locs = [expected_prefix_end_pos, expected_suffix_end_pos];
            let expected_costs = [expected_prefix_cost, expected_suffix_cost];
            for m in matches {
                println!("m: {:?} {:?} {}", m.start, m.end, m.cost);
                for i in 0..expected_locs.len() {
                    if m.end.1 == expected_locs[i] as i32 && m.cost == expected_costs[i] as Cost {
                        found[i] = true;
                    }
                }
            }
            assert!(found[0], "Expected prefix overlap not found");
            assert!(found[1], "Expected suffix overlap not found");
        }
        eprintln!("Passed: {} (skipped: {})", iter - skipped, skipped);
    }

    use pa_types::*;

    #[test]
    fn test_query_trace_path_0_edits() {
        /*
           Q:     ATGC
           T: GGGGATGCGGG
              0123456789*
        */
        let query = b"ATGC";
        let text = b"GGGGATGCGGG";
        let mut searcher = Searcher::<Dna>::new_fwd();
        let matches = searcher.search(query, text, 0);
        let path = matches[0].to_path();
        assert_eq!(path, vec![Pos(0, 4), Pos(1, 5), Pos(2, 6), Pos(3, 7)]);
        // Ends are exclusive
        assert_eq!(matches[0].end, *path.last().unwrap() + Pos(1, 1));
    }

    #[test]
    fn test_query_trace_path_1_edits() {
        let query = b"ATGC";
        let text = b"GGGGATTGCGGG";
        let mut searcher = Searcher::<Dna>::new_fwd();
        let matches = searcher.search(query, text, 1);
        let path = matches[0].to_path();
        assert_eq!(
            path,
            vec![Pos(0, 4), Pos(1, 5), Pos(1, 6), Pos(2, 7), Pos(3, 8)]
        );
        // Ends are exclusive
        assert_eq!(matches[0].end, *path.last().unwrap() + Pos(1, 1));
    }

    #[test]
    fn test_query_trace_path_with_overhang_prefix() {
        let query = b"ATCGATCG";
        let text = b"ATCGGGGGGGGGG"; // half of query removed at start
        let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(0.5);
        searcher.alpha = Some(0.5);
        let matches = searcher.search(query, text, 2);
        let path = matches[0].to_path();
        // This "skips" the first 4 character of the query as they are in the overhang
        assert_eq!(path, vec![Pos(4, 0), Pos(5, 1), Pos(6, 2), Pos(7, 3)]);
        // Ends are exclusive
        assert_eq!(matches[0].end, *path.last().unwrap() + Pos(1, 1));
    }

    #[test]
    fn test_query_trace_path_with_overhang_suffix() {
        let query = b"ATCGATCG";
        let text = b"GGGGGGGATCG"; // half of query removed at end
        let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(0.5);
        searcher.alpha = Some(0.5);
        let matches = searcher.search(query, text, 2);
        println!("matches: {:?}", matches);
        let path = matches[0].to_path();
        assert_eq!(path, vec![Pos(0, 7), Pos(1, 8), Pos(2, 9), Pos(3, 10)]);
        // Ends are exclusive
        assert_eq!(matches[0].end, *path.last().unwrap() + Pos(1, 1));
    }

    #[test]
    fn test_random_queries_60_range() {
        use rand::Rng;
        let mut rng = rand::rng();
        let mut i = 0;

        for _ in 0..10000 {
            let mut searcher = Searcher::<Iupac>::new_rc_with_overhang(0.4);

            // Generate random query of length 126
            let query: Vec<u8> = (0..126).map(|_| b"ACGT"[rng.random_range(0..4)]).collect();

            // Generate random text length between 62-90
            let text_len = rng.random_range(62..91);
            let text: Vec<u8> = (0..text_len)
                .map(|_| b"ACGT"[rng.random_range(0..4)])
                .collect();

            // Use k as half of query length
            let k = query.len() / 2;

            let matches = searcher.search(&query, &text, k);

            // Print every 1000 iterations
            i += 1;
            println!(
                "Iteration {}: Q={}, T={}, k={}, matches={}\nQuery: {}\nText: {}",
                i,
                query.len(),
                text.len(),
                k,
                matches.len(),
                String::from_utf8_lossy(&query),
                String::from_utf8_lossy(&text)
            );

            // Verify matches
            for m in matches {
                assert!(
                    m.cost <= k as Cost,
                    "Match has cost {} > {}: {m:?}\nQuery: {}\nText: {}\n",
                    m.cost,
                    k,
                    String::from_utf8_lossy(&query),
                    String::from_utf8_lossy(&text)
                );
            }
        }
    }

    #[test]
    fn test_case3() {
        let query = b"GTCTTTCATTCTCTCATCATAATCTCTAATACGACACATTGTACATCTGCTTGCGAGCCGGTGTAGCGCCGTCCTGTTATTTCAAGGCTATAATTACGAATTCAATTCCTCCTCTTCCAAAACACG";
        let text = b"AGTGATATCTCAAGGGGCCCTATTGGAAGGAAAGCCGCGATGGGTTCAACGTCAAGTGGATCATTCGATATTCATTAGCCCAACAGAAAC";
        let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(0.4);
        let matches = searcher.search(query, &text, 63);
        println!("Matches: {:?}", matches);
    }

    #[test]
    fn test_case4() {
        let expected_end_pos = 1;
        let expected_cost = 1;
        let query = b"ATC";
        let text = b"CGGGGGG";
        let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(0.5);
        let matches = searcher.search(query, &text, query.len());
        let all_matches = searcher.search_all(query, &text, query.len());

        println!("[ALL MATCHES]");
        let mut all_found = false;
        for m in all_matches {
            println!("\t{}-{} c: {}", m.start, m.end, m.cost);
            if m.end.1 == expected_end_pos && m.cost == expected_cost {
                all_found = true;
                break;
            }
        }

        println!("[LOCAL MATCHES]");
        let mut local_found = false;
        for m in matches {
            println!("\t{}-{} c: {}", m.start, m.end, m.cost);
            if m.end.1 == expected_end_pos && m.cost == expected_cost {
                local_found = true;
                break;
            }
        }

        assert!(
            all_found,
            "No ALL match found ending at {expected_end_pos} with cost {expected_cost}"
        );
        assert!(
            local_found,
            "No LOCAL match found ending at {expected_end_pos} with cost {expected_cost}"
        );
    }

    #[test]
    fn test_match_exact_at_end() {
        let query = b"ATAC".to_vec();
        let text = b"CCCCCCATAC";
        let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(0.5);
        let matches = searcher.search(&query, &text, 0);
        println!("Matches: {:?}", matches);
        let all_matches = searcher.search_all(&query, &text, 0);
        println!("All matches: {:?}", all_matches)
    }

    #[test]
    fn fwd_rc_test_simple() {
        let query = b"ATCATGCTAGC".to_vec();
        let text = b"GGGGGGGGGGATCATGCTAGCGGGGGGGGGGG".to_vec();
        let rc = Iupac::reverse_complement(&query);

        let mut searcher = Searcher::<Iupac>::new_rc_with_overhang(0.5);
        let fwd_matches = searcher.search(&query, &text, 0);
        let rc_matches = searcher.search(&rc, &text, 0);

        assert_eq!(
            fwd_matches.len(),
            rc_matches.len(),
            "Simple test: Forward and RC searches should find the same number of matches"
        );

        for fwd_match in fwd_matches.iter() {
            let matching_rc = rc_matches.iter().find(|rc_match| {
                rc_match.start == fwd_match.start
                    && rc_match.end == fwd_match.end
                    && rc_match.cost == fwd_match.cost
            });
            assert!(
                matching_rc.is_some(),
                "No matching RC match found for forward match: {:?}",
                fwd_match.without_cigar()
            );
        }
    }

    #[test]
    fn fwd_rc_test() {
        let fwd = b"TGAAGCGGCGCACGAAAAACGCGAAAGCGTTTCACGATAAATGCGAAAACNNNNNNNNNNNNNNNNNNNNNNNNGGTTAAACACCCAAGCAGCAATACGTAACTGAACGAAGTACAGGAAAAAAAA";
        let rc = Iupac::reverse_complement(fwd);

        let text = b"TGTTATATTTCCCTGTACTTCGTTCCAGTTATTTTTATGCAAAAAACCGGTGTTTAACCACCACTGCCATGTATCAAAGTACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCAACAGGAAAACTATTTTCTGCAGGCATTTTGCCCGGTTTAGTGCTGGGGCTGGCATTAATAATCTATGCAGTAATCGTCGCCCATAAAAGGCTACGGCGGCCTGCCTAAA".to_vec();

        println!("TEXT LEN: {}", text.len());
        println!("FWD LEN: {}", fwd.len());

        let mut searcher = Searcher::<Iupac>::new_rc_with_overhang(0.5);
        let fwd_matches = searcher.search(fwd, &text, 44);
        let rc_matches = searcher.search(&rc, &text, 44);

        // Print matches for debugging
        println!("Forward matches:");
        for m in fwd_matches.iter() {
            println!("  {:?}", m.without_cigar());
        }
        println!("\nReverse complement matches:");
        for m in rc_matches.iter() {
            println!("  {:?}", m.without_cigar());
        }

        assert_eq!(
            fwd_matches.len(),
            rc_matches.len(),
            "Forward and reverse complement searches should find the same number of matches"
        );

        // For each fwd, there should be matching rc (just strand difference)
        for fwd_match in fwd_matches.iter() {
            let matching_rc = rc_matches.iter().find(|rc_match| {
                rc_match.start == fwd_match.start
                    && rc_match.end == fwd_match.end
                    && rc_match.cost == fwd_match.cost
            });
            assert!(
                matching_rc.is_some(),
                "No matching RC match found for forward match: {:?}",
                fwd_match.without_cigar()
            );
        }
    }
}
