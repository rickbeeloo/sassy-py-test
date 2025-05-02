use crate::delta_encoding::{V, VEncoding};
use pa_types::Cost;
use std::arch::x86_64::_pext_u64;

// Note: also reports minima at the end of the range.
#[allow(unused)] // only for testing
pub fn find_local_minima_slow(query: &[u8], deltas: &[V<u64>], k: Cost) -> Vec<(usize, Cost)> {
    let mut valleys = Vec::new();
    let mut is_decreasing = false;
    let mut prev_cost = query.len() as Cost;
    let mut cur_cost = query.len() as Cost;

    for (word_idx, v) in deltas.iter().enumerate() {
        let (p, m) = v.pm();

        for bit in 0..64 {
            // Calculate cost
            let p_bit = (p >> bit) & 1;
            let m_bit = (m >> bit) & 1;

            cur_cost += (p_bit as Cost) - (m_bit as Cost);
            if cur_cost > prev_cost && is_decreasing {
                if prev_cost <= k {
                    // Going up, but we were going down
                    let value = (word_idx * 64 + bit - 1, prev_cost);
                    // debug!("Push {value:?}");
                    valleys.push(value); // relative prev pos
                }
                is_decreasing = false;
            } else if cur_cost < prev_cost {
                is_decreasing = true;
            }
            prev_cost = cur_cost;
        }
    }

    // If we ended while decreasing, add the final position
    if cur_cost <= k && is_decreasing {
        valleys.push((deltas.len() * 64 - 1, cur_cost));
    }
    valleys
}

pub fn find_local_minima(query: &[u8], deltas: &[V<u64>], k: Cost) -> Vec<(usize, Cost)> {
    let mut prev_cost = query.len() as Cost;
    let mut cur_cost = query.len() as Cost;
    let mut all_valleys = Vec::new();
    let mut is_decreasing = false;

    for (word_idx, v) in deltas.iter().enumerate() {
        let (min, delta) = prefix_min(*v);
        if cur_cost + (min as Cost) <= k {
            let (p, m) = v.pm();
            // Get positions where cost changes occur
            let mut changes = p | m;
            while changes != 0 {
                let pos = changes.trailing_zeros() as usize;
                let p_bit = (p >> pos) & 1;
                let m_bit = (m >> pos) & 1;
                cur_cost += (p_bit as Cost) - (m_bit as Cost);
                if cur_cost > prev_cost && is_decreasing {
                    if prev_cost <= k {
                        all_valleys.push((word_idx * 64 + pos, prev_cost));
                    }
                    is_decreasing = false;
                } else if cur_cost < prev_cost {
                    is_decreasing = true;
                }
                prev_cost = cur_cost;
                // Clear the processed bit
                changes &= changes - 1;
            }
        } else {
            cur_cost += delta as Cost;
            prev_cost = cur_cost;
        }
    }

    // Add valley at right end if still decreasing
    if cur_cost <= k && is_decreasing {
        all_valleys.push((deltas.len() * 64, cur_cost));
    }

    all_valleys
}

pub fn find_below_threshold(
    query: &[u8],
    threshold: Cost,
    deltas: &[V<u64>],
    positions: &mut Vec<usize>,
    costs: &mut Vec<Cost>,
) {
    let mut cur_cost = query.len() as Cost;
    for (i, v) in deltas.iter().enumerate() {
        let (min, delta) = prefix_min(*v);
        if cur_cost + (min as Cost) <= threshold {
            positions.push(i * 64);
            // Cost at start of block
            costs.push(cur_cost as Cost);
        }
        cur_cost += delta as Cost;
    }
}

/// For each byte: (min_cost, end_cost)
/// Each 1 in a byte indicates -1.
/// Each 0 in a byte indicates +1.
const TABLE: [(i8, i8); 256] = {
    let mut table = [(0, 0); 256];

    let mut i = 0;
    while i < 256 {
        let mut min = 0;
        let mut cur = 0;
        let mut j = 0;
        while j < 8 {
            let bit = (i >> j) & 1;
            let delta = if bit == 1 { -1 } else { 1 };
            cur += delta;
            if cur < min {
                min = cur;
            }
            j += 1;
        }
        table[i] = (min, cur);
        i += 1;
    }

    table
};

/// Return the min_cost in the word, and total cost across the word.
#[inline(always)]
pub fn prefix_min(v: V<u64>) -> (i8, i8) {
    // extract only the relevant chars
    let (p, m) = v.pm();
    let delta = p | m;
    let num_p = p.count_ones();
    let num_m = m.count_ones();
    let deltas = unsafe { _pext_u64(m, delta) };
    let mut min = 0;
    let mut cur = 0;
    for i in 0..8 {
        let byte = (deltas >> (i * 8)) as u8 as usize;
        let (min_cost, end_cost) = TABLE[byte];
        min = min.min(cur + min_cost);
        cur += end_cost;
    }

    (min, num_p as i8 - num_m as i8)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::profiles::*;
    use crate::search::*;

    /// Create Vencoding from (position, delta) vec
    fn make_pattern(changes: &[(usize, i8)]) -> V<u64> {
        let mut p = 0u64;
        let mut m = 0u64;
        for &(pos, delta) in changes {
            assert!(pos < 64, "Position must be < 64");
            if delta > 0 {
                p |= 1u64 << pos;
            } else if delta < 0 {
                m |= 1u64 << pos;
            }
        }
        V(p, m)
    }

    #[test]
    fn test_trace_dna() {
        let query = b"ACTG";
        let mock_chunk = [b'G'; 64];

        let mut match_chunk = [b'G'; 64];
        match_chunk[0] = b'A';
        match_chunk[1] = b'C';
        match_chunk[2] = b'T';
        match_chunk[3] = b'G';

        let seq = [
            mock_chunk,
            mock_chunk,
            match_chunk,
            mock_chunk,
            mock_chunk,
            match_chunk,
            mock_chunk,
            mock_chunk,
        ]
        .concat();

        let mut deltas = vec![];
        let mut positions = vec![];
        let mut costs = vec![];
        search_positions::<Dna>(query, &seq, &mut deltas);
        println!("Deltas: {:?}", deltas);

        find_below_threshold(query, 0, &deltas, &mut positions, &mut costs);

        for (pos, cost) in positions.iter().zip(costs.iter()) {
            println!("\nPos: {:?}, Cost: {:?}", pos, cost);
            println!("Position: {}", pos);
            let chunk_pos = pos / 64 + 1;
            println!("Chunk pos: {}", chunk_pos);
            let minima = find_local_minima(query, &deltas[chunk_pos - 1..chunk_pos + 2], 100);
            for (rel_pos, cost) in minima {
                println!("\t-Minima: {:?} -cost: {}", pos + rel_pos, cost);
                println!(
                    "Minima seq: {:?}",
                    String::from_utf8_lossy(&seq[pos + rel_pos - 4..pos + rel_pos])
                );
            }
        }
    }

    #[test]
    fn test_simple_valley() {
        // Pattern: down-down-up-up at start
        let v1 = make_pattern(&[(0, -1), (1, -1), (3, 1), (4, 1)]);
        //Length is 3, -1, -1 = 1 edit at position 2
        let minima = find_local_minima(b"ATG", &[v1], 100);
        assert_eq!(minima, vec![(3, 1)]); // valley till position 2, cost 1
    }

    #[test]
    fn test_cross_boundary_valley() {
        // First word ends with down-down
        let v1 = make_pattern(&[(62, -1), (63, -1)]); // so 64 = 0
        // Second word starts with up-up
        let v2 = make_pattern(&[(1, 1), (2, 1)]); // meaning we go up at the start of 2nd word = 64  + 1 = 65
        let v3 = V(0, 0);
        // Again length is 3, -1, -1 = 1 edit at position 64 (after valley - crossing word)
        let minima = find_local_minima(b"ATG", &[v1, v2, v3], 100);
        assert_eq!(minima, vec![(65, 1)]); // at word boundary
    }

    #[test]
    fn test_multiple_valleys() {
        // Two valleys in same word
        let v1 = make_pattern(&[
            (0, -1), // 2
            (1, -1), // 1 <- first valley
            (2, 1),  // 2
            (3, 1),  // 3
            (4, -1), // 2
            (5, -1), // 1 <- second valley
            (6, 1),  // 2
            (7, 1),  // 3
        ]);
        let minima = find_local_minima(b"ATG", &[v1], 100);
        println!("Minima: {:?}", minima);
        assert_eq!(minima, vec![(2, 1), (6, 1)]);
    }

    #[test]
    fn test_valley_with_plateau() {
        // Valley with flat region in middle
        let v1 = make_pattern(&[
            (10, -1), // 2
            (11, -1), // 1 <- start valley
            // positions 12-14 0 (by default see make_pattern)
            (15, 1), // 1 <- end valley,  so first up position after valley
            (16, 1), // 2
        ]);
        let v2 = V(0, 0);
        let v3 = V(0, 0);

        let minima = find_local_minima(b"ATG", &[v1, v2, v3], 100);
        assert_eq!(minima, vec![(15, 1)]); // valley at end of plateau
    }

    #[test]
    fn test_long_cross_word_valley() {
        // Down at end of first word
        let v1 = make_pattern(&[(62, -1), (63, -1)]);
        // Zeros across entire second word
        let v2 = V(0, 0);
        // Up at start of third word
        let v3 = make_pattern(&[(0, 1), (1, 1)]);
        let minima = find_local_minima(b"ATG", &[v1, v2, v3], 100);
        assert_eq!(minima, vec![(64 * 2, 1)]); // valley at end of second word
    }

    #[test]
    fn test_cost_calculation_simple() {
        let v1 = make_pattern(&[
            (0, -1), // cost becomes 3-1 = 2
            (1, -1), // cost becomes 2-1 = 1
            (2, 1),  // cost becomes 1+1 = 2 (and here after valley)
        ]);
        let v2 = V(0, 0);
        let v3 = V(0, 0);

        let minima = find_local_minima(b"ATG", &[v1, v2, v3], 100);
        assert_eq!(minima, vec![(2, 1)]); // valley at position 1 with cost 8
    }

    #[test]
    fn test_cost_calculation_complex() {
        // Starting at cost 20
        // First word: ends with down(-1), down(-1)  // 20 -> 19 -> 18
        let v1 = make_pattern(&[(62, -1), (63, -1)]);

        // Second word: all zeros (plateau)          // stays at 18
        let v2 = V(0, 0);

        // Third word: starts with up(+1), up(+1)    // 18 -> 19 -> 20
        let v3 = make_pattern(&[(0, 1), (1, 1)]);

        let minima = find_local_minima(&[b'A'; 20], &[v1, v2, v3], 100);
        assert_eq!(minima, vec![(64 * 2, 18)]); // valley at end of second word with cost 18
    }

    #[test]
    fn test_at_right_end() {
        let v1 = make_pattern(&[(0, -1), (1, -1)]); // edits 1 3-2 = 1
        let minima = find_local_minima(b"ATG", &[v1], 100);
        assert_eq!(minima, vec![(64, 1)]); // We end with  a valley, right end true, so still report
    }
}
