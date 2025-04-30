use crate::delta_encoding::{V, VEncoding};
use pa_types::Cost;

pub fn find_local_minima(
    initial_cost: Cost,
    minima_bits: &[V<u64>],
    at_right_end: bool, // If we want to report valleys at the right edge, set to true (i.e. end of seq)
) -> Vec<(usize, Cost)> {
    let mut valleys = Vec::new();
    let mut is_decreasing = false;
    let mut prev_cost = initial_cost;
    let mut cur_cost = initial_cost;

    for (word_idx, v) in minima_bits.iter().enumerate() {
        let (p, m) = v.pm();

        for bit in 0..64 {
            // Calculate cost
            let p_bit = ((p >> bit) & 1) == 1;
            let m_bit = ((m >> bit) & 1) == 1;
            cur_cost += (p_bit as Cost) - (m_bit as Cost);
            if cur_cost > prev_cost && is_decreasing {
                // Going up, but we were going down
                valleys.push((word_idx * 64 + bit - 1, prev_cost)); // relative prev pos
                is_decreasing = false;
            } else if cur_cost < prev_cost {
                is_decreasing = true;
            }
            prev_cost = cur_cost;
        }
    }

    // If we ended while decreasing and at_right_end is true, add the final position
    if at_right_end && is_decreasing {
        valleys.push((minima_bits.len() * 64 - 1, cur_cost));
    }
    valleys
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
        search::<Dna>(query, &seq, &mut deltas);
        println!("Deltas: {:?}", deltas);

        find_below_threshold(query, 0, &deltas, &mut positions, &mut costs);

        for (pos, cost) in positions.iter().zip(costs.iter()) {
            println!("\nPos: {:?}, Cost: {:?}", pos, cost);
            println!("Position: {}", pos);
            let chunk_pos = pos / 64 + 1;
            println!("Chunk pos: {}", chunk_pos);
            let minima = find_local_minima(*cost, &deltas[chunk_pos - 1..chunk_pos + 2], false);
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
        let v2 = V(0, 0);
        let v3 = V(0, 0);

        let minima = find_local_minima(10, &[v1, v2, v3], false);
        assert_eq!(minima, vec![(2, 8)]); // valley at position 2, cost 8
    }

    #[test]
    fn test_cross_boundary_valley() {
        // First word ends with down-down
        let v1 = make_pattern(&[(62, -1), (63, -1)]);
        // Second word starts with up-up
        let v2 = make_pattern(&[(1, 1), (2, 1)]);
        let v3 = V(0, 0);

        let minima = find_local_minima(10, &[v1, v2, v3], false);
        assert_eq!(minima, vec![(64, 8)]); // at word boundary 
    }

    #[test]
    fn test_multiple_valleys() {
        // Two valleys in same word
        let v1 = make_pattern(&[
            (0, -1),
            (1, -1),
            (2, 1),
            (3, 1), // first valley
            (4, -1),
            (5, -1),
            (6, 1),
            (7, 1), // second valley
        ]);
        let v2 = V(0, 0);
        let v3 = V(0, 0);

        let minima = find_local_minima(10, &[v1, v2, v3], false);
        assert_eq!(minima, vec![(1, 8), (5, 8)]);
    }

    #[test]
    fn test_valley_with_plateau() {
        // Valley with flat region in middle
        let v1 = make_pattern(&[
            (10, -1),
            (11, -1), // down
            // positions 12-14 0 (by default see make_pattern)
            (15, 1), // up
            (16, 1),
        ]);
        let v2 = V(0, 0);
        let v3 = V(0, 0);

        let minima = find_local_minima(10, &[v1, v2, v3], false);
        assert_eq!(minima, vec![(14, 8)]); // valley at end of plateau
    }

    #[test]
    fn test_long_cross_word_valley() {
        // Down at end of first word
        let v1 = make_pattern(&[(62, -1), (63, -1)]);
        // Zeros across entire second word
        let v2 = V(0, 0);
        // Up at start of third word
        let v3 = make_pattern(&[(0, 1), (1, 1)]);

        let minima = find_local_minima(10, &[v1, v2, v3], false);
        assert_eq!(minima, vec![(63 * 2 + 1, 8)]); // valley at end of second word
    }

    #[test]
    fn test_cost_calculation_simple() {
        // Starting at cost 10
        let v1 = make_pattern(&[
            (0, -1), // cost becomes 9
            (1, -1), // cost becomes 8 (valley)
            (2, 1),  // cost becomes 9
        ]);
        let v2 = V(0, 0);
        let v3 = V(0, 0);

        let minima = find_local_minima(10, &[v1, v2, v3], false);
        assert_eq!(minima, vec![(1, 8)]); // valley at position 1 with cost 8
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

        let minima = find_local_minima(20, &[v1, v2, v3], false);
        assert_eq!(minima, vec![(127, 18)]); // valley at end of second word with cost 18
    }

    #[test]
    fn test_at_right_end() {
        let v1 = make_pattern(&[(0, -1), (1, -1)]);
        let minima = find_local_minima(10, &[v1], true);
        assert_eq!(minima, vec![(63, 8)]); // We end with  a valley, right end true, so still report
        let minima = find_local_minima(10, &[v1], false);
        assert_eq!(minima, vec![]); // We end with a valley, right end false, so no report
    }
}
