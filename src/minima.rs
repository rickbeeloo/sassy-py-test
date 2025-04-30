use crate::delta_encoding::{V, VEncoding};
use pa_types::Cost;

// TODO: when reaching end of sequence, we can have -1-1-1,0 without going up
pub fn find_local_minima(
    initial_cost: Cost,
    p1: V<u64>,
    p2: V<u64>,
    p3: V<u64>,
) -> Vec<(usize, Cost)> {
    let mut valleys = Vec::new();
    let mut is_decreasing = false;
    let mut prev_cost = initial_cost;
    let mut cur_cost = initial_cost;

    for (word_idx, v) in [p1, p2, p3].iter().enumerate() {
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
    valleys
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::profiles::*;
    use crate::search::*;

    #[test]
    fn test_trace_dna() {
        let query = b"ACTG";
        let mock_chunk = [b'G'; 32];

        let mut match_chunk = [b'G'; 32];
        match_chunk[0] = b'A';
        match_chunk[1] = b'C';
        match_chunk[2] = b'A';
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

        find_below_threshold(query, 1, &deltas, &mut positions, &mut costs);

        for (pos, cost) in positions.iter().zip(costs.iter()) {
            println!("\nPos: {:?}, Cost: {:?}", pos, cost);
            let chunk_pos = pos / 64;
            let minima = find_local_minima(
                *cost,
                deltas[chunk_pos - 1],
                deltas[chunk_pos],
                deltas[chunk_pos + 1],
            );
            println!("\t-Minima: {:?}", minima);
        }
    }
}
