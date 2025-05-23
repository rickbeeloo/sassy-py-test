use crate::{
    delta_encoding::{V, VEncoding},
    search::Deltas,
};
use pa_types::Cost;
use std::{
    arch::x86_64::_pext_u64,
    simd::{
        Mask, Simd,
        cmp::SimdOrd,
        num::{SimdInt, SimdUint},
    },
};

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

pub fn find_local_minima(
    query: &[u8],
    deltas: &mut Deltas,
    k: Cost,
    text_len: usize,
) -> Vec<(usize, Cost)> {
    let mut prev_cost = query.len() as Cost;
    let mut cur_cost = query.len() as Cost;
    let mut all_valleys = Vec::new();
    let mut is_decreasing = false;

    // Handle overhang bits (bits beyond text_len)
    // We set all overhang bits to increasing edits to avoid false valleys
    let overhang = deltas.len() * 64 - text_len;
    if overhang > 0 {
        let mut remaining = overhang;
        for delta in deltas.iter_mut().rev() {
            if remaining >= 64 {
                delta.0 = Cost::MAX;
                delta.1 = V(u64::MAX, 0);
                remaining -= 64;
                if remaining == 0 {
                    break;
                }
            } else {
                // partial overhang
                let mask = (u64::MAX) << (64 - remaining);
                let (mut p, mut m) = delta.1.pm();
                p |= mask;
                m &= !mask;
                delta.1 = V(p, m);
                break;
            }
        }
    }

    for (word_idx, v) in deltas.iter().enumerate() {
        if v.0 == Cost::MAX {
            continue;
        }
        cur_cost = v.0;
        let (min, delta) = prefix_min(v.1.0, v.1.1);
        if cur_cost + (min as Cost) <= k {
            // FIXME?
            let (p, m) = v.1.pm();
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
    deltas: &Deltas,
    positions: &mut Vec<usize>,
    costs: &mut Vec<Cost>,
) {
    let mut cur_cost = query.len() as Cost;
    for (i, v) in deltas.iter().enumerate() {
        let (min, delta) = prefix_min(v.1.0, v.1.1);
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

pub fn find_all_minima(
    query: &[u8],
    deltas: &mut Deltas,
    k: Cost,
    text_len: usize,
) -> Vec<(usize, Cost)> {
    let mut cost = query.len() as Cost;
    let mut all_valleys = Vec::new();

    // Iterate through each block of 64 positions
    for (word_idx, v) in deltas.iter().enumerate() {
        if v.0 == Cost::MAX {
            continue;
        }

        // Reset cost at start of block
        cost = v.0;
        let (p, m) = v.1.pm();
        let base = word_idx * 64;

        // Step through every bit position
        for bit in 0..64 {
            let pos = base + bit;
            if pos > text_len {
                break;
            }

            // Check valley before applying change
            if cost <= k {
                all_valleys.push((pos, cost));
            }

            // Update cost if there's a change at this bit
            let p_bit = ((p >> bit) & 1) as Cost;
            let m_bit = ((m >> bit) & 1) as Cost;
            cost += p_bit;
            cost -= m_bit;
        }
    }

    // Check end-of-text position
    if cost <= k && text_len > 0 {
        all_valleys.push((text_len, cost));
    }

    all_valleys
}

/// Compute any prefix min <= k over 8 bytes via SIMD vectorized DP approach.
#[inline(always)]
pub fn prefix_min(p: u64, m: u64) -> (i8, i8) {
    // extract only the relevant chars
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

// split TABLE into two flat arrays for SIMD gather test
const TABLE_MIN: [Cost; 256] = {
    let mut a = [0 as Cost; 256];
    let mut i = 0;
    while i < 256 {
        a[i] = TABLE[i].0 as Cost;
        i += 1;
    }
    a
};

const TABLE_END: [Cost; 256] = {
    let mut a = [0 as Cost; 256];
    let mut i = 0;
    while i < 256 {
        a[i] = TABLE[i].1 as Cost;
        i += 1;
    }
    a
};

#[inline(always)]
pub fn prefix_min_k(start_cost: Cost, p: u64, m: u64, k: i32) -> (Cost, i8) {
    // We can use the best case scenario for the number of m-bits idea here.
    // after having a certain sum the "best" minima we can reach is always
    // current cost - remaining m-bits.
    // So we can abort early if this is already greater than k.
    let delta = p | m;
    let num_p = p.count_ones() as i8;
    let num_m = m.count_ones() as i8;

    // compress m-bits down into the delta‐positions
    let compressed_m = unsafe { _pext_u64(m, delta) };

    let lanes = Simd::from_array([
        (compressed_m) as u8,
        (compressed_m >> 8) as u8,
        (compressed_m >> 16) as u8,
        (compressed_m >> 24) as u8,
        (compressed_m >> 32) as u8,
        (compressed_m >> 40) as u8,
        (compressed_m >> 48) as u8,
        (compressed_m >> 56) as u8,
    ]);

    // count_ones on each lane in parallel
    let mut byte_pop: Simd<u8, 8> = lanes.count_ones();
    let counts = *byte_pop.as_array();
    byte_pop = Simd::from_array(counts);

    // Build a small suffix-sum table, pass as mut?
    let mut rem_pop = [0u8; 9];
    for i in (0..8).rev() {
        rem_pop[i] = rem_pop[i + 1] + byte_pop[i];
    }

    let mut min_cost = start_cost;
    let mut cur_cost = start_cost;

    for i in 0..8 {
        let byte = (compressed_m >> (i * 8)) as u8 as usize;
        let (tbl_min, tbl_end) = TABLE[byte];
        min_cost = min_cost.min(cur_cost + tbl_min as Cost);
        cur_cost += tbl_end as Cost;

        let remaining_m = rem_pop[i + 1] as i32;

        let min_i = min_cost;
        let best_possible = cur_cost - remaining_m;

        if min_i > k && best_possible > k {
            //  println!("Aborted at byte: {}, best_possible: {}", i, best_possible);
            break;
        }
    }

    (min_cost, num_p - num_m)
}

#[inline(always)]
pub fn prefix_min_k_simd(start_cost: Cost, p: u64, m: u64, k: i32) -> (Cost, i8) {
    let delta = p | m;
    let num_p = p.count_ones() as i8;
    let num_m = m.count_ones() as i8;
    let compressed_m = unsafe { _pext_u64(m, delta) };

    let byte_vec: Simd<u8, 8> = Simd::from_array([
        compressed_m as u8,
        (compressed_m >> 8) as u8,
        (compressed_m >> 16) as u8,
        (compressed_m >> 24) as u8,
        (compressed_m >> 32) as u8,
        (compressed_m >> 40) as u8,
        (compressed_m >> 48) as u8,
        (compressed_m >> 56) as u8,
    ]);
    let byte_pop: [u8; 8] = byte_vec.count_ones().to_array();

    // Fixme: we could pass as mut so we dont have to realloc
    let mut rem_pop = [0i32; 9];
    for i in (0..8).rev() {
        rem_pop[i] = rem_pop[i + 1] + (byte_pop[i] as i32);
    }

    let idxs: Simd<usize, 8> = byte_vec.cast(); // lanes = the raw bytes
    let mins_vec = Simd::gather_or_default(&TABLE_MIN, idxs);
    let ends_vec = Simd::gather_or_default(&TABLE_END, idxs);

    let mins_arr = mins_vec.to_array();
    let ends_arr = ends_vec.to_array();

    // We have to do this loop sequentially
    // (or simd scan prefix idea but that's also a lot of operations and no earlye exit)
    let mut min_cost = start_cost;
    let mut cur_cost = start_cost;

    for i in 0..8 {
        let tbl_min = mins_arr[i];
        let tbl_end = ends_arr[i];

        min_cost = min_cost.min(cur_cost + tbl_min);
        cur_cost += tbl_end;

        let best_possible = cur_cost - rem_pop[i + 1];
        if min_cost > k && best_possible > k {
            break;
        }
    }

    (min_cost, num_p - num_m)
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::Rng;

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

    fn random_pattern(len: usize) -> Vec<(usize, i8)> {
        let mut rng = rand::rng();
        let mut pattern = Vec::new();
        for _ in 0..len {
            let pos = rng.random_range(0..64);
            let delta = if rng.random_bool(0.5) { -1 } else { 1 };
            pattern.push((pos, delta));
        }
        pattern
    }

    fn random_query() -> Vec<u8> {
        let mut rng = rand::rng();
        let len = rng.random_range(1..20); // Random query length between 1 and 20
        let mut query = Vec::with_capacity(len);
        for _ in 0..len {
            query.push(rng.random_range(b'A'..b'Z')); // Random uppercase letters
        }
        query
    }
}

//     #[test]
//     fn test_multiple_valleys() {
//         let v1 = make_pattern(&[
//             (0, -1),
//             (1, -1),
//             (2, 1),
//             (3, 1),
//             (4, -1),
//             (5, -1),
//             (6, 1),
//             (7, 1),
//         ]);
//         let mut deltas = vec![v1];
//         let minima = find_local_minima(b"ATG", &mut deltas, 100, 64);
//         println!("Minima: {:?}", minima);
//         assert_eq!(minima, vec![(2, 1), (6, 1)]);
//     }

//     #[test]
//     fn test_valley_with_plateau() {
//         let v1 = make_pattern(&[(10, -1), (11, -1), (15, 1), (16, 1)]);
//         let v2 = V(0, 0);
//         let v3 = V(0, 0);
//         let mut deltas = vec![v1, v2, v3];
//         let minima = find_local_minima(b"ATG", &mut deltas, 100, 64 * 3);
//         assert_eq!(minima, vec![(15, 1)]); // valley at end of plateau
//     }

//     #[test]
//     fn test_long_cross_word_valley() {
//         let v1 = make_pattern(&[(62, -1), (63, -1)]);
//         let v2 = V(0, 0);
//         let v3 = make_pattern(&[(0, 1), (1, 1)]);
//         let mut deltas = vec![v1, v2, v3];
//         let minima = find_local_minima(b"ATG", &mut deltas, 100, 64 * 3);
//         assert_eq!(minima, vec![(64 * 2, 1)]); // valley at end of second word
//     }

//     #[test]
//     fn test_cost_calculation_simple() {
//         let v1 = make_pattern(&[(0, -1), (1, -1), (2, 1)]);
//         let v2 = V(0, 0);
//         let v3 = V(0, 0);
//         let mut deltas = vec![v1, v2, v3];
//         let minima = find_local_minima(b"ATG", &mut deltas, 100, 64 * 3);
//         assert_eq!(minima, vec![(2, 1)]); // valley at position 1 with cost 8
//     }

//     #[test]
//     fn test_cost_calculation_complex() {
//         let v1 = make_pattern(&[(62, -1), (63, -1)]);
//         let v2 = V(0, 0);
//         let v3 = make_pattern(&[(0, 1), (1, 1)]);
//         let mut deltas = vec![v1, v2, v3];
//         let minima = find_local_minima(&[b'A'; 20], &mut deltas, 100, 64 * 3);
//         assert_eq!(minima, vec![(64 * 2, 18)]); // valley at end of second word with cost 18
//     }

//     #[test]
//     fn test_at_right_end() {
//         let v1 = make_pattern(&[(0, -1), (1, -1)]);
//         let mut deltas = vec![v1];
//         let minima = find_local_minima(b"ATG", &mut deltas, 100, 64);
//         assert_eq!(minima, vec![(64, 1)]); // We end with  a valley, right end true, so still report
//     }
// }
