use crate::{
    delta_encoding::{V, VEncoding},
    search::Deltas,
};
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

pub fn find_local_minima(
    query: &[u8],
    deltas: &mut Deltas,
    k: Cost,
    text_len: usize,
    prefix_min_bytes: &mut [u8; 8],
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
        let (min, delta) = prefix_min(v.1, prefix_min_bytes);
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
    prefix_min_bytes: &mut [u8; 8],
) {
    let mut cur_cost = query.len() as Cost;
    for (i, v) in deltas.iter().enumerate() {
        let (min, delta) = prefix_min(v.1, prefix_min_bytes);
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

static TABLE_MIN: [i8; 256] = {
    let mut a = [0i8; 256];
    let mut i = 0;
    while i < 256 {
        a[i] = TABLE[i].0;
        i += 1;
    }
    a
};
static TABLE_END: [i8; 256] = {
    let mut a = [0i8; 256];
    let mut i = 0;
    while i < 256 {
        a[i] = TABLE[i].1;
        i += 1;
    }
    a
};

#[inline(always)]
pub fn prefix_min(v: V<u64>, bytes: &mut [u8; 8]) -> (i8, i8) {
    let (p, m) = v.pm();
    let delta = p | m;
    let num_p = p.count_ones() as i8;
    let num_m = m.count_ones() as i8;
    let deltas = unsafe { _pext_u64(m, delta) };

    // Split 64 bits in 8 bytes
    //let mut bytes = [0u8; 8];
    bytes.iter_mut().enumerate().for_each(|(i, byte)| {
        *byte = ((deltas >> (i * 8)) & 0xFF) as u8;
    });

    let mins = [
        TABLE_MIN[bytes[0] as usize],
        TABLE_MIN[bytes[1] as usize],
        TABLE_MIN[bytes[2] as usize],
        TABLE_MIN[bytes[3] as usize],
        TABLE_MIN[bytes[4] as usize],
        TABLE_MIN[bytes[5] as usize],
        TABLE_MIN[bytes[6] as usize],
        TABLE_MIN[bytes[7] as usize],
    ];

    let sums = [
        TABLE_END[bytes[0] as usize],
        TABLE_END[bytes[1] as usize],
        TABLE_END[bytes[2] as usize],
        TABLE_END[bytes[3] as usize],
        TABLE_END[bytes[4] as usize],
        TABLE_END[bytes[5] as usize],
        TABLE_END[bytes[6] as usize],
        TABLE_END[bytes[7] as usize],
    ];

    let mut cur = 0;
    let mut min = 0;

    for i in 0..8 {
        min = min.min(cur + mins[i]);
        cur += sums[i];
    }
    (min, num_p - num_m)
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::Rng;
    use std::time::{Duration, Instant};

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
//     fn test_simple_valley() {
//         let v1 = make_pattern(&[(0, -1), (1, -1), (3, 1), (4, 1)]);
//         let mut deltas = vec![v1];
//         let minima = find_local_minima(b"ATG", &mut deltas, 100, 64);
//         assert_eq!(minima, vec![(3, 1)]); // valley till position 2, cost 1
//     }

//     #[test]
//     fn test_cross_boundary_valley() {
//         let v1 = make_pattern(&[(62, -1), (63, -1)]);
//         let v2 = make_pattern(&[(1, 1), (2, 1)]);
//         let v3 = V(0, 0);
//         let mut deltas = vec![v1, v2, v3];
//         let minima = find_local_minima(b"ATG", &mut deltas, 100, 64 * 3);
//         assert_eq!(minima, vec![(65, 1)]); // at word boundary
//     }

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
