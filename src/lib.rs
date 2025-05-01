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

use log::debug;
pub use minima::find_below_threshold;
use minima::find_local_minima;
use pa_types::{Cigar, Cost, Pos};
pub use search::{search_positions, search_positions_bounded};

use profiles::Profile;
use trace::{get_trace, simd_fill};

#[derive(Debug, Clone)]
pub struct Match {
    pub start: Pos,
    pub end: Pos,
    pub cost: Cost,
    pub cigar: Cigar,
}

pub fn search<P: Profile>(query: &[u8], text: &[u8], k: usize) -> Vec<Match> {
    let mut deltas = vec![];
    search_positions::<P>(query, text, &mut deltas);
    debug!("TEXT Len {}", text.len());
    let matches = find_local_minima(query, &deltas, k as Cost);
    debug!("matches {matches:?}");

    if !matches.is_empty() {
        debug!("Num matches: {}", matches.len());
    }

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
            // FIXME: Adjust returned positions for start-index offset.
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
