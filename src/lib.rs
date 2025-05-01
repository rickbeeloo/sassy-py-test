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

pub use minima::find_below_threshold;
use minima::find_local_minima;
use pa_types::{Cost, Pos};
pub use search::{search_positions, search_positions_bounded};

use profiles::Profile;
use trace::{get_trace, simd_fill};

pub fn search<P: Profile>(query: &[u8], text: &[u8], k: usize) -> Vec<Vec<(usize, usize)>> {
    let mut deltas = vec![];
    search_positions::<P>(query, text, &mut deltas);
    let matches = find_local_minima(query.len() as Cost, &deltas);

    let mut traces = Vec::with_capacity(matches.len());

    let fill_len = query.len() + k;
    for matches in matches.chunks(4) {
        let mut text_slices = [[].as_slice(); 4];
        for i in 0..matches.len() {
            let end_post = matches[i].0;
            text_slices[i] = &text[end_post.saturating_sub(fill_len)..end_post];
        }
        // TODO: Reuse allocated costs.
        let costs = simd_fill::<P>(query, text_slices);

        for lane in 0..4 {
            // FIXME: Adjust returned positions for start-index offset.
            traces.push(get_trace::<P>(query, text_slices[lane], &costs[lane]));
        }
    }

    traces
}
