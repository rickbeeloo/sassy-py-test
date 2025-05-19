use crate::edlib_bench::sim_data::Alphabet;
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
pub struct GridConfig {
    pub query_lengths: Vec<usize>,
    pub text_lengths: Vec<usize>,
    pub k: Vec<f32>,
    pub match_fraction: Vec<f64>,
    pub bench_iter: Vec<usize>,
    pub alphabet: Vec<Alphabet>,
    pub profile: Vec<String>,
    pub rc: Vec<String>,
    pub edlib: bool,
}

#[derive(Clone, Debug)]
pub struct ParamSet<'a> {
    pub query_length: usize,
    pub text_length: usize,
    pub k: usize,
    pub match_fraction: f64,
    pub max_edits: usize,
    pub bench_iter: usize,
    pub alphabet: Alphabet,
    pub profile: &'a str,
    pub rc: &'a str,
    pub edlib: bool,
}

impl GridConfig {
    /// Returns an iterator over all parameter combinations as tuples.
    pub fn all_combinations<'a>(&'a self) -> impl Iterator<Item = ParamSet> + 'a {
        self.query_lengths.iter().flat_map(move |&ql| {
            self.text_lengths.iter().flat_map(move |&tl| {
                self.k.iter().flat_map(move |&k| {
                    self.match_fraction.iter().flat_map(move |&mf| {
                        self.bench_iter.iter().flat_map(move |&bi| {
                            self.alphabet.iter().flat_map(move |&a| {
                                self.profile
                                    .iter()
                                    .flat_map(move |p| {
                                        let k = if k < 1.0 {
                                            (k * ql as f32).round() as usize
                                        } else {
                                            k as usize
                                        };
                                        // Only allow matching profile/alphabet pairs
                                        self.rc.iter().map(move |rc| ParamSet {
                                            query_length: ql,
                                            text_length: tl,
                                            k,
                                            match_fraction: mf,
                                            max_edits: k,
                                            bench_iter: bi,
                                            alphabet: a,
                                            profile: p.as_str(),
                                            rc: rc.as_str(),
                                            edlib: self.edlib,
                                        })
                                    })
                                    .filter(|param| {
                                        (param.alphabet == Alphabet::Dna && param.profile == "dna")
                                            || (param.alphabet == Alphabet::Dna
                                                && param.profile == "iupac")
                                            || (param.alphabet == Alphabet::Iupac
                                                && param.profile == "iupac")
                                            || (param.alphabet == Alphabet::Ascii
                                                && param.profile == "ascii")
                                    })
                            })
                        })
                    })
                })
            })
        })
    }
}

pub fn read_grid(path: &str) -> Result<GridConfig, Box<dyn std::error::Error>> {
    let toml_str = fs::read_to_string(path)?;
    let config: GridConfig = toml::from_str(&toml_str)?;
    Ok(config)
}
