use rand::Rng;
use sassy::search::Searcher;
use serde::Deserialize;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(Deserialize)]
struct Config {
    alphas: Vec<f32>,
}

// Generate random dna sequence of length "length" with alphabet
fn generate_random_dna_sequence(length: usize) -> Vec<u8> {
    let mut rng = rand::rng();
    (0..length)
        .map(|_| b"ACGT"[rng.random_range(0..4)])
        .collect()
}

pub fn run(config: &str) {
    let toml_str = fs::read_to_string(config).unwrap();
    let config: Config = toml::from_str(&toml_str).unwrap();

    let query = generate_random_dna_sequence(20);
    let text = generate_random_dna_sequence(100);

    println!(
        "Query (len: {}): {:?}",
        query.len(),
        String::from_utf8(query.clone()).unwrap()
    );
    println!(
        "Text (len: {}): {:?}",
        text.len(),
        String::from_utf8(text.clone()).unwrap()
    );

    // Create output file
    let mut file = File::create("overhang_results.csv").unwrap();
    writeln!(file, "alpha,end_pos,cost").unwrap();

    for alpha in config.alphas {
        let mut searcher = Searcher::<sassy::profiles::Iupac>::new_fwd_with_overhang(alpha);
        let matches = searcher.search_all(&query, &text, query.len() * 2);
        for m in matches {
            // write alpha, end_pos, cost
            writeln!(file, "{},{},{}", alpha, m.end.1, m.cost).unwrap();
        }
    }
}
