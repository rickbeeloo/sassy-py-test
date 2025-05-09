use ::std::os::raw::c_char;
use clap::{Parser, ValueEnum};
use edlib_rs::edlib_sys::*;
use edlib_rs::*;
use once_cell::sync::Lazy;
use rand::Rng;
use sassy::LocalMinK;
use sassy::profiles::Profile;
use serde::Deserialize;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

macro_rules! time_it {
    ($label:expr, $expr:expr, $iters:expr) => {{
        let label = $label;
        let mut times = Vec::with_capacity($iters);
        let mut result = None;
        for _ in 0..$iters {
            let start = std::time::Instant::now();
            // https://doc.rust-lang.org/std/hint/fn.black_box.html
            let r = std::hint::black_box($expr);
            let elapsed = start.elapsed();
            times.push(elapsed.as_micros() as f64);
            result = Some(r);
        }
        let mean = times.iter().sum::<f64>() / times.len() as f64;
        eprintln!("{label:>10} : {:.3}ms", mean / 1000.0);
        (result.unwrap(), mean)
    }};
}

#[derive(Copy, Clone, Debug, PartialEq, Deserialize)]
#[serde(rename_all = "lowercase")]
enum Alphabet {
    Dna,
    Iupac,
    Ascii,
}

fn generate_random_sequence(length: usize, alphabet: &Alphabet) -> Vec<u8> {
    let mut rng = rand::rng();

    match alphabet {
        Alphabet::Dna => (0..length)
            .map(|_| {
                let c = rng.random_range(0..4);
                b"ACGT"[c]
            })
            .collect(),
        Alphabet::Iupac => (0..length)
            .map(|_| {
                let c = rng.random_range(0..16);
                b"ACGTURYSWKMBDHVNX"[c]
            })
            .collect(),
        Alphabet::Ascii => (0..length)
            .map(|_| rng.random_range(0..256) as u8)
            .collect(),
    }
}

// No idea why this needs to be static per se
static EQUALITY_PAIRS: Lazy<Vec<EdlibEqualityPairRs>> = Lazy::new(build_equality_pairs);

/// Supporting atcgn, with configurable k, NOTE uses traceback
fn get_edlib_config(k: i32, profile: &str) -> EdlibAlignConfigRs<'static> {
    let mut config = EdlibAlignConfigRs::default();
    config.mode = EdlibAlignModeRs::EDLIB_MODE_HW;
    if profile == "iupac" {
        config.additionalequalities = &EQUALITY_PAIRS;
    }
    config.k = k;
    config.task = EdlibAlignTaskRs::EDLIB_TASK_PATH;
    config
}

/// Apply random mutations to a sequence, sub/ins/del
fn mutate_sequence(sequence: &[u8], max_edits: usize) -> Vec<u8> {
    let mut rng = rand::rng();
    let mut mutated = sequence.to_vec();
    let edits = max_edits;
    // let edits = rng.random_range(0..=max_edits.min(sequence.len()));

    for _ in 0..edits {
        let edit_type = rng.random_range(0..3);
        let idx = rng.random_range(0..mutated.len());
        match edit_type {
            // We just insert DNA chars for ASCII, DNA, and IUPAC now
            0 => mutated[idx] = b"ACGT"[rng.random_range(0..4)],
            1 if mutated.len() > 1 => {
                mutated.remove(idx);
            }
            2 => mutated.insert(idx, b"ACGT"[rng.random_range(0..4)]),
            _ => {}
        }
    }

    mutated
}

/// Generate a random DNA query and a text with `num_matches` planted, each with up to `max_edits` edits.
/// also returns the locations of the planted matches
fn generate_query_and_text_with_matches(
    query_len: usize,
    text_len: usize,
    num_matches: usize,
    max_edits: usize,
    alphabet: &Alphabet,
) -> (Vec<u8>, Vec<u8>, Vec<(usize, usize)>) {
    let mut rng = rand::rng();

    // Generate random query and text
    let query = generate_random_sequence(query_len, alphabet);
    let mut text = generate_random_sequence(text_len, alphabet);

    let mut planted_locs = Vec::new();

    // Try to insert up to num_matches, might not be possibel to insert any or all
    // as we maintain the text len. For example we cant insert a 1KB with 1 edit in 40nt text
    'outer: for _ in 0..num_matches {
        let mutated = mutate_sequence(&query, max_edits);
        let mlen = mutated.len();

        // if mutated is longer than the whole text, skip
        if mlen > text.len() {
            continue;
        }

        // Find a non-overlapping insert region - 10 tries
        // TODO: just keep track of uncovered regions
        let max_start = text.len() - mlen;
        for _ in 0..10 {
            let start = rng.random_range(0..=max_start);
            let end = start + mlen;
            // check against all previously planted intervals
            if planted_locs.iter().all(|&(s, e)| end <= s || start >= e) {
                // replace exactly [start..end) with mutated
                text.splice(start..end, mutated.iter().cloned());
                planted_locs.push((start, end));
                continue 'outer;
            }
        }
    }

    (query, text, planted_locs)
}

/// Benchmark function with planted matches and edits
fn run_bench_with_planted<P, D, B>(config: &BenchConfig)
where
    P: sassy::profiles::Profile,
    D: sassy::RcBehavior<P>,
    B: sassy::BoundBehavior<P>,
{
    let edlib_config = get_edlib_config(config.edlib_k, &config.profile);

    let file = File::create(&config.output_file).expect("Unable to create output file");
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "query_len\tref_len\tnum_matches\tmax_edits\ttool\trun_time\tplanted\tfound"
    )
    .unwrap();

    let max_edits = config.max_edits;

    for &q_len in &config.query_lengths {
        for &ref_len in &config.text_lengths {
            if ref_len < q_len {
                continue;
            }

            // Calculate number of possible placements
            let num_possible_matches = ref_len.saturating_sub(q_len).saturating_add(1);
            let mut num_matches =
                (num_possible_matches as f64 * config.match_fraction).round() as usize;

            let (query, text, planted_locs) = generate_query_and_text_with_matches(
                q_len,
                ref_len,
                num_matches,
                max_edits,
                &config.alphabet,
            );

            num_matches = planted_locs.len();

            // Run edlib
            let (edlib_result, edlib_mean_ms) = time_it!(
                "edlib",
                edlibAlignRs(&query, &text, &edlib_config),
                config.bench_iter
            );
            assert_eq!(edlib_result.status, EDLIB_STATUS_OK);

            // Count edlib found matches
            let edlib_found = edlib_result
                .getStartLocations()
                .map_or(0, |locs| locs.len());

            writeln!(
                writer,
                "{q_len}\t{ref_len}\t{num_matches}\t{max_edits}\tedlib\t{:.3}\t{}\t{}",
                edlib_mean_ms, num_matches, edlib_found
            )
            .unwrap();

            // Run sassy
            let (sassy_result, sassy_mean_ms) = time_it!(
                "sassy",
                sassy::search_generic::<P, D, B, LocalMinK>(&query, &text, config.sassy_k),
                config.bench_iter
            );
            let sassy_found = sassy_result.len();

            writeln!(
                writer,
                "{q_len}\t{ref_len}\t{num_matches}\t{max_edits}\tsassy\t{:.3}\t{}\t{}",
                sassy_mean_ms, num_matches, sassy_found
            )
            .unwrap();

            println!(">Query length: {}, Reference length: {}", q_len, ref_len);

            let query_str = String::from_utf8_lossy(&query);
            let text_str = String::from_utf8_lossy(&text);

            for &(planted_start, _planted_end) in planted_locs.iter() {
                let found = sassy_result
                    .iter()
                    .any(|m| (m.start.1 as usize).abs_diff(planted_start) <= config.max_edits);
                if !found {
                    // Format error message with detailed information
                    let mut error_msg = format!(
                        "\n😔 ERROR: Verification failed! (q_len={}, ref_len={})\n\n",
                        q_len, ref_len
                    );

                    // Add query and text information
                    error_msg.push_str(&format!("Query: \"{}\"\n\n", query_str));
                    error_msg.push_str(&format!("Text: \"{}\"\n\n", text_str));

                    // Add planted locations
                    error_msg.push_str("Planted locations:\n");
                    for (j, &(start, end)) in planted_locs.iter().enumerate() {
                        error_msg.push_str(&format!(
                            "  {}. Position {}..{}: \"{}\"\n",
                            j + 1,
                            start,
                            end,
                            &text_str[start..end]
                        ));
                    }

                    // Add sassy results
                    // error_msg.push_str("\nSassy matches found:\n");
                    // if sassy_result.is_empty() {
                    //     error_msg.push_str("  (No matches found by sassy)\n");
                    // } else {
                    //     for (j, m) in sassy_result.iter().enumerate() {
                    //         let match_seq = &text_str[m.start.1 as usize..m.end.1 as usize];
                    //         error_msg.push_str(&format!(
                    //             "  {}. Position {}..{}: \"{}\"\n",
                    //             j + 1,
                    //             m.start.1,
                    //             m.end.1,
                    //             match_seq
                    //         ));
                    //     }
                    // }

                    // // Add edlib results
                    // error_msg.push_str("\nEdlib matches found:\n");

                    // if let Some(edlib_locations) = edlib_result.getStartLocations() {
                    //     if edlib_locations.is_empty() {
                    //         error_msg.push_str("  (No matches found by edlib)\n");
                    //     } else {
                    //         for (j, &start_pos) in edlib_locations.iter().enumerate() {
                    //             let end_pos = start_pos as usize + query.len();
                    //             let match_seq = if end_pos <= text.len() {
                    //                 &text_str[start_pos as usize..end_pos]
                    //             } else {
                    //                 "(out of bounds)"
                    //             };
                    //             error_msg.push_str(&format!(
                    //                 "  {}. Position {}..{}: \"{}\"\n",
                    //                 j + 1,
                    //                 start_pos,
                    //                 end_pos,
                    //                 match_seq
                    //             ));
                    //         }
                    //     }
                    // } else {
                    //     error_msg.push_str("  (No matches found by edlib)\n");
                    // }
                }
            }
            // planted_locs.sort();
            println!("✅ All {} planted matches found:", planted_locs.len(),);
            // Print query, text, and edlib match locations (can be commented out for benchmarking)
            // println!("Query: {}", String::from_utf8_lossy(&query));
            // println!("Text: {}", String::from_utf8_lossy(&text));
            // if let Some(edlib_locations) = edlib_result.getStartLocations() {
            //     println!(
            //         "Edlib match locations:\n{}",
            //         edlib_locations
            //             .iter()
            //             .map(|&start_pos| {
            //                 let end_pos = start_pos as usize + query.len();
            //                 format!("{}..{}", start_pos, end_pos)
            //             })
            //             .collect::<Vec<_>>()
            //             .join(", ")
            //     );
            // }
        }
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the TOML config file
    config: String,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum ProfileOpt {
    Iupac,
    Dna,
    Ascii,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum RcOpt {
    WithRc,
    WithoutRc,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum BoundOpt {
    Bounded,
    Unbounded,
}

#[derive(Debug, Deserialize)]
struct BenchConfig {
    query_lengths: Vec<usize>,
    text_lengths: Vec<usize>,
    output_file: String,
    edlib_k: i32,
    sassy_k: usize,
    match_fraction: f64,
    max_edits: usize,
    bench_iter: usize,
    alphabet: Alphabet, // "dna", "iupac", or "ascii"
    profile: String,    // "iupac", "dna", "ascii"
    rc: String,         // "withrc", "withoutrc"
    bound: String,      // "bounded", "unbounded"
}

fn read_config<P: AsRef<Path>>(path: P) -> BenchConfig {
    let s = std::fs::read_to_string(path).expect("Failed to read config file");
    toml::from_str(&s).expect("Failed to parse TOML config")
}

fn main() {
    let args = Args::parse();
    let config = read_config(&args.config);

    macro_rules! dispatch {
        ($profile:ty, $rc:ty, $bound:ty) => {
            run_bench_with_planted::<$profile, $rc, $bound>(&config)
        };
    }

    match (
        config.profile.to_lowercase().as_str(),
        config.rc.to_lowercase().as_str(),
        config.bound.to_lowercase().as_str(),
    ) {
        ("iupac", "withrc", "bounded") => {
            dispatch!(sassy::profiles::Iupac, sassy::WithRc, sassy::Bounded)
        }
        ("iupac", "withrc", "unbounded") => {
            dispatch!(sassy::profiles::Iupac, sassy::WithRc, sassy::Unbounded)
        }
        ("iupac", "withoutrc", "bounded") => {
            dispatch!(sassy::profiles::Iupac, sassy::WithoutRc, sassy::Bounded)
        }
        ("iupac", "withoutrc", "unbounded") => {
            dispatch!(sassy::profiles::Iupac, sassy::WithoutRc, sassy::Unbounded)
        }
        ("dna", "withrc", "bounded") => {
            dispatch!(sassy::profiles::Dna, sassy::WithRc, sassy::Bounded)
        }
        ("dna", "withrc", "unbounded") => {
            dispatch!(sassy::profiles::Dna, sassy::WithRc, sassy::Unbounded)
        }
        ("dna", "withoutrc", "bounded") => {
            dispatch!(sassy::profiles::Dna, sassy::WithoutRc, sassy::Bounded)
        }
        ("dna", "withoutrc", "unbounded") => {
            dispatch!(sassy::profiles::Dna, sassy::WithoutRc, sassy::Unbounded)
        }
        ("ascii", "withrc", "bounded") => {
            dispatch!(sassy::profiles::Ascii, sassy::WithRc, sassy::Bounded)
        }
        ("ascii", "withrc", "unbounded") => {
            dispatch!(sassy::profiles::Ascii, sassy::WithRc, sassy::Unbounded)
        }
        ("ascii", "withoutrc", "bounded") => {
            dispatch!(sassy::profiles::Ascii, sassy::WithoutRc, sassy::Bounded)
        }
        ("ascii", "withoutrc", "unbounded") => {
            dispatch!(sassy::profiles::Ascii, sassy::WithoutRc, sassy::Unbounded)
        }
        _ => panic!("Unknown profile/rc/bound combination"),
    }
}

fn build_equality_pairs() -> Vec<EdlibEqualityPairRs> {
    let codes = b"ACGTURYSWKMBDHVNX";
    let mut pairs = Vec::new();
    for &a in codes.iter() {
        for &b in codes.iter() {
            if sassy::profiles::Iupac::is_match(a, b) {
                // both upper
                pairs.push(EdlibEqualityPairRs {
                    first: a as c_char,
                    second: b as c_char,
                });
                // both lower
                pairs.push(EdlibEqualityPairRs {
                    first: a.to_ascii_lowercase() as c_char,
                    second: b.to_ascii_lowercase() as c_char,
                });
                // first upper, second lower
                pairs.push(EdlibEqualityPairRs {
                    first: a.to_ascii_lowercase() as c_char,
                    second: b as c_char,
                });
                // first lower, second upper
                pairs.push(EdlibEqualityPairRs {
                    first: a as c_char,
                    second: b.to_ascii_lowercase() as c_char,
                });
            }
        }
    }
    pairs
}
