use crate::edlib_bench::edlib::*;
use crate::edlib_bench::grid::*;
use crate::edlib_bench::sim_data::*;
use sassy::search::{Match, Searcher};
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};

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
pub fn run(grid_config: &str) {
    // Read grid config file for benching
    let grid = read_grid(grid_config).expect("Invalid grid config");

    // Open output file and write header
    let file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open("results.csv")
        .expect("Unable to open results.csv");
    let mut writer = BufWriter::new(file);

    // Write header
    writeln!(
        writer,
        "query_length,text_length,k,match_fraction,max_edits,bench_iter,alphabet,profile,rc,bound,edlib_ms,sassy_ms"
    ).unwrap();

    // Get combinations
    for param_set in grid.all_combinations() {
        println!("Param set: {:?}", param_set);

        // Aboslute number of matches
        let num_matches = (param_set.match_fraction * param_set.text_length as f64).ceil() as usize;

        // Generating random data
        let (q, t, locs) = generate_query_and_text_with_matches(
            param_set.query_length,
            param_set.text_length,
            1,
            param_set.max_edits,
            &param_set.alphabet,
        );

        // Number of  bench iterations;
        let bench_iter = 50;

        // Running Edlib
        let edlib_config = get_edlib_config(param_set.k as i32, &param_set.alphabet);
        let (edlib_result, edlib_mean_ms) =
            time_it!("edlib", run_edlib(&q, &t, &edlib_config), bench_iter);

        // Get the correct search function (not timed)
        let search_fn = get_search_fn(&param_set);

        // Now time the search
        let (sassy_result, sassy_mean_ms) =
            time_it!("sassy", search_fn(&q, &t, param_set.k), bench_iter);

        // Print number of matches to validate
        let edlib_matches = edlib_result.startLocations.unwrap_or(vec![]).len();
        let sassy_matches = sassy_result.len();
        println!("Edlib matches: {:?}", edlib_matches);
        println!("Sassy matches: {:?}", sassy_matches);

        // Write row to CSV
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{:.6},{:.6}",
            param_set.query_length,
            param_set.text_length,
            param_set.k,
            param_set.match_fraction,
            param_set.max_edits,
            param_set.bench_iter,
            format!("{:?}", param_set.alphabet),
            param_set.profile,
            param_set.rc,
            edlib_mean_ms,
            sassy_mean_ms
        )
        .unwrap();
    }
    // Ensure all data is written
    writer.flush().unwrap();
}

type SearchFn = fn(&[u8], &[u8], usize) -> Vec<Match>;

fn get_search_fn(param_set: &ParamSet) -> SearchFn {
    match (param_set.profile, param_set.rc, &param_set.alphabet) {
        // IUPAC profile
        ("iupac", "withrc", _) => {
            |q, t, k| Searcher::<sassy::profiles::Iupac, true, false>::new().search(&q, &t, k)
        }
        ("iupac", "withoutrc", _) => {
            |q, t, k| Searcher::<sassy::profiles::Iupac, false, false>::new().search(&q, &t, k)
        }

        // DNA profile
        ("dna", "withrc", _) => {
            |q, t, k| Searcher::<sassy::profiles::Dna, true, false>::new().search(&q, &t, k)
        }
        ("dna", "withoutrc", _) => {
            |q, t, k| Searcher::<sassy::profiles::Dna, false, false>::new().search(&q, &t, k)
        }

        // ASCII profile
        ("ascii", "withrc", _) => {
            |q, t, k| Searcher::<sassy::profiles::Ascii, true, false>::new().search(&q, &t, k)
        }
        ("ascii", "withoutrc", _) => {
            |q, t, k| Searcher::<sassy::profiles::Ascii, false, false>::new().search(&q, &t, k)
        }

        _ => panic!(
            "Unsupported combination: {:?} {:?} {:?}",
            param_set.profile, param_set.rc, param_set.alphabet
        ),
    }
}
