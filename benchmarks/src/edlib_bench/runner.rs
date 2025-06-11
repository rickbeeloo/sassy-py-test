use crate::edlib_bench::edlib::*;
use crate::edlib_bench::grid::*;
use crate::edlib_bench::sim_data::*;
use sassy::search::{Match, Searcher};
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};

macro_rules! time_it {
    ($label:expr, $expr:expr, $iters:expr, $pairs:expr) => {{
        let label = $label;
        let mut times = Vec::with_capacity($iters);
        let mut final_result = None;

        for (q, t, _) in $pairs.iter() {
            let start = std::time::Instant::now();
            let r = std::hint::black_box($expr(q, t));
            let elapsed = start.elapsed();
            times.push(elapsed.as_micros() as f64);
            final_result = Some(r); // So just the last result, assuming on random data it's approximately the same across pairs
        }

        let mean = times.iter().sum::<f64>() / times.len() as f64;
        eprintln!("{label:>10} : {:.3}ms", mean / 1000.0);
        (final_result.unwrap(), mean)
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
        .open(grid.output_file())
        .expect(&format!("Unable to open {}", grid.output_file()));

    let mut writer = BufWriter::new(file);

    // Write header
    writeln!(
        writer,
        "query_length,text_length,k,edlib_matches,sassy_matches,max_edits,bench_iter,alphabet,profile,rc,edlib_ms,sassy_ms"
    ).unwrap();

    // Get combinations
    for param_set in grid.all_combinations() {
        println!("Param set: {:?}", param_set);

        let num_matches = param_set.matches;
        let bench_iter = param_set.bench_iter;

        // Generate all query-text pairs upfront
        let pairs: Vec<_> = (0..bench_iter)
            .map(|_| {
                generate_query_and_text_with_matches(
                    param_set.query_length,
                    param_set.text_length,
                    num_matches,
                    param_set.max_edits,
                    param_set.max_edits,
                    &param_set.alphabet,
                )
            })
            .collect();

        // Running Edlib
        let (edlib_matches, edlib_mean_ms) = if param_set.edlib {
            let edlib_config = get_edlib_config(param_set.k as i32, &param_set.alphabet);
            let (r, ms) = time_it!(
                "edlib",
                |q, t| run_edlib(q, t, &edlib_config),
                bench_iter,
                &pairs
            );
            let edlib_matches = r.startLocations.unwrap_or(vec![]);
            (edlib_matches, ms)
        } else {
            (vec![], 0.0)
        };

        // Get the correct search function (not timed)
        let mut search_fn = get_search_fn(&param_set);

        // Now time the search
        let (sassy_matches, sassy_mean_ms) = time_it!(
            "sassy",
            |q, t| search_fn(q, t, param_set.k),
            bench_iter,
            &pairs
        );

        if param_set.edlib {
            println!("Edlib matches: {:?}", edlib_matches.len());
        }
        println!("Sassy matches: {:?}", sassy_matches.len());

        if param_set.verbose {
            println!("Edlib matches");
            for loc in edlib_matches.iter() {
                println!("{}", loc);
            }
            println!("Sassy matches");
            for loc in sassy_matches.iter() {
                println!("{:?}", loc);
            }
        }

        // Write row to CSV
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{},{:.6},{:.6}",
            param_set.query_length,
            param_set.text_length,
            param_set.k,
            edlib_matches.len(),
            sassy_matches.len(),
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

type SearchFn = Box<dyn FnMut(&[u8], &[u8], usize) -> Vec<Match>>;

fn get_search_fn(param_set: &ParamSet) -> SearchFn {
    let rc = match param_set.rc {
        "withrc" => true,
        "withoutrc" => false,
        x => panic!("Unsupported rc config: {x}"),
    };
    match param_set.profile {
        // IUPAC profile
        "iupac" => {
            let mut searcher = Searcher::<sassy::profiles::Iupac>::new(rc, None);
            Box::new(move |q, t, k| searcher.search(&q, &t, k))
        }

        // DNA profile
        "dna" => {
            let mut searcher = Searcher::<sassy::profiles::Dna>::new(rc, None);
            Box::new(move |q, t, k| searcher.search(&q, &t, k))
        }

        // ASCII profile
        "ascii" => {
            let mut searcher = Searcher::<sassy::profiles::Ascii>::new(rc, None);
            Box::new(move |q, t, k| searcher.search(&q, &t, k))
        }

        _ => panic!(
            "Unsupported combination: {:?} {:?} {:?}",
            param_set.profile, param_set.rc, param_set.alphabet
        ),
    }
}
