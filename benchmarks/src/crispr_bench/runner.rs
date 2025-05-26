use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::Duration;

use crate::crispr_bench::tools::*;
use serde::Deserialize;

#[derive(Deserialize)]
pub struct Config {
    pub sassy_path: String,
    pub swofinder_path: String,
    pub chopoff_path: String,
    pub target_file: String,
    pub out_dir: String,
    pub dists: Vec<usize>,
    pub threads: usize,
    pub chopoff_db_path: Option<String>,
    pub guides_file: String,
}

fn load_config(path: &str) -> Config {
    let s = std::fs::read_to_string(path).expect("Failed to read config TOML");
    toml::from_str(&s).expect("Failed to parse config TOML")
}

struct BenchmarkResult {
    tool_name: String,
    operation: String,
    duration: Duration,
    distance: usize,
    threads: usize,
    guide_seq: String,
}

// fn load_guide_sequences(path: &str) -> Result<Vec<String>, String> {
//     File::open(path)
//         .map_err(|e| format!("Unable to open {}: {}", path, e))
//         .and_then(|file| {
//             BufReader::new(file)
//                 .lines()
//                 .collect::<Result<_, _>>()
//                 .map_err(|e| format!("Error reading guides: {}", e))
//         })
// }

fn ensure_directory_exists(path: &str) -> Result<(), String> {
    if !Path::new(path).exists() {
        fs::create_dir_all(path).map_err(|e| format!("Failed to create {}: {}", path, e))?;
    }
    Ok(())
}

fn format_duration(d: Duration) -> String {
    let secs = d.as_secs_f64();
    if secs < 60.0 {
        format!("{:.3} s", secs)
    } else {
        let m = (secs / 60.0).floor();
        let s = secs % 60.0;
        format!("{} m {:.3} s", m, s)
    }
}

fn build_chopoff_db(tool: &Chopoff, cfg: &Config) -> Option<BenchmarkResult> {
    if let Some(ref db_path) = cfg.chopoff_db_path {
        if !Path::new(db_path).exists() {
            println!("Building Chopoff DB at {}...", db_path);
            let dur = tool
                .build(
                    cfg.target_file.as_str(),
                    db_path,
                    cfg.dists[0],
                    &format!("human_dist{}", cfg.dists[0]),
                    cfg.threads,
                )
                .expect("Build failed");

            return Some(BenchmarkResult {
                tool_name: tool.name().into(),
                operation: "build".into(),
                duration: dur,
                distance: cfg.dists[0],
                threads: cfg.threads,
                guide_seq: String::new(),
            });
        }
    }
    None
}

fn run_chopoff_search(tool: &Chopoff, cfg: &Config) -> Option<BenchmarkResult> {
    if let Some(ref db_path) = cfg.chopoff_db_path {
        println!("Running Chopoff search for {}...", cfg.guides_file);
        match tool.run(
            cfg.dists[0], // For now we just do dist loop in main, but nicer here, fixme
            cfg.guides_file.as_str(),
            db_path,
            &format!("{}/chopoff_{}.txt", cfg.out_dir, cfg.guides_file),
            cfg.threads,
        ) {
            Ok(dur) => Some(BenchmarkResult {
                tool_name: tool.name().into(),
                operation: "search".into(),
                duration: dur,
                distance: cfg.dists[0],
                threads: cfg.threads,
                guide_seq: cfg.guides_file.clone(),
            }),
            Err(e) => {
                eprintln!("Chopoff error: {}", e);
                None
            }
        }
    } else {
        eprintln!("No DB path for Chopoff, skipping");
        None
    }
}

fn run_other_tool(tool: &dyn Tool, cfg: &Config) -> Option<BenchmarkResult> {
    println!("Running {} for {}...", tool.name(), cfg.guides_file);
    match tool.run(
        cfg.dists[0],
        cfg.guides_file.as_str(),
        cfg.target_file.as_str(),
        &format!("{}/{}.txt", cfg.out_dir, tool.name().to_lowercase()),
        cfg.threads,
    ) {
        Ok(dur) => Some(BenchmarkResult {
            tool_name: tool.name().into(),
            operation: "run".into(),
            duration: dur,
            distance: cfg.dists[0],
            threads: cfg.threads,
            guide_seq: cfg.guides_file.clone(),
        }),
        Err(e) => {
            eprintln!("{} error: {}", tool.name(), e);
            None
        }
    }
}

fn run_benchmark(tools: &[ToolVariant], cfg: &Config) -> Vec<BenchmarkResult> {
    ensure_directory_exists(&cfg.out_dir).expect("Can't create out dir");
    let mut results = Vec::new();

    for tool in tools {
        match tool {
            ToolVariant::Chop(chopoff) => {
                if let Some(r) = build_chopoff_db(chopoff, cfg) {
                    results.push(r);
                }
                //for seq in guide_seqs {
                if let Some(r) = run_chopoff_search(chopoff, cfg) {
                    results.push(r);
                }
                //}
            }

            ToolVariant::Sassy(sassy) => {
                if let Some(r) = run_other_tool(sassy, cfg) {
                    results.push(r);
                }
            }

            ToolVariant::Swo(swo) => {
                if let Some(r) = run_other_tool(swo, cfg) {
                    results.push(r);
                }
            }
        }
    }
    results
}

fn print_results(results: &[BenchmarkResult]) {
    println!("\n=== Results ===");
    println!("Tool	Op	Dist	Th	Seq	Time");
    for r in results {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            r.tool_name,
            r.operation,
            r.distance,
            r.threads,
            r.guide_seq,
            format_duration(r.duration)
        );
    }
}

pub fn run(config_path: &str) {
    let cfg = load_config(config_path);
    let tools = vec![
        ToolVariant::Sassy(SassyTool::new(&cfg.sassy_path)),
        ToolVariant::Swo(Swofinder::new(&cfg.swofinder_path)),
        ToolVariant::Chop(Chopoff::new(&cfg.chopoff_path)),
    ];
    let mut all_results = Vec::new();
    for &dist in &cfg.dists {
        // Generate a unique DB path for each distance
        let db_path = format!(
            "{}_dist{}.db",
            cfg.chopoff_db_path.as_deref().unwrap_or("chopoff_db"),
            dist
        );
        let bench_cfg = Config {
            target_file: cfg.target_file.clone(),
            out_dir: cfg.out_dir.clone(),
            dists: vec![dist],
            threads: cfg.threads,
            chopoff_db_path: Some(db_path),
            sassy_path: cfg.sassy_path.clone(),
            swofinder_path: cfg.swofinder_path.clone(),
            chopoff_path: cfg.chopoff_path.clone(),
            guides_file: cfg.guides_file.clone(),
        };
        let results = run_benchmark(&tools, &bench_cfg);
        all_results.extend(results);
    }
    print_results(&all_results);
}
