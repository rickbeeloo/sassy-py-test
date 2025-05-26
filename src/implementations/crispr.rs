use crate::search::Match;
use crate::{profiles::Iupac, search::Searcher, search::StaticText, search::Strand};
use pa_types::CigarOp;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;
use std::{
    io::{BufWriter, Write},
    path::PathBuf,
    sync::Mutex,
};

#[derive(clap::Parser)]
pub struct CrisprArgs {
    /// Path to file with guide sequences (including PAM)
    #[arg(long, short = 'g')]
    guide: String,

    /// Report matches up to (and including) this distance threshold (excluding PAM).
    #[arg(long, short = 'k')]
    k: usize,

    /// Fasta file to search. May be gzipped.
    #[arg(long, short = 't')]
    target: PathBuf,

    /// Require the first N bases of the guide to be exact matches.
    #[arg(long, short = 'p')]
    exact_prefix: Option<usize>,

    /// Require the last N bases of the guide to be exact matches.
    #[arg(long, short = 's')]
    exact_suffix: Option<usize>,

    /// Whether to include matches of the reverse-complement string.
    #[arg(long, short = 'r')]
    rc: bool,

    /// Number of threads to use. All CPUs by default.
    #[arg(short = 'j', long)]
    threads: Option<usize>,

    /// Allow at most max_n_frac of N characters in the target sequence.
    #[arg(long, short = 'n')]
    max_n_frac: Option<f32>,

    /// Output file path.
    #[arg(short = 'o', long)]
    output: PathBuf,
}

fn check_edit_free(m: &Match, target: isize) -> bool {
    // We assume PAMs are always at the ends so we can either check if the first X are
    // Match or the last X
    let is_negative = target < 0;
    let to_check: &pa_types::CigarElem = match is_negative {
        true => m.cigar.ops.last().unwrap(),
        false => m.cigar.ops.first().unwrap(),
    };
    to_check.op == CigarOp::Match && to_check.cnt >= target.abs() as i32
}

fn check_n_frac(m: &Match, max_n_frac: f32, target_text: &[u8]) -> bool {
    let start = m.start.1 as usize;
    let end = m.end.1 as usize;
    let slice_len = end - start;
    let n_count = target_text[start..end]
        .iter()
        .filter(|c| (**c & 0xDF) == b'N') // Convert to uppercase check against N
        .count() as f32;
    let n_frac = n_count / slice_len as f32;
    n_frac <= max_n_frac
}

fn print_and_check_params(args: &CrisprArgs, guide_sequence: &[u8]) -> (Option<isize>, f32) {
    // Only allow one of prefix/suffix
    if let (Some(_), Some(_)) = (args.exact_prefix, args.exact_suffix) {
        eprintln!("[crispr] Error: cannot specify both exact prefix and suffix");
        std::process::exit(1);
    }

    let edit_free: Option<isize> = match (args.exact_prefix, args.exact_suffix) {
        (Some(prefix), None) => Some(prefix as isize),
        (None, Some(suffix)) => Some(-(suffix as isize)),
        (None, None) => None,
        _ => None,
    };

    let max_n_frac = args.max_n_frac.unwrap_or(100.0); // Allow all to be N by default

    // Print info
    if let Some(v) = edit_free {
        let guide_str = String::from_utf8_lossy(guide_sequence);
        if v < 0 {
            let pam_start = guide_sequence.len() - (-v) as usize;
            let prefix = &guide_str[..pam_start];
            let pam = &guide_str[pam_start..];
            println!("[PAM] Edit-free region in brackets: {}[{}]", prefix, pam);
        } else {
            let pam = &guide_str[..v as usize];
            let suffix = &guide_str[v as usize..];
            println!("[PAM] Edit-free region in brackets: [{}]{}", pam, suffix);
        }
    } else {
        println!("[PAM] Edits are allowed");
    }

    if args.max_n_frac.is_some() {
        println!(
            "[N-chars] Allowing up to {}% N characters",
            max_n_frac * 100.0
        );
    } else {
        println!("[N-chars] No N-character filtering");
    }

    (edit_free, max_n_frac)
}

fn pass(
    m: &Match,
    edit_free: Option<isize>,
    edit_free_value: isize,
    max_n_frac: f32,
    text: &[u8],
) -> bool {
    let pam_ok = if edit_free.is_some() {
        check_edit_free(m, edit_free_value)
    } else {
        true
    };
    let n_ok = if max_n_frac < 100.0 {
        check_n_frac(m, max_n_frac, text)
    } else {
        true
    };
    pam_ok && n_ok
}

pub fn read_guide_sequences(path: &str) -> Vec<Vec<u8>> {
    let file = File::open(path).expect("Failed to open guide file");
    let reader = BufReader::new(file);
    reader
        .lines()
        .map(|l| l.unwrap().as_bytes().to_vec())
        .collect::<Vec<_>>()
        .into_iter()
        .filter(|seq| !seq.is_empty())
        .collect()
}

pub fn crispr(args: CrisprArgs) {
    let guide_sequences = read_guide_sequences(&args.guide);
    println!("[GUIDES] Found {} guides", guide_sequences.len());

    if !guide_sequences.is_empty() {
        // Read the first record from the FASTA file for benchmarking
        println!("Creating output file with path: {}", args.output.display());
        let file = File::create(&args.output).expect("Failed to create output file");
        let writer = Mutex::new(BufWriter::new(file));
        let reader = Mutex::new(needletail::parse_fastx_file(args.target.clone()).unwrap());

        let (edit_free, max_n_frac) = print_and_check_params(&args, &guide_sequences[0]);
        let edit_free_value = edit_free.unwrap_or(0);

        let total_found = Arc::new(AtomicUsize::new(0));
        let edits_in_pam = Arc::new(AtomicUsize::new(0));

        let num_threads = args.threads.unwrap_or_else(num_cpus::get);
        println!("[Threads] Using {} threads", num_threads);

        let start = Instant::now();
        std::thread::scope(|scope| {
            for _ in 0..num_threads {
                scope.spawn(|| {
                    while let Ok(mut guard) = reader.lock()
                        && let Some(record) = guard.next()
                    {
                        // Get fasta record
                        let record = record.unwrap();
                        let id = String::from_utf8(record.id().to_vec()).unwrap();
                        let text = &record.seq().into_owned();

                        // Create static text by precomputing reverse
                        let static_text = StaticText::new(text);

                        // Searcher, IUPAC and always reverse complement
                        let mut searcher = Searcher::<Iupac, true, true>::new();

                        // Search for each guide sequence
                        guide_sequences.iter().for_each(|guide_sequence| {
                            let matches = searcher.search(guide_sequence, &static_text, args.k);

                            total_found.fetch_add(matches.len(), Ordering::Relaxed);

                            let mut writer_guard = writer.lock().unwrap();
                            for m in matches {
                                if pass(&m, edit_free, edit_free_value, max_n_frac, text) {
                                    let cost = m.cost;
                                    let start = m.start.1 as usize;
                                    let end = m.end.1 as usize;

                                    let slice = match m.strand {
                                        Strand::Fwd => {
                                            String::from_utf8_lossy(&text[start..end]).to_string()
                                        }
                                        Strand::Rc => String::from_utf8_lossy(
                                            text[text.len() - end..text.len() - start]
                                                .iter()
                                                .rev()
                                                .copied()
                                                .collect::<Vec<_>>()
                                                .as_slice(),
                                        )
                                        .to_string(),
                                    };
                                    let cigar = m.cigar.to_string();
                                    let strand = match m.strand {
                                        Strand::Fwd => "+",
                                        Strand::Rc => "-",
                                    };
                                    writeln!(
                                        writer_guard,
                                        "{id}\t{cost}\t{strand}\t{start}\t{end}\t{slice}\t{cigar}"
                                    )
                                    .unwrap();
                                } else {
                                    edits_in_pam.fetch_add(1, Ordering::Relaxed);
                                }
                            }
                        });
                    }
                });
            }
        });

        println!("\nSummary");
        println!(
            "  Total targets found:   {}",
            total_found.load(Ordering::Relaxed)
        );
        println!(
            "  Discarded (edits + N's): {}",
            edits_in_pam.load(Ordering::Relaxed)
        );
        println!(
            "  Total targets passed:  {}",
            total_found.load(Ordering::Relaxed) - edits_in_pam.load(Ordering::Relaxed)
        );
        println!("  Time taken: {:?}", start.elapsed());
    }
}
