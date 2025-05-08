use std::fs::File;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::{
    io::{BufWriter, Write},
    path::PathBuf,
    sync::Mutex,
};

use crate::Match;
use crate::{Strand, profiles::Iupac, search_maybe_rc};
use pa_types::CigarOp;

#[derive(clap::Parser)]
pub struct CrisprArgs {
    /// Guide sequence to search for (with PAM)
    #[arg(long, short)]
    guide: String,

    /// Report matches up to (and including) this distance threshold (excluding PAM).
    #[arg(long, short)]
    k: usize,

    /// Fasta file to search. May be gzipped.
    /// #[arg(long)]
    path: PathBuf,

    /// Require the first N bases of the guide to be exact matches.
    #[arg(long)]
    exact_prefix: Option<usize>,

    /// Require the last N bases of the guide to be exact matches.
    #[arg(long)]
    exact_suffix: Option<usize>,

    /// Whether to include matches of the reverse-complement string.
    #[arg(long, short)]
    rc: bool,

    /// Number of threads to use. All CPUs by default.
    #[arg(short = 'j', long)]
    threads: Option<usize>,

    /// Output file path.
    #[arg(short = 'o', long)]
    output: PathBuf,
}

fn check_edit_free(m: &Match, target: isize) -> bool {
    // We assume PAMs are always at the ends so we can either check if the first X are
    // Match or the last X
    let is_negative = target < 0;
    let to_check = match (m.strand, is_negative) {
        (Strand::Fwd, false) | (Strand::Rc, true) => m.cigar.ops.first().unwrap(),
        (Strand::Fwd, true) | (Strand::Rc, false) => m.cigar.ops.last().unwrap(),
    };
    // Check if all are matches and that we have at least abs(target) matches
    to_check.op == CigarOp::Match && to_check.cnt >= target.abs() as i32
}

pub fn crispr(args: CrisprArgs) {
    let guide_sequence = args.guide.as_bytes();

    // Create a buffered writer for output
    let file = File::create(&args.output).expect("Failed to create output file");
    let writer = Mutex::new(BufWriter::new(file));

    let edit_free: Option<isize> = match (args.exact_prefix, args.exact_suffix) {
        (Some(_), Some(_)) => {
            eprintln!("[crispr] Error: cannot specify both exact prefix and suffix");
            std::process::exit(1);
        }
        (Some(prefix), None) => Some(prefix as isize),
        (None, Some(suffix)) => Some(-(suffix as isize)),
        (None, None) => None,
    };
    let edit_free_value = match edit_free {
        Some(v) => {
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
            v
        }
        None => {
            println!("[PAM] Edits are allowed");
            0
        }
    };

    let reader = Mutex::new(needletail::parse_fastx_file(args.path).unwrap());

    let total_found = Arc::new(AtomicUsize::new(0));
    let edits_in_pam = Arc::new(AtomicUsize::new(0));

    // Launch as many threads as there are cores.
    let num_threads = args.threads.unwrap_or_else(num_cpus::get);
    std::thread::scope(|scope| {
        for _ in 0..num_threads {
            scope.spawn(|| {
                while let Ok(mut guard) = reader.lock()
                    && let Some(record) = guard.next()
                {
                    let record = record.unwrap();
                    let id = String::from_utf8(record.id().to_vec()).unwrap();
                    let text = &record.seq().into_owned();

                    // Using IUPAC by default as most PAMs contain IUPAC characters.
                    let matches = search_maybe_rc::<Iupac>(guide_sequence, text, args.k, args.rc);

                    total_found.fetch_add(matches.len(), Ordering::Relaxed);

                    let mut writer_guard = writer.lock().unwrap();
                    for m in matches {
                        let pass = if edit_free.is_some() {
                            check_edit_free(&m, edit_free_value)
                        } else {
                            true
                        };
                        if pass {
                            let cost = m.cost;
                            let start = m.start.1 as usize;
                            let end = m.end.1 as usize;
                            let slice = str::from_utf8(&text[start..end]).unwrap();
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
                }
            });
        }
    });

    println!("\n[crispr] Summary");
    println!(
        "  Total targets found:   {}",
        total_found.load(Ordering::Relaxed)
    );
    if edit_free.is_some() {
        println!(
            "  Discarded (edit-free): {}",
            edits_in_pam.load(Ordering::Relaxed)
        );
        println!(
            "  Total targets passed:  {}",
            total_found.load(Ordering::Relaxed) - edits_in_pam.load(Ordering::Relaxed)
        );
    } else {
        println!("  Discarded (edit-free): 0");
        println!(
            "  Total targets passed:  {}",
            total_found.load(Ordering::Relaxed)
        );
    }
}
