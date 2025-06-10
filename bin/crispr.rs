use sassy::{
    profiles::Iupac, profiles::Profile, search::Searcher, search::StaticText, search::Strand,
};
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

    /// Whether to include matches of the reverse-complement string.
    #[arg(long, short = 'r')]
    rc: bool,

    /// Allow edits in PAM sequence 
    #[arg(long, short = 'e')]
    allow_pam_edits: bool,

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

fn check_n_frac(max_n_frac: f32, match_slice: &[u8]) -> bool {
    let n_count = match_slice
        .iter()
        .filter(|c| (**c & 0xDF) == b'N') // Convert to uppercase check against N
        .count() as f32;
    let n_frac = n_count / match_slice.len() as f32;
    n_frac <= max_n_frac
}



fn print_and_check_params(args: &CrisprArgs, guide_sequences: &[Vec<u8>]) -> (String, f32) {
    let max_n_frac = args.max_n_frac.unwrap_or(1.0); // Allow all to be N by default

    // Check if n frac is within valid range
    if !(0.0..=1.0).contains(&max_n_frac) {
        eprintln!("[N-chars] Error: max_n_frac must be between 0 and 100");
        std::process::exit(1);
    }

    // If no guide sequences throw error
    if guide_sequences.is_empty() {
        eprintln!("[PAM] Error: No guide sequences provided, please check your input file (one guide sequence per line)");
        std::process::exit(1);
    }

    // We have at least one guide sequence, extract the PAM sequence
    let pam = if !guide_sequences.is_empty() {
        let guide = &guide_sequences[0];
        let pam = &guide[guide.len() - 3..];
        println!("[PAM] Sequence: [{}]", String::from_utf8_lossy(pam));
        println!("[PAM] If the above PAM is incorrect, please make sure that the guide sequence ENDs with the PAM-sequence, i.e. XXXXXGGN (not it's reverse complement)");
        pam
    } else {
        unreachable!("No guide sequences provided");
    };

    // If we have multiple guide sequences, ensure they have the same PAM sequence.
    // not per se a requirement for the code to work but now we define that as a fixed PAM in the closure below
    if guide_sequences.len() > 1 {
        for guide_sequence in guide_sequences {
            let guide_pam = &guide_sequence[guide_sequence.len() - 3..];
            if pam != guide_pam {
                eprintln!("[PAM] One of the guide sequences has a PAM different than the provided PAM");
                eprintln!("[PAM] provided PAM {}, detected PAM {}", String::from_utf8_lossy(pam), String::from_utf8_lossy(guide_pam));
                std::process::exit(1);
            }
        }
    }

    println!("[PAM] PAM used to filter: {}", String::from_utf8_lossy(pam));
    println!("[PAM] Edits in PAM are allowed: {}", args.allow_pam_edits);
    println!("[N-chars] Allowing up to {}% N characters", max_n_frac * 100.0);

    (String::from_utf8_lossy(pam).into_owned(), max_n_frac)
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

pub fn matching_seq<P: Profile>(seq1: &[u8], seq2: &[u8]) -> bool {
    for (c1, c2) in seq1.iter().zip(seq2.iter()) {
        if !P::is_match(*c1, *c2) {
            return false;
        }
    }
    true
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

        let (pam, max_n_frac) = print_and_check_params(&args, &guide_sequences);
        let pam = pam.as_bytes();
        let pam_compl = Iupac::complement(pam);
        let pam_compl = pam_compl.as_slice();

        let total_found = Arc::new(AtomicUsize::new(0));

        let num_threads = args.threads.unwrap_or_else(num_cpus::get);
        println!("[Threads] Using {num_threads} threads");

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
                        let mut searcher = Searcher::<Iupac>::new_rc();

                        // Only used in loop if no edits are allowed, otherwise falls back to default report all search
                        let filter_fn = 
                        move |_q: &[u8], text_up_to_end: &[u8], strand: Strand| {
                            let pam_slice = &text_up_to_end[text_up_to_end.len() - 3..];
                            if strand == Strand::Fwd {
                                matching_seq::<Iupac>(pam_slice, pam)
                            } else {
                                matching_seq::<Iupac>(pam_slice, pam_compl)
                            }
                        };
                        

                        // Search for each guide sequence
                        guide_sequences.iter().for_each(|guide_sequence| {
                        
                            let matches = if !args.allow_pam_edits {
                                searcher.search_with_fn(guide_sequence, &static_text, args.k, true, filter_fn)
                            } else {
                                searcher.search_all(guide_sequence, &static_text, args.k)
                            };

                            total_found.fetch_add(matches.len(), Ordering::Relaxed);

                            let mut writer_guard = writer.lock().unwrap();

                            for m in matches {
                                let start = m.start.1 as usize;
                                let end = m.end.1 as usize;
                                let slice = &text[start..end];

                                // Check if satisfies user max N cut off
                                let n_ok = if max_n_frac < 100.0 {
                                    check_n_frac(max_n_frac, slice)
                                } else {
                                    true
                                };

                                if !n_ok {
                                    // println!("N-frac too high: {}", String::from_utf8_lossy(slice));
                                    continue;
                                }

                                total_found.fetch_add(1, Ordering::Relaxed);

                                // If the match we found is in reverse complement, we also rc the matched text
                                // to make it easier to spot the PAM in the output file
                                let matched_slice = if m.strand == Strand::Rc {
                                    let rc = <Iupac as Profile>::reverse_complement(slice);
                                    String::from_utf8_lossy(&rc).into_owned()
                                } else {
                                    String::from_utf8_lossy(slice).into_owned()
                                };
                                let cost = m.cost;
                                let cigar = m.cigar.to_string();
                                let strand = match m.strand {
                                    Strand::Fwd => "+",
                                    Strand::Rc => "-",
                                };
                                writeln!(
                                    writer_guard,
                                    "{id}\t{cost}\t{strand}\t{start}\t{end}\t{matched_slice}\t{cigar}"
                                )
                                .unwrap();
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
        println!("  Time taken: {:?}", start.elapsed());
    }
}


