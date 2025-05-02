#![feature(let_chains)]
use std::{path::PathBuf, sync::Mutex};

use clap::Parser;
use sassy::{
    Strand,
    profiles::{Ascii, Dna, Iupac},
    search_maybe_rc,
};

#[derive(clap::Parser)]
struct Args {
    /// Pattern to search for.
    query: String,
    /// Report matches up to (and including) this distance threshold.
    k: usize,
    /// Fasta file to search. May be gzipped.
    path: PathBuf,

    /// The alphabet to use. DNA=ACTG, or IUPAC=ACTG+NYR...
    #[arg(long, value_enum)]
    alphabet: Alphabet,

    /// Whether to include matches of the reverse-complement string.
    #[arg(long)]
    rc: bool,

    /// Number of threads to use. All CPUs by default.
    #[arg(short = 'j', long)]
    threads: Option<usize>,
}

#[derive(clap::ValueEnum, Default, Clone)]
enum Alphabet {
    Ascii,
    #[default]
    Dna,
    Iupac,
}

fn main() {
    let args = Args::parse();

    env_logger::init();

    let query = args.query.as_bytes();

    let reader = Mutex::new(needletail::parse_fastx_file(args.path).unwrap());
    let write_lock = Mutex::new(());

    // Launch as many threads as there are cores.
    let num_threads = args.threads.unwrap_or_else(|| num_cpus::get());
    std::thread::scope(|scope| {
        for _ in 0..num_threads {
            scope.spawn(|| {
                while let Ok(mut guard) = reader.lock()
                    && let Some(record) = guard.next()
                {
                    let record = record.unwrap();
                    let id = String::from_utf8(record.id().to_vec()).unwrap();
                    let text = &record.seq().into_owned();

                    // Release the lock.
                    drop(record);
                    drop(guard);

                    // eprintln!("Searching {id}");
                    let matches = match args.alphabet {
                        Alphabet::Ascii => {
                            search_maybe_rc::<Ascii<true>>(query, text, args.k, args.rc)
                        }
                        Alphabet::Dna => search_maybe_rc::<Dna>(query, text, args.k, args.rc),
                        Alphabet::Iupac => search_maybe_rc::<Iupac>(query, text, args.k, args.rc),
                    };

                    let _write_lock = write_lock.lock().unwrap();
                    for m in matches {
                        let cost = m.cost;
                        let start = m.start.1 as usize;
                        let end = m.end.1 as usize;
                        let slice = str::from_utf8(&text[start..end]).unwrap();
                        let cigar = m.cigar.to_string();
                        let strand = match m.strand {
                            Strand::Fwd => "+",
                            Strand::Rc => "-",
                        };
                        println!("{id}\t{cost}\t{strand}\t{start}\t{end}\t{slice}\t{cigar}");
                    }
                }
            });
        }
    });
}
