#![feature(let_chains)]
use std::{path::PathBuf, sync::Mutex};

use clap::Parser;
use sassy::{
    profiles::{Ascii, Dna, Iupac},
    search,
};

#[derive(clap::Parser)]
struct Args {
    /// Pattern to search for.
    query: String,
    /// Report matches up to this threshold.
    k: usize,
    /// Fasta file to search. May be gzipped.
    path: PathBuf,

    /// The alphabet to use. DNA=ACTG, or IUPAC=ACTG+NYR...
    #[arg(long, value_enum)]
    alphabet: Alphabet,

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
                        Alphabet::Ascii => search::<Ascii<true>>(query, text, args.k),
                        Alphabet::Dna => search::<Dna>(query, text, args.k),
                        Alphabet::Iupac => search::<Iupac>(query, text, args.k),
                    };

                    let _write_lock = write_lock.lock().unwrap();
                    for m in matches {
                        let cost = m.0;
                        let start = m.1.first().unwrap().1;
                        let end = m.1.last().unwrap().1;
                        let slice = str::from_utf8(&text[start..end]).unwrap();
                        println!("{id}\t{cost}\t{start}\t{end}\t{slice}");
                    }
                }
            });
        }
    });
}
