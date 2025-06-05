use sassy::{
    profiles::{Ascii, Dna, Iupac},
    search::{Searcher, Strand},
};
use std::{path::PathBuf, sync::Mutex};

//FIXME: either adjust CLI to match cirispr (named arguments) or other way around

#[derive(clap::Parser)]
pub struct SearchArgs {
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
pub enum Alphabet {
    Ascii,
    #[default]
    Dna,
    Iupac,
}

pub fn search(args: SearchArgs) {
    let query = args.query.as_bytes();
    let reader = Mutex::new(needletail::parse_fastx_file(args.path).unwrap());
    let write_lock = Mutex::new(());

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

                    let matches = match args.alphabet {
                        Alphabet::Ascii => {
                            Searcher::<Ascii, false>::new_fwd().search(query, &text, args.k)
                        }
                        Alphabet::Dna => {
                            Searcher::<Dna, false>::new(args.rc).search(query, &text, args.k)
                        }
                        Alphabet::Iupac => {
                            Searcher::<Iupac, false>::new(args.rc).search(query, &text, args.k)
                        }
                    };

                    let _write_lock = write_lock.lock().unwrap();
                    for m in matches {
                        let cost = m.cost;

                        // We have to adjust the start and end based on reverse complement
                        // as we reverse the text these should be adjusted based on text length
                        let (start, end) = if m.strand == Strand::Rc {
                            (
                                text.len() - m.end.1 as usize,
                                text.len() - m.start.1 as usize,
                            )
                        } else {
                            (m.start.1 as usize, m.end.1 as usize)
                        };

                        let slice = &text[start..end];
                        let slice_str = String::from_utf8_lossy(slice);
                        let cigar = m.cigar.to_string();
                        let strand = match m.strand {
                            Strand::Fwd => "+",
                            Strand::Rc => "-",
                        };
                        println!("{id}\t{cost}\t{strand}\t{start}\t{end}\t{slice_str}\t{cigar}");
                    }
                }
            });
        }
    });
}
