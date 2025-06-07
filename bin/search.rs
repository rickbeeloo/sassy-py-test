use sassy::{
    profiles::{Ascii, Dna, Iupac},
    search::{Searcher, Strand},
};
use std::{path::PathBuf, sync::Mutex};

//FIXME: either adjust CLI to match cirispr (named arguments) or other way around

#[derive(clap::Parser, Clone)]
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

    /// Disable reverse complement search
    #[arg(long)]
    no_rc: bool,

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
                            Searcher::<Ascii>::new_fwd().search(query, &text, args.k)
                        }
                        Alphabet::Dna => {
                            Searcher::<Dna>::new(!args.no_rc).search(query, &text, args.k)
                        }
                        Alphabet::Iupac => {
                            Searcher::<Iupac>::new(!args.no_rc).search(query, &text, args.k)
                        }
                    };

                    let _write_lock = write_lock.lock().unwrap();
                    for m in matches {
                        let cost = m.cost;
                        let start = m.start.1 as usize;
                        let end = m.end.1 as usize;
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

mod test {

    use super::*;
    use rand::Rng;
    use std::io::Write;

    fn random_dna_string(len: usize) -> Vec<u8> {
        let mut rng = rand::rng();
        (0..len).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
    }

    fn rc(dna: &[u8]) -> Vec<u8> {
        dna.iter()
            .rev()
            .map(|c| match c {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => panic!("Invalid DNA character"),
            })
            .collect()
    }

    #[test]
    fn end_to_end_search() {
        eprintln!("WARNING: Run this test with -- --nocapture to see the output");
        // Create file at data/test.fasta, with two fwd, one reverse complement match
        // insert at 10, 50, and 100
        let dna = random_dna_string(1000);
        let mut text = dna.clone();
        let query = b"TAGCTAGAC";
        text.splice(10..10, query.iter().copied());
        text.splice(50..50, query.iter().copied());
        text.splice(100..100, rc(query).iter().copied());
        // Write to file
        let mut file = std::fs::File::create("data/test.fasta").unwrap();
        writeln!(file, ">test").unwrap();
        writeln!(file, "{}", String::from_utf8_lossy(&text)).unwrap();

        let mut args: SearchArgs = SearchArgs {
            query: String::from_utf8(query.to_vec()).unwrap(),
            k: 0,
            path: PathBuf::from("data/test.fasta"),
            alphabet: Alphabet::Dna,
            no_rc: true,
            threads: Some(1),
        };

        // FIXME: capture output and assert or write to file for easy check
        // anyway for now run with -- --nocapture and check for 10,50,100
        println!("Search without RC");
        search(args.clone());

        println!("Search with RC");
        args.no_rc = false;
        search(args);
    }
}
