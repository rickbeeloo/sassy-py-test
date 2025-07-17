use needletail::parse_fastx_file;
use sassy::rec_iter::{Query, RecordIterator};
use sassy::{
    profiles::{Ascii, Dna, Iupac, Profile},
    search::{Match, OwnedStaticText, SearchAble, Searcher, Strand},
};
use std::fs::File;
use std::sync::Arc;
use std::{
    io::{BufWriter, Write},
    path::PathBuf,
    sync::Mutex,
};

#[derive(clap::Parser, Clone)]
pub struct SearchArgs {
    // Named args
    /// Pattern to search for.
    #[arg(short = 'p', long, conflicts_with = "pattern_fasta")]
    pattern: Option<String>,

    /// Pattern file
    #[arg(short = 'f', long, conflicts_with = "pattern")]
    pattern_fasta: Option<String>,

    /// Report matches up to (and including) this distance threshold.
    #[arg(short, long)]
    k: usize,

    // Flags
    /// The alphabet to use. DNA=ACTG, or IUPAC=ACTG+NYR...
    #[arg(long, short = 'a', value_enum)]
    alphabet: Alphabet,

    /// Disable reverse complement search
    #[arg(long)]
    no_rc: bool,

    // Optional
    /// Number of threads to use. All CPUs by default.
    #[arg(short = 'j', long)]
    threads: Option<usize>,

    /// Output file
    #[arg(short = 'o', long)]
    output: Option<String>,

    // Positional
    /// Fasta file to search. May be gzipped.
    path: PathBuf,
}

#[derive(clap::ValueEnum, Default, Clone, PartialEq)]
pub enum Alphabet {
    Ascii,
    #[default]
    Dna,
    Iupac,
}

fn get_queries(args: &SearchArgs) -> Vec<Query> {
    if let Some(p) = &args.pattern {
        // Single inline pattern; give it a dummy id
        vec![Query {
            id: "query".to_string(),
            seq: p.as_bytes().to_vec(),
        }]
    } else if let Some(pattern_fasta) = &args.pattern_fasta {
        // Pull all sequences and ids from the FASTA file
        let mut reader = parse_fastx_file(pattern_fasta).expect("valid path/file");
        let mut queries = Vec::new();
        while let Some(record) = reader.next() {
            let seqrec = record.expect("invalid record");
            let id = String::from_utf8(seqrec.id().to_vec()).unwrap_or_else(|_| "query".into());
            queries.push(Query {
                id,
                seq: seqrec.seq().to_vec(),
            });
        }
        queries
    } else {
        unreachable!("No pattern or pattern_fasta provided - clap should handle");
    }
}

fn get_output_writer(args: &SearchArgs) -> Box<dyn Write + Send> {
    if let Some(output_path) = &args.output {
        Box::new(BufWriter::new(File::create(output_path).unwrap())) as Box<dyn Write + Send>
    } else {
        Box::new(std::io::stdout()) as Box<dyn Write + Send>
    }
}

fn ensure_valid_profile(args: &mut SearchArgs) {
    if args.alphabet == Alphabet::Ascii && !args.no_rc {
        eprintln!(
            "WARNING: Reverse complement search is not supported for ASCII profile, disabling"
        );
        args.no_rc = true;
    }
}

fn as_output_line(
    m: &Match,
    pat_id: &str,
    text_id: &str,
    owned_static_text: &OwnedStaticText,
    alphabet: &Alphabet,
) -> String {
    let cost = m.cost;
    let start = m.start.1 as usize;
    let end = m.end.1 as usize;
    let slice = &owned_static_text.text[start..end];

    // If we match reverse complement, reverse complement the slice to make it easier to read
    let slice_str = if m.strand == Strand::Rc {
        match alphabet {
            Alphabet::Dna => {
                String::from_utf8_lossy(&<Dna as Profile>::reverse_complement(slice)).into_owned()
            }
            Alphabet::Iupac => {
                String::from_utf8_lossy(&<Iupac as Profile>::reverse_complement(slice)).into_owned()
            }
            Alphabet::Ascii => unreachable!("no rc for ascii"), // Guarded against above
        }
    } else {
        String::from_utf8_lossy(slice).into_owned()
    };

    let cigar = m.cigar.to_string();
    let strand = match m.strand {
        Strand::Fwd => "+",
        Strand::Rc => "-",
    };
    format!("{pat_id}\t{text_id}\t{cost}\t{strand}\t{start}\t{end}\t{slice_str}\t{cigar}\n")
}

// To create a fixed type, cost of match should be neglible
#[derive(Clone)]
enum SearchWrapper {
    Ascii(Searcher<Ascii>),
    Dna(Searcher<Dna>),
    Iupac(Searcher<Iupac>),
}

impl SearchWrapper {
    fn search<I: SearchAble>(&mut self, query: &[u8], input: &I, k: usize) -> Vec<Match> {
        match self {
            SearchWrapper::Ascii(s) => s.search(query, input, k),
            SearchWrapper::Dna(s) => s.search(query, input, k),
            SearchWrapper::Iupac(s) => s.search(query, input, k),
        }
    }
}

pub fn search(args: &mut SearchArgs) {
    // Get queries based on `pattern` or `pattern_fasta`
    let queries = get_queries(args);
    assert!(!queries.is_empty(), "No query sequences found");

    // Create output writer, stdout by default (or user provided), and share it safely across threads
    let output_writer = Arc::new(Mutex::new(get_output_writer(args)));

    // Auto-disable rc search for ASCII profile as it does not make sense, let the user know
    ensure_valid_profile(args);

    let k = args.k;
    let rc_enabled = !args.no_rc;
    let alphabet = args.alphabet.clone();

    let num_threads = args.threads.unwrap_or_else(num_cpus::get);
    let iter = Arc::new(RecordIterator::new(&args.path, queries, None));
    std::thread::scope(|s| {
        for _ in 0..num_threads {
            let it = iter.clone();
            let out = output_writer.clone();
            let alphabet = alphabet.clone();

            s.spawn(move || {
                // Each thread has own searcher here
                let mut searcher: SearchWrapper = match alphabet {
                    Alphabet::Ascii => SearchWrapper::Ascii(Searcher::<Ascii>::new(false, None)),
                    Alphabet::Dna => SearchWrapper::Dna(Searcher::<Dna>::new(rc_enabled, None)),
                    Alphabet::Iupac => {
                        SearchWrapper::Iupac(Searcher::<Iupac>::new(rc_enabled, None))
                    }
                };

                while let Some(batch) = it.next_batch() {
                    for item in batch {
                        let pat_id = &item.query.id;
                        let pat_seq = &item.query.seq;
                        let rec = item.record.as_ref();
                        let thread_id = thread_id::get();
                        eprintln!("Thread {thread_id} q: {pat_id}, against text id: {}", rec.0,);
                        let matches = searcher.search(pat_seq, &rec.1, k);
                        for m in matches {
                            let line = as_output_line(&m, pat_id, &rec.0, &rec.1, &alphabet);
                            out.lock().unwrap().write_all(line.as_bytes()).unwrap();
                        }
                    }
                }
            });
        }
    });
}

#[cfg(test)]
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
        let mut file = std::fs::File::create("test.fasta").unwrap();
        writeln!(file, ">test").unwrap();
        writeln!(file, "{}", String::from_utf8_lossy(&text)).unwrap();

        let mut args: SearchArgs = SearchArgs {
            pattern: Some(String::from_utf8(query.to_vec()).unwrap()),
            k: 0,
            path: PathBuf::from("test.fasta"),
            alphabet: Alphabet::Dna,
            no_rc: true,
            threads: Some(1),
            pattern_fasta: None,
            output: None,
            all: false,
        };

        // FIXME: capture output and assert or write to file for easy check
        // anyway for now run with -- --nocapture and check for 10,50,100
        println!("Search without RC");
        search(&mut args.clone());

        println!("Search with RC");
        args.no_rc = false;
        search(&mut args);
    }
}
