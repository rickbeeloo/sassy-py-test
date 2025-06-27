use sassy::{
    profiles::{Ascii, Dna, Iupac, Profile},
    search::{Searcher, Strand},
};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
    sync::{Arc, Mutex},
};

#[derive(clap::Parser, Clone)]
pub struct QueryArgs {
    /// FASTA file containing queries to search for. May be gzipped.
    queries: PathBuf,

    /// Report matches up to (and including) this distance threshold.
    k: usize,

    /// FASTA file to search against. May be gzipped.
    target: PathBuf,

    /// The alphabet to use. DNA=ACTG, or IUPAC=ACTG+NYR...
    #[arg(long, value_enum)]
    alphabet: Alphabet,

    /// Disable reverse complement search
    #[arg(long)]
    no_rc: bool,

    /// Number of threads to use. All CPUs by default.
    #[arg(short = 'j', long)]
    threads: Option<usize>,

    /// Output file for results. Defaults to stdout.
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Overhang cost
    #[arg(long)]
    overhang_cost: Option<f32>,
}

#[derive(clap::ValueEnum, Default, Clone, PartialEq)]
pub enum Alphabet {
    Ascii,
    #[default]
    Dna,
    Iupac,
}

// As we cant build "generic" searchers since they implement different traits, we need to
// build a macro to generate the code for each alphabet.The least ugly option I guess
macro_rules! run_query_worker {
    ($profile:ty, $make_searcher:expr, $args:ident, $query_reader:ident, $write_lock:ident, $target_path:ident) => {{
        let mut searcher = $make_searcher;
        loop {
            // We lock the query reader, pull a seq, and instantly unlock so other threads can pull queries
            let query_work = {
                let mut reader = $query_reader.lock().unwrap();
                reader.next().map(|r| {
                    let record = r.unwrap();
                    let query_id = String::from_utf8(record.id().to_vec()).unwrap();
                    let query = record.seq().into_owned();
                    (query_id, query)
                })
            };
            if query_work.is_none() {
                break;
            }
            let (query_id, query) = query_work.unwrap();
            let mut target_reader = needletail::parse_fastx_file(&$target_path).unwrap();
            while let Some(target_record) = target_reader.next() {
                let target_record = target_record.unwrap();
                let target_id = String::from_utf8(target_record.id().to_vec()).unwrap();
                let text = target_record.seq().into_owned();
                let matches = searcher.search(&query, &text, $args.k);

                let mut results = Vec::new();
                for m in matches {
                    let cost = m.cost;
                    let start = m.start.1 as usize;
                    let end = m.end.1 as usize;
                    let slice = &text[start..end];
                    let slice_str = if m.strand == Strand::Rc {
                        String::from_utf8_lossy(&<$profile as Profile>::reverse_complement(slice)).into_owned()
                    } else {
                        String::from_utf8_lossy(slice).into_owned()
                    };
                    let cigar = m.cigar.to_string();
                    let strand = match m.strand {
                        Strand::Fwd => "+",
                        Strand::Rc => "-",
                    };
                    results.push(format!("{query_id}\t{target_id}\t{cost}\t{strand}\t{start}\t{end}\t{slice_str}\t{cigar}"));
                }

                if !results.is_empty() {
                    let mut writer = $write_lock.lock().unwrap();
                    for result in results {
                        writeln!(writer, "{}", result).unwrap();
                    }
                }
            }
        }
    }};
}

pub fn query(args: &mut QueryArgs) {
    // Auto-disable rc search for ASCII profile as it does not make sense, let the user know
    if args.alphabet == Alphabet::Ascii && !args.no_rc {
        eprintln!(
            "WARNING: Reverse complement search is not supported for ASCII profile, disabling"
        );
        args.no_rc = true;
    }
    let overhang_cost = args.overhang_cost;
    if overhang_cost.is_some() {
        println!(
            "Using overhang cost: {} for alignments extending beyond sequence bounds",
            overhang_cost.unwrap()
        );
    }

    let num_threads = args.threads.unwrap_or_else(num_cpus::get);

    // Create output writer, shared and needs to be locked to be thread safe
    let output_writer = if let Some(output_path) = &args.output {
        Box::new(BufWriter::new(File::create(output_path).unwrap())) as Box<dyn Write + Send>
    } else {
        Box::new(std::io::stdout()) as Box<dyn Write + Send>
    };

    let write_lock = Arc::new(Mutex::new(output_writer));

    // Create a shared query reader, all threads pull queries from this reader
    let query_reader = Arc::new(Mutex::new(
        needletail::parse_fastx_file(&args.queries).unwrap(),
    ));

    std::thread::scope(|scope| {
        for _ in 0..num_threads {
            let query_reader = Arc::clone(&query_reader);
            let write_lock = Arc::clone(&write_lock);
            let thread_args = args.clone();
            let target_path = args.target.clone();

            scope.spawn(move || match thread_args.alphabet {
                Alphabet::Ascii => run_query_worker!(
                    Ascii,
                    Searcher::<Ascii>::new(false, overhang_cost),
                    thread_args,
                    query_reader,
                    write_lock,
                    target_path
                ),
                Alphabet::Dna => run_query_worker!(
                    Dna,
                    Searcher::<Dna>::new(!thread_args.no_rc, overhang_cost),
                    thread_args,
                    query_reader,
                    write_lock,
                    target_path
                ),
                Alphabet::Iupac => run_query_worker!(
                    Iupac,
                    Searcher::<Iupac>::new(!thread_args.no_rc, overhang_cost),
                    thread_args,
                    query_reader,
                    write_lock,
                    target_path
                ),
            });
        }
    });
}

#[allow(unused)]
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
    fn end_to_end_query() {
        eprintln!("WARNING: Run this test with -- --nocapture to see the output");

        // Create queries file with multiple queries
        let mut queries_file = std::fs::File::create("data/queries.fasta").unwrap();
        let query1 = b"TAGCTAGAC";
        let query2 = rc(query1);
        writeln!(queries_file, ">query1").unwrap();
        writeln!(queries_file, "{}", String::from_utf8_lossy(query1)).unwrap();
        writeln!(queries_file, ">query2").unwrap();
        writeln!(queries_file, "{}", String::from_utf8_lossy(&query2)).unwrap();

        // Create target file with sequences containing matches
        let mut target_file = std::fs::File::create("data/target.fasta").unwrap();

        // First target sequence
        let dna1 = random_dna_string(1000);
        let mut text1 = dna1.clone();
        text1.splice(10..10, query1.iter().copied());
        text1.splice(50..50, rc(query1).iter().copied());
        writeln!(target_file, ">target1").unwrap();
        writeln!(target_file, "{}", String::from_utf8_lossy(&text1)).unwrap();

        // Second target sequence
        let dna2 = random_dna_string(1000);
        let mut text2 = dna2.clone();
        text2.splice(20..20, query2.iter().copied());
        text2.splice(80..80, rc(&query2).iter().copied());
        writeln!(target_file, ">target2").unwrap();
        writeln!(target_file, "{}", String::from_utf8_lossy(&text2)).unwrap();

        let mut args: QueryArgs = QueryArgs {
            queries: PathBuf::from("data/queries.fasta"),
            k: 0,
            target: PathBuf::from("data/target.fasta"),
            alphabet: Alphabet::Dna,
            no_rc: true,
            threads: Some(1),
            output: None,
            overhang_cost: None,
        };

        println!("Query without RC");
        query(&mut args.clone());

        println!("Query with RC");
        args.no_rc = false;
        query(&mut args);
    }
}
