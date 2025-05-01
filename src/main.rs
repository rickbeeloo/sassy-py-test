use std::path::PathBuf;

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

    /// When set, use the extended ACTG+NYR... alphabet.
    #[arg(long, value_enum)]
    alphabet: Alphabet,
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

    let query = args.query.as_bytes();

    let mut reader = needletail::parse_fastx_file(args.path).unwrap();

    while let Some(record) = reader.next() {
        let record = record.unwrap();
        let matches = match args.alphabet {
            Alphabet::Ascii => search::<Ascii<true>>(query, &record.seq(), args.k),
            Alphabet::Dna => search::<Dna>(query, &record.seq(), args.k),
            Alphabet::Iupac => search::<Iupac>(query, &record.seq(), args.k),
        };
        let id = str::from_utf8(record.id()).unwrap();
        for m in matches {
            let cost = m.0;
            let start = m.1.first().unwrap().0;
            let end = m.1.last().unwrap().0;
            println!("{id}\t{cost}\t{start}\t{end}");
        }
    }
}
