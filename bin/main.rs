mod crispr;
mod query;
mod search;

use clap::Parser;
use {
    crispr::{CrisprArgs, crispr},
    query::{QueryArgs, query},
    search::{SearchArgs, search},
};

#[derive(clap::Parser)]
#[command(author, version, about)]
enum Args {
    /// Default search behavior
    Search(SearchArgs),
    /// CRISPR-specific search with PAM and edit-free region
    Crispr(CrisprArgs),
    /// Search multiple queries from a FASTA file against a target FASTA file
    Query(QueryArgs),
}

fn main() {
    let args = Args::parse();
    env_logger::init();

    match args {
        Args::Search(search_args) => search(&mut search_args.clone()),
        Args::Crispr(crispr_args) => crispr(crispr_args),
        Args::Query(query_args) => query(&mut query_args.clone()),
    }
}
