mod crispr;
mod search;

use clap::Parser;
use {
    crispr::{CrisprArgs, crispr},
    search::{SearchArgs, search},
};

#[derive(clap::Parser)]
#[command(author, version, about)]
enum Args {
    /// Search a single sequence or multi-fasta in a multi-fasta text
    Search(SearchArgs),
    /// CRISPR-specific search with PAM and edit-free region
    Crispr(CrisprArgs),
}

fn main() {
    let args = Args::parse();
    env_logger::init();

    match args {
        Args::Search(search_args) => search(&mut search_args.clone()),
        Args::Crispr(crispr_args) => crispr(crispr_args),
    }
}
