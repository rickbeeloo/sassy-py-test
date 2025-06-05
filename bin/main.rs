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
    /// Default search behavior
    Search(SearchArgs),
    /// CRISPR-specific search with PAM and edit-free region
    Crispr(CrisprArgs),
}

fn main() {
    let args = Args::parse();
    env_logger::init();

    match args {
        Args::Search(search_args) => search(search_args),
        Args::Crispr(crispr_args) => crispr(crispr_args),
    }
}
