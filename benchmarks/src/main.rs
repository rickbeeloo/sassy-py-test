use clap::{Parser, Subcommand};

mod edlib_bench;
use edlib_bench::runner as edlib_runner;

mod crispr_bench;
use crispr_bench::runner as crispr_runner;

#[derive(Parser)]
#[command(author, version, about)]
struct Args {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run the edlib grid benchmark
    Edlib {
        /// Path to the grid config TOML file
        #[arg(long)]
        config: String,
    },
    /// Run the CRISPR benchmark
    Crispr {
        /// Path to the CRISPR config TOML file
        #[arg(long)]
        config: String,
    },
}

fn main() {
    let args = Args::parse();
    match args.command {
        Commands::Edlib { config } => {
            println!("Running edlib grid");
            edlib_runner::run(&config);
        }
        Commands::Crispr { config } => {
            println!("Running CRISPR benchmark");
            crispr_runner::run(&config);
        }
    }
}
