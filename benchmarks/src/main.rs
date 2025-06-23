use clap::{Parser, Subcommand};

mod edlib_bench;
use edlib_bench::runner as edlib_runner;

mod crispr_bench;
use crispr_bench::runner as crispr_runner;

mod overhang;
use overhang::runner as overhang_runner;

mod profiles;
use profiles::runner as profiles_runner;

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

    /// Run the overhang benchmark
    Overhang {
        /// Path to the overhang config TOML file
        #[arg(long)]
        config: String,
    },

    /// Run the profiles benchmark
    Profiles {
        /// Path to the profiles config TOML file
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
        Commands::Overhang { config } => {
            println!("Running overhang benchmark");
            overhang_runner::run(&config);
        }
        Commands::Profiles { config } => {
            println!("Running profiles benchmark");
            profiles_runner::run(&config);
        }
    }
}
