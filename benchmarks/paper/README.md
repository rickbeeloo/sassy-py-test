

## Running the bench 
First make sure to create the "benchmarks" executable:

`cargo build --release -p benchmarks`

This creates the `benchmarks` executable in `/target/release/`.
This is what the python scripts use to run the benchmarks. 


### Tool comparison (0 & 1 match)

#### Generating the data
Should be used from within **this folder** otherwise change `bench_cmd = ""` in `run_tool_bench.py` to the correct folder
```bash
python3 run_tool_bench.py
```
This creates a folder called `benchmarks` which has the tool configurations and several `results_*.csv` files which have the bench results.

#### Creating throughput (0-match) plot
```bash
python3 plot_tool_bench.py
```
Creates `figs` folder if it does not exist yet and writes the plot to `tool_comparison_high_res.svg` 

#### Throughput stats
To get GB/s for the different cut-offs, and the speed up compared to Edlib use 
```bash
python3 throughput_stats.py results_*.csv
```

#### Creating trace cost plot 
```bash
python3 plot_trace_bench.py
```

### CIRSPR off-traget
We used the benchmark of the [ChopOff paper](https://www.biorxiv.org/content/10.1101/2025.01.06.603201v1.full.pdf) from [their gitlab](https://git.app.uib.no/valenlab/chopoff-benchmark/-/tree/master?ref_type=heads) with [our fork here](https://github.com/rickbeeloo/sassy-crispr-bench):

This workflows uses conda/mamba so you would need to install that if you don't have it already.

```
mamba env create --name sassy-benchmark --file environment.yaml
mamba activate sassy-benchmark
snakemake --cores 16
```
Then the output of the tools is in `out_dir` and the timings in `summary.txt`.

Modifications we made to the Chopoff benchmark code:
- **Force serial execution** as the code did not use `threads:` in each rule it would execute multiple tools simultaneously that all use the maximum 
CPU's thereby competing with each other. 
- **Time indexes**, the Chopoff index construction was not timed, we added timings for the construction time for all edit cut-offs.
- **Added Sassy**
- **Remove older tools** that were shown to perform worse (CRISPRITz, Cas-OFFinder)
- Used a genome without N characters as tools might deal with that different ([chm13](https://github.com/marbl/CHM13))

