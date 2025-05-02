# Sassy: SIMD Approximate String Searching

Sassy is a library and tool for approximately searching short patterns in texts,
the problem that goes by many names:
- approximate string matching,
- pattern searching,
- fuzzy searching.

The motivating application is searching short (~20bp) DNA fragments in a human
genome (3GB), but it also works well for longer patterns up to ~1Kbp, and
shorter texts.

## Usage

```
> cargo run -r -- --help

Usage: sassy [OPTIONS] --alphabet <ALPHABET> <QUERY> <K> <PATH>

Arguments:
  <QUERY>  Pattern to search for
  <K>      Report matches up to (and including) this distance threshold
  <PATH>   Fasta file to search. May be gzipped

Options:
      --alphabet <ALPHABET>  The alphabet to use. DNA=ACTG, or IUPAC=ACTG+NYR... [possible values: ascii, dna, iupac]
      --rc                   Whether to include matches of the reverse-complement string
  -j, --threads <THREADS>    Number of threads to use. All CPUs by default
  -h, --help                 Print help
```

You can also first build with `cargo build -r` and then do
`./target/release/sassy <args>`.

## Example

To find matches of a pattern with up to 1 edit:

``` rust
cargo run -r -- --alphabet dna --rc ACTGCTACTGTACA 1 hg.fa
```

Output is written as tab-separated values to stdout, containing the sequence id,
distance, strand, start and end position, matched substring, and cigar string.
(For matches to the reverse-complement strand, the query is reversed and matches
are reported in the forward direction.)

```
chr1	2	+	74462851	74462868	ATCGGTGTCATCAATAA	 =D11=X4=
chr1	2	+	97381917	97381934	ACTCGGTGTCCTCATAA	 10=X2=D4=
chr1	2	+	196285921	196285938	actggtgtcatcggtaa	 3=D10=X3=
chr1	2	-	199471583	199471601	TTATCGATGACACTGAAT	 13=X2=X=
chr1	2	-	229999068	229999085	ATATCGATGACACCAGT	 X13=D3=
chr1	2	-	231082126	231082144	ttatcaatgacaacgagt	 5=X6=X5=
```

## Alphabets
Three alphabets are supported:
- *ASCII*: only equal characters match.
- *DNA*: Only `ACTG` and `actg` characters are expected, and treated case-insensitively.
- *IUPAC*: On top of the DNA characters, also supports `NYR` and furter
  characters (again, case insensitive), so that `A` matches `N`.
