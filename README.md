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
> cargo run -r --  --help

Usage: sassy <COMMAND>

Commands:
  search  Default search behavior
  crispr  CRISPR-specific search with PAM and edit-free region
  query   Search multiple queries from a FASTA file against a target FASTA file
  help    Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

```

## Example

### Search
To search a single pattern (i.e. `ACTGCTACTGTACA`) in a fasta file

To find matches of a pattern with up to 1 edit:

``` rust
cargo run -r -- search --alphabet dna ACTGCTACTGTACA 1 hg.fa
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

### Query
When searching a multi-fasta against another (multi) fasta file. For example to search a list of protospacers in a contig collection.

To find all spacer matches with up to 3 edits

```rust
cargo run -r --  query  --alphabet iupac spacers.fasta 3 contigs.fasta
```


## Alphabets
Three alphabets are supported:
- *ASCII*: only equal characters match.
- *DNA*: Only `ACTG` and `actg` characters are expected, and treated case-insensitively.
- *IUPAC*: On top of the DNA characters, also supports `NYR` and furter
  characters (again, case insensitive), so that `A` matches `N`.
