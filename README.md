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
  search  Search a single sequence or multi-fasta in a multi-fasta text
  crispr  CRISPR-specific search with PAM and edit-free region
  help    Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

```

## Examples

#### Search single pattern (--pattern)
To search the pattern `ATGAGCA` in the fasta file `text.fasta` allowing up to `1` edit:
```bash 
cargo run -r -- sassy search --pattern "ATGAGCA" --alphabet dna -k 1 text.fasta
```
This will print the output to `stdout`, if you want to save it to a file use `--output matches.txt`. 
For alphabets see [alphabets section](#alphabets).

#### Search with multi Fasta (--pattern-fasta)
If you have more than one pattern to search, you can use `--pattern-fasta` instead of `--pattern`:
```bash 
cargo run -r -- sassy search --pattern-fasta patterns.fasta --alphabet dna -k 1 text.fasta
```

### CRISPR off-target 
To search a list of sgRNAs flanked by a PAM sequence. 

```bash 
cargo run -r -- sassy crispr --guide guides.txt --k 1 text.fasta
```
*Note* to stick to common format the input for `--guide` is a .txt file, not a fasta file, with a 
guide per line.
If you want to limit the matches with `N` charactres, you can use `--max-n-frac`, and if you 
do allow edits in the PAM sequence you can use the `--allow-pam-edits` flag.


### Output
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

## Evals
For the evals see [evals/README.md](evals/README.md).

## Python bindings
For python bindings see [python/README.md](python/README.md).
