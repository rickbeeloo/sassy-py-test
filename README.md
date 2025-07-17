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


### Search

---
**Search single pattern (--pattern)**


To search the pattern `ATGAGCA` in the fasta file `text.fasta` allowing up to `1` edit:
```bash 
cargo run -r -- sassy search --pattern "ATGAGCA" --alphabet dna -k 1 text.fasta
```
This will print the output to `stdout`, if you want to save it to a file use `--output matches.txt`. 
For alphabets see [alphabets section](#alphabets).




---
**Search with multi Fasta (--pattern-fasta)**


If you have more than one pattern to search, you can use `--pattern-fasta` instead of `--pattern`:
```bash 
cargo run -r -- sassy search --pattern-fasta patterns.fasta --alphabet dna -k 1 text.fasta
```
---

### Off-target (CRISPR)
To search a list of sgRNAs flanked by a PAM sequence. 

```bash 
cargo run -r -- sassy crispr --guide guides.txt --k 1 text.fasta
```
*Note* to stick to common format the input for `--guide` is a .txt file, not a fasta file, with a 
guide per line.
If you want to limit the matches with `N` characters, you can use `--max-n-frac`, and if you 
do allow edits in the PAM sequence you can use the `--allow-pam-edits` flag.


## Alphabets
Three alphabets are supported:
- *ASCII*: only equal characters match.
- *DNA*: Only `ACTG` and `actg` characters are expected, and treated case-insensitively.
- *IUPAC*: On top of the DNA characters, also supports `NYR` and furter
  characters (again, case insensitive), so that `A` matches `N`. See [here](https://www.bioinformatics.org/sms/iupac.html) for full table of IUPAC codes. 

## Evals
For the evals see [evals/README.md](evals/README.md).

## Python bindings
For python bindings see [python/README.md](python/README.md).
