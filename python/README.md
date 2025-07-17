# simd-sassy Python Bindings

This directory contains the Python interface for the simd-sassy sequence search library, powered by Rust and PyO3.

## Installation

You can install the package using [maturin](https://github.com/PyO3/maturin):

```bash
maturin develop --features python
```

Or install from PyPI (once published):

```bash
pip install simd-sassy
```

## Usage

A simple usage is as follows:

``` python
query = b"ATCGATCG"
text = b"GGGGATCGATCGTTTT"
# alphabet: ascii, dna, uipac
searcher = sassy.PySearcher("dna")
matches = searcher.search(query, text, k=1)
for i, match in enumerate(matches):
    print(f"Match {i+1}:")
    print(f"    Start: {match.text_start}")
    print(f"    End: {match.text_end}")
    print(f"    Cost: {match.cost}")
    print(f"    Strand: {match.strand}")
    print(f"    CIGAR: {match.cigar}")
```

This finds 3 matches:

``` text
Match 1:
    Start: 4
    End: 12
    Cost: 0
    Strand: +
    CIGAR: 8=
Match 2:
    Start: 6
    End: 14
    Cost: 1
    Strand: -
    CIGAR: 6=X=
Match 3:
    Start: 2
    End: 10
    Cost: 1
    Strand: -
    CIGAR: X7=
```

Further options are `sassy.PySearcher(alpha=0.5)` to allow overhang alignments,
and `sassy.PySearcher("dna", rc=False)` to disable reverse complements for DNA
or IUPAC strings.

See [example.py](sassy/example.py) for a larger example.
