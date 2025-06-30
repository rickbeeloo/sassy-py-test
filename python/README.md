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

```python
import sassy
searcher = sassy.PySearcher('dna', False)
results = searcher.search(b'ACGT', b'ACGTTT', 1)
for m in results:
    print(m.start, m.end, m.cost, m.strand, m.cigar)
```

See the main project README for more details and advanced usage. 