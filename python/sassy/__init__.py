# sassy Python package 

# Import the Rust bindings
from . import sassy as _sassy

# Re-export the main classes and functions
PySearcher = _sassy.PySearcher
PyMatch = _sassy.PyMatch
search_sequence = _sassy.search_sequence
search_sequences = _sassy.search_sequences

__all__ = ['PySearcher', 'PyMatch', 'search_sequence', 'search_sequences'] 