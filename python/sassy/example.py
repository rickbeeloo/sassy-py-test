#!/usr/bin/env python3
"""
Example usage of the sassy Python bindings.
"""

import sassy

# Random text of length 1000000000
import random
import time

sassy.features()

print("gen")
n = 100000
m = 20
k = 1
text = bytes(random.choices(b"ACGT", k=n))
pattern = bytes(random.choices(b"ACGT", k=m))
# time instant now

searcher = sassy.Searcher("dna", rc=False)
for _ in range(10):
    start = time.time()
    matches = searcher.search(pattern, text, k=k)
    print(f"GB/s: {n / (time.time()-start) / 10**9}")

# Example 1: Simple DNA search
print("=== DNA Search Example ===")
pattern = b"ATCGATCG"
text = b"GGGGATCGATCGTTTT"

# Search with DNA alphabet
matches = sassy.Searcher("dna", alpha=0.5).search(pattern, text, k=0)

print(f"Pattern: {pattern.decode()}")
print(f"Text:  {text.decode()}")
print(f"Found {len(matches)} matches:")

for i, match in enumerate(matches):
    print(f"  Match {i+1}:")
    print(f"    Start: {match.text_start}")
    print(f"    End: {match.text_end}")
    print(f"    Cost: {match.cost}")
    print(f"    Strand: {match.strand}")
    print(f"    CIGAR: {match.cigar}")

# Example 2: Reverse complement search
print("\n=== Reverse Complement Example ===")
searcher = sassy.Searcher("dna")

pattern2 = b"GCTAGCTA"
text2 = b"AAAAAGCTAGCTAAAAA"

matches2 = searcher.search(pattern2, text2, k=1)

print(f"Pattern: {pattern2.decode()}")
print(f"Text:  {text2.decode()}")
print(f"Found {len(matches2)} matches with k=1:")

for i, match in enumerate(matches2):
    print(f"  Match {i+1}: cost={match.cost}, strand={match.strand}")

# Example 3: ASCII search
print("\n=== ASCII Search Example ===")
pattern3 = b"hello"
text3 = b"world hello there"

matches3 = sassy.Searcher("ascii").search(pattern3, text3, k=0)

print(f"Pattern: {pattern3.decode()}")
print(f"Text:  {text3.decode()}")
print(f"Found {len(matches3)} matches:")

for i, match in enumerate(matches3):
    print(
        f"  Match {i+1}: start={match.text_start}, end={match.text_end}, cost={match.cost}"
    )
