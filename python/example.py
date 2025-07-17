#!/usr/bin/env python3
"""
Example usage of the sassy Python bindings.
"""

import sassy


# Example 1: Simple DNA search
print("=== DNA Search Example ===")
query = b"ATCGATCG"
text = b"GGGGATCGATCGTTTT"

# Search with DNA alphabet
matches = sassy.PySearcher("dna", alpha=0.5).search(query, text, k=0)

print(f"Query: {query.decode()}")
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
searcher = sassy.PySearcher("dna")

query2 = b"GCTAGCTA"
text2 = b"AAAAAGCTAGCTAAAAA"

matches2 = searcher.search(query2, text2, k=1)

print(f"Query: {query2.decode()}")
print(f"Text:  {text2.decode()}")
print(f"Found {len(matches2)} matches with k=1:")

for i, match in enumerate(matches2):
    print(f"  Match {i+1}: cost={match.cost}, strand={match.strand}")

# Example 3: ASCII search
print("\n=== ASCII Search Example ===")
query3 = b"hello"
text3 = b"world hello there"

matches3 = sassy.PySearcher("ascii").search(query3, text3, k=0)

print(f"Query: {query3.decode()}")
print(f"Text:  {text3.decode()}")
print(f"Found {len(matches3)} matches:")

for i, match in enumerate(matches3):
    print(
        f"  Match {i+1}: start={match.text_start}, end={match.text_end}, cost={match.cost}"
    )
