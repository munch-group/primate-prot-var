#!/usr/bin/env python3
"""Fetch all human protein-coding genes and cache them."""

from primate_aa_variants import EnsemblClient, get_protein_coding_genes

print("Fetching all human protein-coding genes...")
print("This will take approximately 10-15 minutes.")
print()

client = EnsemblClient()
genes = get_protein_coding_genes(client, "human", cache_file=".cache/human_genes.pkl")

print()
print("="*60)
print("SUMMARY")
print("="*60)
print(f"Total genes fetched: {len(genes)}")
print(f"Expected: ~19,000-20,000")
print()

# Count by chromosome
from collections import Counter
chr_counts = Counter(g["chromosome"] for g in genes)
print("Genes per chromosome:")
for chrom in [str(i) for i in range(1, 23)] + ["X", "Y"]:
    count = chr_counts.get(chrom, 0)
    print(f"  Chr {chrom:2s}: {count:5d} genes")

print()
if len(genes) >= 19000:
    print("✓ SUCCESS: Found expected number of genes!")
else:
    print(f"⚠ WARNING: Only found {len(genes)} genes (expected ~19,000)")
