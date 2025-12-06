#!/usr/bin/env python3
"""Test gene fetching for a single chromosome."""

import sys
sys.path.insert(0, '.')

from primate_aa_variants import EnsemblClient, get_chromosome_info, get_protein_coding_genes

# Test chromosome info
print("Testing chromosome info fetch...")
client = EnsemblClient()
chrom_info = get_chromosome_info(client, "human")
print(f"\nChromosomes found: {len(chrom_info)}")
for chrom, length in sorted(chrom_info.items(), key=lambda x: (len(x[0]), x[0])):
    print(f"  Chr {chrom}: {length:,} bp")

# Test gene fetching for chromosome 21 (smallest autosome)
print("\n" + "="*60)
print("Testing gene fetch for chromosome 21 only...")
print("="*60)

# Temporarily modify the function to only query chr21
import primate_aa_variants
original_func = primate_aa_variants.get_chromosome_info

def chr21_only(client, species):
    full_info = original_func(client, species)
    return {"21": full_info["21"]}

primate_aa_variants.get_chromosome_info = chr21_only

genes = get_protein_coding_genes(client, "human", cache_file=".cache/human_chr21_test.pkl")
print(f"\nTotal genes on chr21: {len(genes)}")
print("\nFirst 10 genes:")
for gene in genes[:10]:
    print(f"  {gene['gene_name']:15s} {gene['gene_id']:20s} {gene['start']:10,}-{gene['end']:10,}")
