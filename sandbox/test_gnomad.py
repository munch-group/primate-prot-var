#!/usr/bin/env python3
"""Test gnomAD API query."""

import requests
import json

gnomad_url = "https://gnomad.broadinstitute.org/api/"

print("="*60)
print("TEST 1: Simple query (working)")
print("="*60)

simple_query = """
query {
    gene(gene_symbol: "TP53", reference_genome: GRCh38) {
        gene_id
        symbol
    }
}
"""

try:
    response = requests.post(
        gnomad_url,
        json={"query": simple_query},
        headers={"Content-Type": "application/json"},
        timeout=30
    )
    print(f"Status: {response.status_code}")
    print(f"Response: {response.text[:200]}")
except Exception as e:
    print(f"Error: {e}")

print("\n" + "="*60)
print("TEST 2: Query with variables (from script)")
print("="*60)

# Query from the script
variant_query = """
query GeneVariants($geneSymbol: String!, $dataset: DatasetId!) {
    gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
        gene_id
        symbol
        variants(dataset: $dataset) {
            variant_id
            pos
            ref
            alt
            exome {
                af
                ac
                an
            }
            genome {
                af
                ac
                an
            }
            transcript_consequence {
                gene_symbol
                transcript_id
                consequence_terms
                amino_acids
                codons
                hgvsp
                protein_position
                is_canonical
            }
        }
    }
}
"""

try:
    response = requests.post(
        gnomad_url,
        json={
            "query": variant_query,
            "variables": {
                "geneSymbol": "TP53",
                "dataset": "gnomad_r4"
            }
        },
        headers={"Content-Type": "application/json"},
        timeout=60
    )
    print(f"Status: {response.status_code}")
    print(f"Response: {response.text[:1000]}")

    if response.status_code == 200:
        data = response.json()
        if "errors" in data:
            print("\nERRORS:")
            print(json.dumps(data["errors"], indent=2))
        else:
            gene = data.get("data", {}).get("gene", {})
            variants = gene.get("variants", [])
            print(f"\nSuccess! Found {len(variants)} total variants for TP53")
            if variants:
                print(f"First variant: {variants[0].get('variant_id')}")

except Exception as e:
    print(f"Error: {e}")
