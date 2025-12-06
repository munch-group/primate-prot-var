#!/usr/bin/env python
"""Test gnomAD API with known good genes"""

import requests
import json
import pandas as pd

def query_gnomad_gene(gene_symbol):
    """Query gnomAD API for a specific gene"""
    url = "https://gnomad.broadinstitute.org/api"

    query = """
    query GeneVariants($geneSymbol: String!, $datasetId: DatasetId!) {
      gene(gene_symbol: $geneSymbol, reference_genome: GRCh37) {
        gene_id
        symbol
        variants(dataset: $datasetId) {
          variant_id
          consequence
          lof
          exome {
            ac
            af
          }
        }
      }
    }
    """

    variables = {
        "geneSymbol": gene_symbol,
        "datasetId": "gnomad_r2_1"
    }

    response = requests.post(
        url,
        json={'query': query, 'variables': variables},
        timeout=30
    )

    if response.status_code == 200:
        data = response.json()
        if 'errors' in data:
            print(f"❌ {gene_symbol}: {data['errors']}")
            return None

        if 'data' in data and data['data']['gene']:
            gene_data = data['data']['gene']
            n_variants = len(gene_data.get('variants', []))
            lof_variants = [v for v in gene_data['variants'] if v.get('lof') == 'HC']
            print(f"✅ {gene_symbol}: {n_variants} total variants, {len(lof_variants)} LoF variants")
            return gene_data
        else:
            print(f"❌ {gene_symbol}: Gene not found")
            return None
    else:
        print(f"❌ API error: {response.status_code}")
        return None

# Test with known cancer genes
test_genes = ['BRCA1', 'BRCA2', 'TP53', 'PTEN', 'ATM']

print("Testing gnomAD API with known genes:")
print("=" * 60)

for gene in test_genes:
    query_gnomad_gene(gene)
