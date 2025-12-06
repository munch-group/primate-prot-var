#!/usr/bin/env python3
"""Check output columns and sample data."""

import pandas as pd

df = pd.read_csv('results/test_variants_human_gnomad.csv')

print(f'Columns ({len(df.columns)}):')
for col in df.columns:
    print(f'  {col}')

print()
print('Sample data (first row, transposed):')
print(df.iloc[0].to_string())
