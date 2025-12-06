

## Gnomad common (>10%) variants hg38

Retrieve all common missense aa variants (>10%)

  python scripts/primate_aa_variants.py --method gnomad --output results/common_missense

produces

  results/common_missense_human_gnomad.csv

## Gnomad loss-of-function variants hg38

Retrieve predicted loss of function mutations for all human protein coding genes

  python scripts/gnomad-optimized-pipeline.py --method tsv --version v2.1.1 --output-dir results/gnomad_lof

produces:

  results/gnomad_lof/gnomad_lof_all_genes_20251206_003756.csv
  results/gnomad_lof/gnomad_lof_all_genes_20251206_003756.parquet
  results/gnomad_lof/gnomad_lof_summary_20251206_003756.csv


## AlphaMissense missense variants hg38

Download AlphaMissense data and create hdf5 file searchable on the 'uniprot_id' column:

    cd data
    wget https://zenodo.org/records/8208688/files/AlphaMissense_aa_substitutions.tsv.gz?download=1
    mv AlphaMissense_aa_substitutions.tsv.gz?download=1 AlphaMissense_aa_substitutions.tsv.gz
    gzip -d --stdout AlphaMissense_aa_substitutions.tsv.gz | grep -v '#' > alpha_missense_hg38.tsv
    python -c 'import pandas ; pandas.read_csv("alpha_missense_hg38.tsv", sep="\t").to_parquet("alpha_missense_hg38.parquet");'
    rm -f alpha_missense_hg38.tsv




