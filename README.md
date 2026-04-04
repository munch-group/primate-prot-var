
# Primate Protein Variation

Comparative analysis of protein-coding variation across primates. Integrates human population variant data (gnomAD, AlphaMissense) with primate ortholog alignments to study amino acid variation in an evolutionary context.

## Overview of scripts

### `primate_aa_variants.py`
Compiles common (>10% MAF) missense amino acid variants in human, chimpanzee, and gorilla. Supports three data retrieval methods: Ensembl REST API (`--method ensembl`), Ensembl BioMart (`--method biomart`), and gnomAD GraphQL API (`--method gnomad`). The gnomAD method retrieves comprehensive annotations including population frequencies, in silico predictor scores (CADD, REVEL, AlphaMissense), and LoF flags for human genes. Produces per-species CSV files.

### `gnomad-optimized-pipeline.py`
Scalable pipeline for retrieving predicted loss-of-function (LoF) variants across all human protein-coding genes from gnomAD. Offers three processing backends: gnomAD GraphQL API (`--method api`), pre-built TSV constraint files (`--method tsv`), and Parquet files (`--method parquet`). Outputs full variant tables (CSV + Parquet) and a per-gene summary with variant counts and mean allele frequencies.

### `gnomad-lof-retrieval.py`
Collection of approaches for retrieving LoF variants with population-specific allele frequencies from gnomAD. Implements four methods: gnomAD GraphQL API, `pynoma` package, `gnomad-db` SQLite database, and direct VCF streaming via `pysam`/tabix. Includes a batch wrapper to process gene lists and save results to CSV.

### `gnomad-all-lof-download.py`
Utilities for downloading and processing the full set of gnomAD LoF variants. Includes downloading pre-processed constraint files, streaming VCF files chromosome-by-chromosome via `pysam`, setting up the gnomAD SQLite database (~100 GB), and downloading pre-filtered high-confidence pLoF files.

### `gnomad-all-genes-lof.py`
Full-scale pipeline for extracting LoF variants from gnomAD VCF files for all protein-coding genes with gene names, CDS positions, mutation types, and population frequencies. Supports processing via Hail, parallel VCF processing across chromosomes, or sequential processing. Produces a per-variant table and a per-gene summary.

### `get_alignments.py`
Fetches ortholog CDS sequences for a given human gene across a specified primate taxon using the Ensembl REST API. Converts gene symbols to Ensembl/UniProt IDs, retrieves ortholog coding sequences, appends human CDS variants, and runs codon-aware alignment with MACSE. Produces protein and CDS alignment FASTA files and per-sequence alignment statistics.

    python scripts/get_alignments.py -t hominidae TTLL10

produces `TTLL10_ENSG00000162571_Q6ZVT0_hominidae_protein.fa`, `*_cds.fa`, `*_stats.csv`.

### `extract_variants.py`
Extracts amino acid variants from baboon (papAnu4) population genomic data stored in Zarr format (via sgkit). Parses a GTF annotation to identify CDS coordinates, reconstructs per-sample haplotype CDS sequences, translates them, and compares to the reference protein to identify nonsynonymous changes and in-frame stop codons. Outputs a CSV of per-sample protein variants.

### `annotate_var.py`
Looks up AlphaMissense pathogenicity scores for a given gene from a pre-built HDF5 store. Takes a gene/UniProt ID and the HDF5 file path as arguments and prints the matching records.

### `build_alpha_missense_hdf5.py`
Converts the AlphaMissense amino acid substitutions TSV (gzipped) into an HDF5 file indexed by UniProt ID for fast per-protein lookups. Each UniProt entry is stored under a hierarchical key (`<first_letter>/<uniprot_id>`).

### `alphamissense_downloader.py`
Downloads AlphaMissense pathogenicity data from Zenodo and extracts scores for a specific protein (default: TTLL10 / Q6ZVM7). Includes multiple download strategies (direct streaming, `zenodo_get`, local file extraction, HegedLab API) and a basic pathogenicity analysis that categorizes variants as likely benign, ambiguous, or likely pathogenic.

### `sgkit_workflow.py`
GWF workflow script that converts merged BCF files from the baboon diversity project into per-chromosome Zarr stores using `bcftools` and `vcf2zarr`. Reads metadata to identify VCF files and chromosome regions, then maps conversion jobs across chromosomes and subspecies.

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




