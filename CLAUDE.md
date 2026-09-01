# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Research analysis repo (not an installable package) for comparative analysis of
protein-coding variation across primates: human population data (gnomAD, VEP,
AlphaMissense) integrated with primate ortholog alignments and baboon population
genomics. Runs on the GenomeDK cluster (SLURM). Compute nodes have no internet —
do downloads on the login node; keep large data on project storage, not `$HOME`.

The `pyproject.toml` packaging metadata (`src/` layout, console script) and the
Quarto book config (`_quarto.yml`, `docs/`, "projectname"/"Joanna Doh") are
munch-group project-template boilerplate — there is no `src/` package and no test
suite. The real content is `scripts/`, `vep_data/`, `long-read-lof/`, and
`notebooks/`.

## Environments (pixi)

Two independent pixi environments; always run tools through `pixi run` from the
directory owning the environment (linux-64 only):

- **Repo root** — defined in `pyproject.toml` `[tool.pixi.*]`: pandas, biopython,
  sgkit, pysam, MACSE, ensembl-vep, duckdb, etc. Used by `scripts/` and notebooks.
- **`vep_data/`** — own `pixi.toml`: gwf 2.x, ensembl-vep 115, bcftools, samtools,
  cyvcf2, pyarrow, duckdb. Used only by the VEP annotation workflow.

## The three strands

### 1. `scripts/` — gnomAD / AlphaMissense / alignments (repo root)

Standalone scripts, each documented in the root `README.md` (read it for per-script
details). The main entry points:

```bash
# common (>10% MAF) missense variants via gnomAD GraphQL
pixi run python scripts/primate_aa_variants.py --method gnomad --output results/common_missense

# predicted LoF variants for all human protein-coding genes
pixi run python scripts/gnomad-optimized-pipeline.py --method tsv --version v2.1.1 --output-dir results/gnomad_lof

# codon-aware primate ortholog alignments (MACSE) for one gene
pixi run python scripts/get_alignments.py -t hominidae TTLL10
```

AlphaMissense data lives in `data/` as `alpha_missense_hg38.parquet` (built from the
Zenodo TSV; recipe in root README) plus an HDF5 keyed by UniProt ID
(`build_alpha_missense_hdf5.py` → `annotate_var.py` for lookups).
`extract_variants.py` works on baboon (papAnu4) Zarr stores produced by
`sgkit_workflow.py`.

### 2. `vep_data/` — VEP annotation pipeline (gwf)

Self-contained gwf workflow (`vep_data/workflow.py`) that annotates the full
Ensembl human variation VCFs (release 115, GRCh38) with VEP and emits a
Hive-partitioned Parquet dataset for DuckDB queries.

Per-chromosome pipeline, all intermediates under `vep_data/steps/`:
filter (`bcftools view -e 'ALT="<.>"'`, drops structural variants VEP can't parse)
→ split into 10 Mb chunks → VEP per chunk (offline cache in `steps/vep_cache/`,
`--fork 4`) → `bcftools concat` per chromosome → parse CSQ with cyvcf2 into
`steps/parquet/chrom=chrN/part-0.parquet` (a `.done` sentinel per partition marks
completion for gwf).

```bash
cd vep_data
pixi run gwf status            # SLURM backend; state in .gwf/
pixi run gwf run               # submit all targets (or: gwf run RunVEP_chr1_1)
pixi run gwf logs <target>
```

Input VCFs (`steps/vcf/homo_sapiens-chrN.vcf.gz`) and the VEP cache tarball are
fetched manually by rsync from Ensembl FTP — commands in `vep_data/README.md`.
SLURM account is set at the top of `workflow.py` (`ACCOUNT`).

Legacy versions of this pipeline exist and are superseded by the gwf workflow:
`scripts/vep_sbatch.sh` + `scripts/vcf2parquet.py` + `README_vep.md` (sbatch
array over chromosomes, cache in `vep_data/.vep/`), and a copy of the sbatch
script inside `vep_data/`. Prefer `vep_data/workflow.py` for new work.

Root-level `vep_data.py` is a separate demo of VEP query patterns (Ensembl REST
for ≤200 variants, FTP bulk + Parquet/DuckDB/tabix); root `workflow.py` is a gwf
*tutorial example* rendered into the Quarto book, not a real pipeline.

### 3. `long-read-lof/` — primate gene-loss intersections

Has its own detailed `CLAUDE.md`; read it before working in that directory.

## Gotchas

- Ensembl human VCFs name chromosomes *without* the `chr` prefix; the workflow
  strips it for `bcftools -r`, while output Parquet partitions and filenames use
  `chr`-prefixed names.
- VEP requires the reference FASTA to be **bgzip**-compressed (plain gzip fails
  for HGVS); hence the gunzip→bgzip→`samtools faidx` dance in `PrepareFasta`.
- Query the Parquet dataset with DuckDB hive partitioning, e.g.
  `read_parquet('steps/parquet/**/*.parquet', hive_partitioning=true)`; rows are
  sorted by `pos` so range predicates use row-group stats.
- Large data (`data/`, `steps/`, `gnomad_cache/`, VCFs, caches) is untracked —
  don't add it to git; regenerate from the recipes in the READMEs.
