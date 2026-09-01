# CLAUDE.md

Project context for an assistant working in this directory on GenomeDK.

## What this is

`intersect_primate_lof.py` cross-references a **query gene list** (e.g. the
meiosis/recombination/sex-chromosome set in `meiosis_repair_genes.py`) against
lineage-specific loss-of-function / gene-status data from long-read primate
genome resources, and writes per-branch intersections plus an UpSet-ready matrix.

Driving question: *which genes in a curated set are lost / disrupted on which
branches of the primate phylogeny, and how do those losses distribute by
chromosome (esp. chrX, for the meiotic-drive / sex-chromosome-conflict angle).*

Two input modes:

| Mode | Flag | Scope | Notes |
|---|---|---|---|
| Single workbook | `--table` | one curated supplement | Mao 2024 Data S2, Yoo 2025 Suppl. Table VIII.34, or any gene+lineage[+status] table |
| TOGA tree | `--toga-dir` | every gene × every assembly | genome-wide; the real "all genes, all species" path |

The paper supplements are curated (SV-disruption / gain-loss summaries). The
genome-wide per-gene × per-species intact/lost/missing/paralog matrix comes from
**TOGA** (Kirilenko et al. 2023). Use `--toga-dir` when the task is "all genes,
all species".

Two companion scripts extend this to humans and to a unified view:
- `annotate_sv_lof.py` — turns the **1000 Genomes long-read SV catalog**
  (1KG_ONT_VIENNA) into a human gene-LoF `--table` source (`data/human_1kgp_sv_lof.tsv`
  → `results/human_1kgp_sv/lof.*`). This is the **only human loss source** (TOGA/Mao
  use human as the reference). LoF tiers `Human_anyLoF/commonLoF/homozygousLoF` +
  superpops; **polymorphic**, not fixed. See `METHODS.md` §2F.
- `merge_loss_matrices.py` — outer-joins the per-source `loss_matrix`/`matrix` files
  into one source-prefixed gene × lineage matrix (`results/merged/`) plus a coverage
  mask. See `METHODS.md` §3 (Cross-source combination).

## Environment

GenomeDK, SLURM. The tool is pure Python + pandas; tiny and CPU/RAM-light (runs
on a login node or a short interactive `srun`). Only **downloads** need internet
(login / transfer node). Keep the downloaded TOGA tree and big supplements on
project storage (`/faststorage/project/<proj>/...`), not `$HOME`.

Deps: `python >=3.10`, `pandas`, `openpyxl` (for `.xlsx`). gz inputs are read
transparently (no manual gunzip).

```bash
pixi init . && pixi add "python>=3.10" pandas openpyxl
# or: mamba create -n lofx "python>=3.10" pandas openpyxl -y && mamba activate lofx
```

## The TOGA data (already downloaded here)

Pulled from the Senckenberg mirror with:

```bash
wget -e robots=off --recursive --no-parent \
  https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates
```

This grabs the **whole** tree including multi-GB codon/protein FASTAs. To refresh
just the status files next time, restrict to the one file the script needs:

```bash
wget -e robots=off -r -np -nH --cut-dirs=4 -A 'loss_summ_data.tsv.gz' \
  https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates/
```

Layout: `…/Primates/<Genus_species>__<Common_name>__<AssemblyID>/`. Important
facts the walker accounts for:

- **Everything is gzipped** (`*.tsv.gz`, `*.bed.gz`, `*.fa.gz`).
- **Two file-naming conventions coexist** in the same download:
  - set A: `loss_summ_data.tsv.gz`, `orthologsClassification.tsv.gz`,
    `geneAnnotation.bed.gz`/`.gtf.gz`, `geneInactivatingMutations.tsv.gz`,
    `processedPseudogeneAnnotation.bed.gz`, `codonAlignments[.allCESARexons].fa.gz`,
    `proteinAlignment(s)[.allCESARexons].fa.gz`.
  - set B (the `…T` / standard-TOGA dirs): `loss_summ_data.tsv.gz`,
    `orthology_classification.tsv.gz`, `inact_mut_data.txt.gz`,
    `query_annotation.bed.gz`, `codon.fasta.gz`, `prot.fasta.gz`, `<dir>.md5.txt`.
  - **`loss_summ_data.tsv.gz` is the only file present in every dir under both
    conventions** — the walker keys on it and ignores the rest.
- **Multiple assemblies per species** (e.g. `panTro6`, `HLpanTro7`, `HLpanTroT`
  for *Pan troglodytes*). `--toga-label species` collapses them (most-severe per
  gene); `--toga-label assembly` keeps each separate.
- Common-name field may be a placeholder (`Primates`, `-`); species names may be
  trinomial (`Cebus_capucinus_imitator`). Dir parsing handles both.

`loss_summ_data.tsv` is headerless, 3 tab-separated columns:
`entry_type {GENE|TRANSCRIPT|PROJECTION}  identifier  status`. Status letters:
`I` intact, `PI` partial-intact, `UL` uncertain-loss, `L` lost, `M` missing,
`PG` paralogous-projection. The walker keeps `GENE` rows by default (`--toga-level`).
**Check the identifier type once** — `zcat <dir>/loss_summ_data.tsv.gz | head` —
because in the hg38 reference these may be Ensembl IDs, not symbols. If so, supply
`--annotation` to translate to symbols (otherwise the symbol query won't match).

## Chromosome + symbol mapping

`loss_summ_data.tsv` has **no chromosome** and may not use gene symbols. Both come
from a one-off reference table passed via `--annotation` (TSV/CSV/BED, gz ok):

```
id<TAB>gene<TAB>chrom         # --annotation-cols id,gene,chrom  (use '-' to skip a col)
```

Build it from the hg38 annotation TOGA used (`TOGAInput/human_hg38` in
github.com/hillerlab/TOGA2) or any hg38 gene table; the `id` column must match the
`loss_summ` identifier, `chrom` gives `chrX`/`chr10`/… `--chrom` in `--toga-dir`
mode requires this (the script errors otherwise).

## Commands

```bash
# genome-wide, chrX only, status-aware, collapsed to species
python intersect_primate_lof.py --query meiosis_repair_genes.txt \
    --toga-dir genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates \
    --annotation hg38_toga_genes.tsv --annotation-cols id,gene,chrom \
    --chrom chrX --out meiosis_chrX

# keep each assembly separate (no species collapse), all chromosomes
python intersect_primate_lof.py --query meiosis_repair_genes.txt \
    --toga-dir .../Primates --toga-label assembly --out meiosis_perasm

# single curated workbook: discover columns first, then run
python intersect_primate_lof.py --query meiosis_repair_genes.txt \
    --table mao2024_suppl/mmc2.xlsx --list-columns
python intersect_primate_lof.py --query meiosis_repair_genes.txt \
    --table mao2024_suppl/mmc2.xlsx --sheet S2 \
    --gene-col Gene --lineage-col Lineage --chrom-col Chr --chrom chrX --out mao_meiosis_chrX
```

## Outputs

- `<out>.toga_long.tsv` — (`--toga-dir`) tidy long table `gene, chrom, lineage,
  status, assembly, raw_id`; the pre-collapse audit trail. Reuse it with `--table`.
- `<out>.intersections.tsv` — one row per (gene, lineage): chrom, most-severe
  status, and `all_statuses` (shows cross-assembly disagreement).
- `<out>.status_matrix.tsv` — gene × lineage → status code.
- `<out>.loss_matrix.tsv` — gene × lineage → 0/1 (1 = loss); **feed to `geneset.py`
  UpSet**.
- (`--table` without `--status-col`: a single `<out>.matrix.tsv` 0/1 presence.)

## Conventions / gotchas

- Symbols normalised before matching: uppercased, trailing `.N` stripped,
  multi-gene cells split on `,;/|`/space.
- Status collapse per (gene, lineage) takes the **most severe** via `STATUS_RANK`
  in the script — so `--toga-label species` reports "lost in ≥1 assembly". Use
  `--toga-label assembly` for per-assembly resolution; `all_statuses` shows the spread.
- `--loss-status` overrides which codes count as loss (default = TOGA loss set;
  `I`/`PI` are not loss).
- `--alias alias.tsv` (`alias<TAB>canonical`) harmonises symbol drift; build from
  HGNC. RefSeq-vs-GENCODE naming differs between the paper workbooks (Yoo 2025 used
  RefSeq GCF_000001405.40-RS_2023_10).
- The walker ignores everything except `loss_summ_data.tsv*`; the big FASTAs in the
  tree are not read.
- Mao 2024 Data S2 is SV-*disruption* presence (no TOGA status) → use `--table`
  without `--status-col`.

## Data provenance / cite

- Mao et al. 2024, *Cell* — Structurally divergent and recurrently mutated regions
  of primate genomes. doi:10.1016/j.cell.2024.01.052 (PMC10947866).
- Yoo et al. 2025, *Nature* 641:401–418 — Complete sequencing of ape genomes.
  doi:10.1038/s41586-025-08816-3 (PMC12058530).
- Kirilenko et al. 2023, *Science* 380:eabn3107 — TOGA. github.com/hillerlab/TOGA2.
  Data: Senckenberg Comparative Genomics, human_hg38_reference.
- Mao et al. 2021, *Nature* (s41586-021-03519-x) — bonobo; origin of the CASK / DMD /
  CACNA1C ILS+LoF seed candidates.

## Research context (for relevance, not required to run)

The chrX-restricted, status-aware run is the one that probes the meiotic-drive /
sex-chromosome-conflict question directly (CASK and DMD are X-linked; MEIG1 is a
spermiogenesis gene flagged among Mao 2024's disrupted set). `loss_matrix.tsv` is
shaped to drop straight into the group's `geneset.py` UpSet visualisation.
