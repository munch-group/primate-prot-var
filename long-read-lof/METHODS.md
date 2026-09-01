# Methods — primate loss-of-function / gene-status intersection analysis

How every input file in this analysis was produced, and how they are combined
into the per-branch intersections and UpSet-ready matrices. Companion to
`CLAUDE.md` (orientation) and `results/README.md` (output manifest). The engine
is `intersect_primate_lof.py`; the genome-wide driver is `run_toga_walk.sh`.

**Driving question.** Which genes in a curated set are lost / disrupted on which
branches of the primate phylogeny, and how do those losses distribute by
chromosome (especially chrX, for the meiotic-drive / sex-chromosome-conflict
angle)?

---

## 0. Background — the analysis in plain terms

*(This section assumes no familiarity with the source papers or tools.)*

**Genes can be lost over evolutionary time.** A gene that was present and working
in a common ancestor can become non-functional in some descendant lineages — it
may be deleted outright, broken by a frameshift or premature stop codon, have its
splice sites destroyed, or be overwritten by a large structural rearrangement. The
gene's remnant may still be visible in the genome (a "pseudogene") or may be gone
entirely. Because every living primate descends from shared ancestors arranged in
a known family tree (the **phylogeny**), we can ask, for any gene, *on which
branches of that tree did it break?* The pattern is informative: a gene lost
repeatedly and independently in several lineages, or lost on a branch where it was
thought essential, points to relaxed selection, a change in biology, or genetic
conflict — rather than random accident.

**How gene loss is detected.** You compare each species' genome against a
well-annotated reference (here, the human genome, build hg38). For every human
gene you ask whether that species has an intact, working copy in the expected
place. The answer is not just yes/no: a copy can be **intact**, **partially
intact**, an **uncertain loss**, clearly **lost**, **missing** simply because that
part of the genome was not assembled (a sequencing gap, not real biology), or
present only as a **paralog** (a copy that landed elsewhere). Distinguishing "truly
lost" from "we couldn't see it" is the central difficulty — assembly gaps in
lower-quality genomes masquerade as losses — and it is why this analysis is careful
about *which* of these categories it counts as a real loss (see §5).

**Three complementary datasets feed the analysis.** No single method sees every
loss, so we combine three independent, recently published primate-genome resources:

- **TOGA** (the backbone). A computational tool that takes the human gene set and,
  using whole-genome alignments plus a trained classifier, projects every human
  gene onto every other primate genome and assigns it one of the status categories
  above. This gives a uniform, automated, *genome-wide* read — **all genes × all
  species** — and is the main source here (169 genome assemblies, 130 species).
- **Mao 2024** (structural-variant view). A study of large structural changes
  (big insertions, deletions, and rearrangements) across primate genomes that
  curated a list of genes **disrupted by fixed structural variants** in particular
  non-human-primate lineages, including genes hit *recurrently* in several
  lineages. It catches breakages caused by large rearrangements that an
  alignment-based method can score differently, so it cross-checks TOGA.
- **Yoo 2025** (the gains, for contrast). Gap-free ("telomere-to-telomere")
  sequencing of ape genomes that identified **lineage-specific new genes** —
  genes *gained* in particular ape lineages. This is the mirror image of loss and
  is included only as a contrast; it is treated as simple presence, not loss.

**Why the focus on chrX, meiosis, and recombination.** The biological motivation
is **genetic conflict over inheritance**. During meiosis (the cell division that
makes eggs and sperm), the two copies of each chromosome are supposed to be passed
on with equal odds, but "selfish" genetic elements can cheat this lottery — a
phenomenon called **meiotic drive**. The X chromosome is a special arena for such
conflict: it is carried differently in the two sexes (males have one X, females
two), so a variant that biases transmission in its own favour can distort the sex
ratio and trigger an evolutionary arms race. Genes that run meiosis and DNA
**recombination**, and especially those on the **X chromosome**, are therefore
expected to evolve fast and to be gained, duplicated, or lost more often than
average. Asking *which meiosis/recombination genes are disrupted on which primate
branches, and whether chrX is enriched* is a direct way to look for the genomic
footprints of this conflict.

**What the pipeline does with all this.** Each dataset is reduced to the same
simple shape — for every (gene, lineage) pair, what is the gene's status? — keyed
on the gene's standardised symbol so the three sources line up. A user supplies a
**query gene list** (e.g. a curated set of recombination genes); the tool reports,
for just those genes, which lineages show a loss/disruption, optionally restricted
to one chromosome (e.g. chrX). Where a species has several genome assemblies, the
most severe call is kept ("lost in at least one assembly"). The outputs are tidy
tables and **0/1 gene × lineage matrices** ready for an UpSet plot (which shows how
losses are shared across lineages). Because the three sources are independent,
agreement between them — the same gene flagged by TOGA's status *and* Mao's
structural-variant list in overlapping lineages — is the strongest signal.

---

## 1. Software environment

- GenomeDK, SLURM. Pure Python + pandas; CPU/RAM-light (runs on a login node or a
  short `srun`; the genome-wide walk takes ~25–30 min single-threaded).
- `python >=3.10`, `pandas`, `openpyxl` (`.xlsx`), `pyarrow` + **`pd-lfs`**
  ([munch-group.org/pd-lfs](https://munch-group.org/pd-lfs)) for the multi-file
  parquet outputs. gz inputs are read transparently (no manual gunzip).
- Environment used here: pixi default env
  (`../.pixi/envs/default/bin/python`, Python 3.12, pandas 2.3, pyarrow + pd-lfs).
- Only the original **downloads** need internet; everything else is offline.

```bash
pixi add "python>=3.10" pandas openpyxl pyarrow pd-lfs
# or: mamba create -n lofx "python>=3.10" pandas openpyxl pyarrow -y && pip install pd-lfs
```

---

## 2. Input data sources and how each file was produced

Four kinds of input feed the analysis: (A) the genome-wide TOGA gene-status tree,
(B) curated paper supplements (Mao 2024, Yoo 2025), (C) a reference annotation
that gives TOGA's gene IDs a symbol and chromosome, and (D) query gene lists
(+ an alias map). Each is described below with the exact processing applied.

### 2A. TOGA gene-status (genome-wide) → `toga_long` + matrices

**Source.** TOGA precomputed primate data (Kirilenko et al. 2023), hg38-referenced,
from the Senckenberg Comparative Genomics mirror. Downloaded with:

```bash
wget -e robots=off -r -np -nH --cut-dirs=4 -A 'loss_summ_data.tsv.gz' \
  https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates/
```

The `-A 'loss_summ_data.tsv.gz'` filter avoids pulling the multi-GB codon/protein
FASTAs. Layout: `…/Primates/<Genus_species>__<Common_name>__<AssemblyID>/`. The
tree here holds **169 assemblies across 130 species** (multiple assemblies per
species are common, e.g. `panTro6`, `HLpanTro7`, `HLpanTroT` for *Pan troglodytes*).

**File read.** `loss_summ_data.tsv(.gz)` is the only file present in every
assembly directory under both of TOGA's file-naming conventions, so the walker
keys on it and ignores everything else. It is **headerless, 3 tab-separated
columns**: `entry_type {GENE|TRANSCRIPT|PROJECTION}  identifier  status`. In the
hg38 reference the identifiers are bare Ensembl gene IDs (`ENSG…`). Status letters:

| code | meaning | counts as loss? |
|---|---|---|
| `I`  | intact | no |
| `PI` | partially intact | no |
| `UL` | uncertain loss | yes (default set) |
| `PG` | paralogous projection | yes (default set) |
| `PM` | partial missing | yes (default set) |
| `L`  | lost | **yes** (strict set) |
| `M`  | missing (assembly gap) | yes (default set) |

**The walk** (`intersect_primate_lof.py --toga-dir`, driven by `run_toga_walk.sh`):
glob `*/loss_summ_data.tsv(.gz)`; parse each directory name into
species / common-name / assembly (robust to `__Primates__` / `__-__` placeholders
and trinomial names); keep `GENE` rows (`--toga-level`, default `GENE`); stamp the
species (`--toga-label species`, the default — collapses a species' assemblies) or
the assembly. This yields the tidy long table **`toga_long`** with columns
`gene, chrom, lineage, status, assembly, raw_id` (~3.3M rows genome-wide).

### 2B. Reference annotation → `data/hg38_toga_genes.tsv`

`loss_summ_data.tsv` carries **no chromosome** and uses `ENSG…` rather than
symbols, so both are supplied by a one-off reference table (`--annotation`,
`--annotation-cols id,gene,chrom`):

- Built from the **Ensembl GRCh38 release-110 GTF** (`gene` features) as
  `id<TAB>gene<TAB>chrom`. Header `id  gene  chrom`; **62,754 genes**.
- `chrom` is stored bare (`1`, `X`); the engine's `norm_chrom` strips any `chr`
  prefix on both sides, so a `--chrom chrX` filter matches `X`.
- Coverage of TOGA's 19,464 genes: **99.9 % chromosome, 99.5 % symbol**
  (843 on chrX). 91 genes have no symbol and remain as `ENSG…` in the outputs.

The `id` column must match the `loss_summ` identifier exactly; the engine uses it
to translate every `ENSG…` to a symbol (`id2sym`) and chromosome (`id2chrom`).

### 2C. Mao et al. 2024 (Cell) — fixed-SV gene disruption

**Source.** `mao2024_suppl.xlsx` (Mao 2024 Data S2; provided directly). Two sheets
are used:

- **`XVI_gene_disrupted`** → `data/mao_gene_disrupted.tsv`. Cleaning:
  - keep only well-formed chromosomes: `CHR` matched against `^chr[0-9XYM]+$`
    (drops malformed footer rows such as `1//1`, `Baboon_SRS660354`);
  - fix the lone case variant `OWL_monkey` → `Owl_monkey`;
  - **explode multi-valued lineage cells** (split on `[;,]+`) to one
    `(gene, lineage)` per row — necessary because the engine's `--table` path
    splits *gene* cells but **not** lineage cells.
  - **canonicalise `GENE_TYPE`** whitespace → underscore (`protein coding` →
    `protein_coding`, `processed pseudogene` → `processed_pseudogene`, …); the
    supplement mixes both forms, which would otherwise split one biotype in two
    on any downstream filter (27 → 23 distinct biotypes; `protein_coding` = 1,334).
  - Result columns: `GENE, GENE_ID, CHR, LINEAGE, CONSEQUENCE, GENE_TYPE`;
    9,375 rows, **4,868 unique genes across 8 NHP lineages** (Bonobo, Chimpanzee,
    Gibbon, Gorilla, Macaque_indian, Marmoset, Orangutan, Owl_monkey); 247 on chrX.
    Coordinates are hg38.
- **`XXXIV_recurrect disrupted genes`** → `data/mao_recurrent_disrupted.tsv`
  (`GENE, LINEAGE` + `GENE_ID, GENE_TYPE` **joined in** from
  `mao_gene_disrupted.tsv` on `GENE` — all 139 are present there; 324 rows,
  **139 genes**) and a flat list `data/mao_recurrent_genes.txt` (candidate curated
  query set).

Mao Data S2 records **SV-disruption presence**, not a TOGA-style status, so it
drives a 0/1 presence `matrix` (no `status_matrix`/`loss_matrix`).

**Carrying gene metadata through.** Mao's list spans all biotypes, so symbol-less
GENCODE loci appear under clone-accession names (`AC007993`, 1,792 of them, ~99.7%
non-coding lncRNA/pseudogene). These are kept (not dropped) and made
joinable/filterable by passing `--meta-cols GENE_ID:gene_id,GENE_TYPE:gene_type`,
which writes `gene_id` (stable Ensembl id, 1:1 with gene) and `gene_type` (biotype)
as leading columns after `chrom` in every Mao output. The 0/1 lineage cells are
identical to the pre-`--meta-cols` matrices (verified). The same flag is applied to
the `results/recombination_genes/mao[_chrX].*` cuts.

### 2D. Yoo et al. 2025 (Nature) — lineage-specific novel genes

**Source.** `data/yoo2025/41586_2025_8816_MOESM4_ESM.xlsx`, **sheet 40** →
`data/yoo_lineage_specific.tsv` (`GENE, LINEAGE, CHR_query, Species,
Amino acid length`; 1,983 rows, **1,075 genes × 9 lineages** = 6 apes + Pan/Pongo
clades + HSA). Two caveats carried into the analysis:

- these are **gene gains, not LoF** — included for contrast, treated as presence;
- `CHR_query` is the **query-assembly contig** (e.g. `CM055450.2`), *not* hg38,
  so **no chrX filter is meaningful** for this source.

### 2E. Query gene lists and the alias map

- `data/all_toga_genes.txt` (103,654) / `data/all_query_genes.txt` (105,999, a
  superset that also includes Mao/Yoo symbols) — used to materialise the
  genome-wide outputs (every gene becomes a "query", so the matrices span the
  whole gene set).
- `recombination_genes.txt` (user-curated, 104 lines) — the curated recombination /
  meiosis query. Resolution against TOGA:
  - **95 genes resolved**;
  - **5 stray non-gene tokens dropped** (phrase fragments left in the list):
    `DSB`, `formation`, `invasion`, `processing`, `Strand`;
  - **13 alias/typo fixes** applied via `--alias` (see below);
  - **3 unresolved, deliberately left out**: `BLMM1A`, `SPO16` (a gene distinct
    from `SHOC1`, which is already in the list), `REDIC1`.
- `data/recombination_alias.tsv` — `alias<TAB>canonical` (no header), built from
  HGNC/standard usage: `CTIP→RBBP8, FANCJ→BRIP1, FIRRM→C1orf112, HEI10→CCNB1IP1,
  HOP2→PSMC3IP, NBS1→NBN, RAD21L→RAD21L1, RAD54→RAD54L, SIX6OS1→C14orf39,
  SWS1→SWSAP1, TO6BL→C11orf80(TOP6BL), XPF→ERCC4, PPCH2→TRIP13`. Each target was
  verified present in TOGA.

**Pre-cut intermediates** (avoid re-walking; produced from `toga_long` by simple
filters): `data/toga_chrX_long.tsv` (chrX rows of `toga_long`, 140,955 rows) and
`data/recomb_toga_long.tsv` (the 95 recombination genes, 15,745 rows), each fed
back through the engine via `--table`.

### 2F. 1000 Genomes long-read SV — human gene loss-of-function

The only **human** loss source (TOGA/Mao use human as the reference; Yoo's `HSA` is
gains). Unlike the others these SVs are **polymorphic** — segregating within humans,
not fixed inter-species differences.

**Source.** 1KG_ONT_VIENNA (Schloissnig et al. 2025; IGSR `data_collections/
1KG_ONT_VIENNA/release/v1.1/giggles-genotyping/`), file
`giggles-genotypes-biallelic.ac0-filtered.vcf.gz` (409 MB, hg38, 967 samples,
sequence-resolved biallelic SVs with genotypes + `AC_Hom`/`AC_Het`). Downloaded to
`data/1kg_ont/` (gitignored). Superpopulation labels from the 2504-sample panel
(`integrated_call_samples_v3.20130502.ALL.panel`) → `data/1kg_ont/superpop_panel.tsv`
(covers 908/967 samples; global AF uses all 967 via INFO, only per-superpop AF drops
the 59 uncovered).

**SV → gene LoF (`annotate_sv_lof.py`).** Pure-python + cyvcf2 + a flat GTF parse
(no bedtools/pyranges). One streaming pass over `tmp/Homo_sapiens.GRCh38.110.gtf.gz`
builds a cached CDS-exon interval model keyed to the same rel-110 gene_ids as
`hg38_toga_genes.tsv` (`data/hg38_cds_intervals.tsv.gz`, `data/hg38_gene_cds_span.tsv`;
281,570 canonical/non-canonical CDS exons over 20,481 coding genes; `canonical` =
MANE_Select/Ensembl_canonical). Each SV is overlapped against CDS via per-chromosome
100 kb binning. Records are sequence-resolved (REF/ALT sequences; type from the ID
token / REF-ALT lengths; reference span = len(REF)); sub-50 bp records skipped.

**Status rules → tiers.** `HC_LOF` = DEL removing a canonical CDS exon / whole-gene
CDS deletion / COMPLEX with ≥50 bp net coding loss (and INV breakpoint-in-CDS for
symbolic inputs); `LC_LOF` = INS-in-CDS / DUP-over-CDS / COMPLEX. Per (gene, SV) the
**strongest** SV (status then AF) is the representative; `--loss-status HC_LOF` is the
default loss set. Human "lineages": `Human_anyLoF`, `Human_commonLoF` (AF≥1%),
`Human_homozygousLoF` (`AC_Hom`>0), and `AFR/AMR/EAS/EUR/SAS_LoF` (allele observed in
that superpop, from genotypes×panel). Output `data/human_1kgp_sv_lof.tsv` (GENE,
GENE_ID, CHR, LINEAGE, STATUS, SVTYPE, CONSEQUENCE, AF, AC, AN, N_HOM_ALT, N_HET,
SV_ID) → engine `--table` → `results/human_1kgp_sv/lof.*`.

```
python annotate_sv_lof.py --build-cds-cache          # one-time, from the GTF
python annotate_sv_lof.py \
  --vcf data/1kg_ont/giggles-genotypes-biallelic.ac0-filtered.vcf.gz \
  --superpop-map data/1kg_ont/superpop_panel.tsv --out data/human_1kgp_sv_lof.tsv
python intersect_primate_lof.py --query data/all_query_genes.txt \
  --table data/human_1kgp_sv_lof.tsv --gene-col GENE --lineage-col LINEAGE \
  --status-col STATUS --chrom-col CHR --loss-status HC_LOF \
  --meta-cols GENE_ID:gene_id,SVTYPE:svtype,CONSEQUENCE:consequence,AF:af,N_HOM_ALT:n_hom_alt \
  --out results/human_1kgp_sv/lof
```

**Validation.** 202,291 SVs → 6,567 hit a CDS → 4,544 genes; **1,733 genes carry a
homozygous HC-LoF**, matching gnomAD's ~1,800 genes with homozygous pLoF. Essential
genes (RAD51, TP53, BRCA2) show no homozygous LoF. **Caveat:** this biallelic
graph-genotyped callset under-ascertains segmental-duplication–mediated common gene
deletions (GSTM1/GSTT1/UGT2B17/CYP2D6/AMY1A absent or rare-only — confirmed by direct
VCF inspection); GSTT1/LILRA3 are also outside the TOGA gene universe.

---

## 3. How the files are combined — the intersection engine

`intersect_primate_lof.py` reduces any source to a common shape: a mapping
`store[gene][lineage] → {set of status codes}` plus `chrom_of[gene]`. The shared
join key across all sources is the **normalised gene symbol**.

**Symbol / status / chromosome normalisation.**
- Symbols: uppercased, trailing version (`.N`) stripped; multi-gene cells split on
  `[,;/|\s]`. Then mapped through the **alias map** (`canon()`), applied to *both*
  the query and the data so drift between RefSeq/GENCODE/TOGA naming is harmonised.
- Statuses: uppercased, spaces/hyphens → `_` (so `partially intact` ≡ `PI`).
- Chromosomes: lowercased, leading `chr` stripped, re-uppercased (`chrX` ≡ `X`).

**TOGA path (`--toga-dir`).** For every walked row: `raw_id` → symbol via the
annotation `id2sym` → alias `canon()`; chromosome via `id2chrom`; optional
`--chrom` filter; accumulate the status into `store[gene][lineage]`. All rows are
also emitted verbatim as `toga_long` (the reusable audit trail).

**Table path (`--table`).** Read a workbook / TSV / **or a pd-lfs parquet
directory** (`read_parquet`, so a previously written `toga_long.parquet/` can be
re-cut without re-walking). For each row: optional `--chrom` filter; explode the
gene cell; take the lineage from `--lineage-col` (not split); status from
`--status-col` (or the literal `HIT` when no status column exists, i.e. presence
mode); accumulate.

**Collapse to one call per (gene, lineage).** Multiple assemblies of a species (or
multiple SV rows) give a set of statuses; the engine keeps the **most severe** via
`STATUS_RANK`:

```
I (0) < PI (1) < UL ≈ PG (2) < PM (3) < L (4) < M (5)
```

So `--toga-label species` reports "lost in ≥1 assembly"; the full per-assembly
spread is preserved in `toga_long.assembly` and in the `all_statuses` column of
`intersections`. Use `--toga-label assembly` for per-assembly resolution.

**Loss call.** A (gene, lineage) is a "loss" iff its most-severe status is in the
loss set. Default set = `{UL, PG, PM, L, M}` (and spelled-out synonyms);
`I`/`PI` are never loss. Override with `--loss-status` (e.g. `--loss-status L` for
true loss only — see the caveat in §5).

**Outputs** (computed identically from `store`, restricted to query ∩ data):
- `intersections` — one row per (gene, lineage): chrom, most-severe status,
  `all_statuses`;
- `status_matrix` — gene × lineage → most-severe status code (when status-aware);
- `loss_matrix` — gene × lineage → 0/1 (1 = loss); UpSet-ready, feeds `geneset.py`;
- `matrix` — gene × lineage → 0/1 presence (the `--table`, no-status fallback,
  e.g. Mao/Yoo);
- `toga_long` — only in `--toga-dir` mode.

**Output format.** `--toga-dir` writes each output as a **pd-lfs multi-file
parquet directory** `<out>.<name>.parquet/` (part-`*`.parquet shards ≤50 MB +
`_manifest.json`) so the genome-wide tables stay under GitHub's file-size limit;
read them back with `pd_lfs.parquet.read_parquet` (local path or https URL, dtypes
restored from the manifest). `--table` writes plain `<name>.tsv`. (Implementation
note: string columns are cast to pandas `"string"` dtype before writing so
high-cardinality columns such as `gene`/`raw_id` get a concrete Arrow schema
rather than being mis-inferred as `null` from an empty slice.)

**Output location.** Outputs for one curated query gene list are grouped in a
**per-query folder** `results/<query-stem>/`: a **bare** `--out` label (no `/`) is
written to `results/<query-stem>/<label>.*`, so every cut of `recombination_genes.txt`
(`toga`, `toga_strictL`, `toga_chrX`, `mao`, …) collects in `results/recombination_genes/`.
An `--out` value **containing `/`** (or an absolute path) is written verbatim, flat —
used for the genome-wide walk (`results/toga_primates_species.*`, kept at the `results/`
root) and for the full-source reformats grouped by source (`results/mao2024/`,
`results/yoo2025/`). `--results-dir` changes the base directory (default `results`).

**Cross-source combination (`merge_loss_matrices.py`).** The sources use different,
non-commensurable lineage vocabularies (TOGA binomial species, Mao common names,
`Human_*`/superpop tiers), so they are combined as a **per-gene** union, never a
cell-by-cell average. The merge outer-joins each source's `loss_matrix`/`matrix` on
the `gene` key and **prefixes** every lineage column with its source (`toga:`,
`mao:`, `human:`) — provenance stays explicit and columns never collide. Genes
absent from a source are filled 0; a companion **coverage mask**
(`<out>.coverage.tsv`, `<src>:_tested`) distinguishes "tested & not lost" from "gene
not in that source" (for the human source the tested universe is *all* coding genes,
via `--universe human=data/hg38_gene_cds_span.tsv`). Yoo is excluded by default
(gains, non-hg38 coordinates; add with `--include-yoo`). Output:
`results/merged/loss_matrix.tsv` (genome-wide: 23,095 genes × 130 `toga:` + 8 `mao:`
+ 8 `human:` columns) and `results/recombination_genes/merged.*`. Agreement across
sources is the corroboration signal — e.g. on the recombination set IHO1/MEI4/SYCP1
are homozygous-LoF in humans **and** lost in a TOGA primate.

```bash
python merge_loss_matrices.py \
  --matrix toga=results/toga_primates_species.loss_matrix.parquet \
  --matrix mao=results/mao2024/gene_disrupted.matrix.tsv \
  --matrix human=results/human_1kgp_sv/lof.loss_matrix.tsv \
  --universe human=data/hg38_gene_cds_span.tsv \
  --out results/merged/loss_matrix.tsv
```

---

## 4. Reproducing

```bash
# genome-wide, species-collapsed, parquet outputs (the toga_primates_species.* set)
./run_toga_walk.sh data/all_toga_genes.txt results/toga_primates_species

# a curated subset: a bare --out label nests under results/<query-stem>/
#   -> results/recombination_genes/toga_chrX_strictL.{intersections,loss_matrix,...}
./run_toga_walk.sh recombination_genes.txt toga_chrX_strictL \
    --alias data/recombination_alias.tsv --chrom chrX --loss-status L

# re-cut an existing genome-wide table WITHOUT re-walking (reads the parquet dir);
# bare --out 'toga_strictL' -> results/recombination_genes/toga_strictL.*
../.pixi/envs/default/bin/python intersect_primate_lof.py \
    --query recombination_genes.txt --alias data/recombination_alias.tsv \
    --table results/toga_primates_species.toga_long.parquet \
    --gene-col gene --lineage-col lineage --status-col status --chrom-col chrom \
    --loss-status L --out toga_strictL

# a curated workbook (discover columns first, then run)
../.pixi/envs/default/bin/python intersect_primate_lof.py \
    --query recombination_genes.txt --table mao2024_suppl.xlsx --list-columns
```

See `run_toga_walk.sh --help` (run with no args) for env overrides
(`PYTHON`, `TOGA_DIR`, `ANNOTATION`, `ANNO_COLS`, `LABEL`) and the SLURM hint.

---

## 5. Caveats

- **Loss definition dominates the result for essential-gene sets.** For the
  recombination set the TOGA status mix is ~79.7 % `I`, 10.6 % `PI`, 5.0 % `UL`,
  0.9 % `L`, 0.5 % `M`. The **default** loss set (`UL+PG+PM+L+M`) flags *all 95*
  genes as "lost somewhere", dominated by `M` (assembly gaps) and `UL` in
  low-quality assemblies (Rhinopithecus/Nasalis top the total-loss ranking). Use
  `--loss-status L` for real LoF (51/95 genes); even then, true-`L` calls for
  essential genes (RAD51, BRCA1…) deserve manual scrutiny — they can be
  paralog/assembly artefacts. `status_matrix` keeps the raw codes for any threshold.
- **Species collapse** = most-severe across that species' assemblies; per-assembly
  detail is in `toga_long.assembly` / `--toga-label assembly`.
- **Symbol drift.** Mao/Yoo use GENCODE/RefSeq symbols; supply `--alias` to
  harmonise against TOGA's hg38 symbols where needed (done for the recombination set).
- **Yoo coordinates are query-assembly contigs**, not hg38 — no chrX filtering.
- **Mao is SV-disruption presence**, not status — interpret its `matrix` as
  "disrupted in this lineage", independent of TOGA's loss calls.

---

## 6. Provenance / cite

- **TOGA** — Kirilenko et al. 2023, *Science* 380:eabn3107. Data: Senckenberg
  Comparative Genomics, human_hg38_reference. github.com/hillerlab/TOGA2.
- **Mao 2024** — Mao et al. 2024, *Cell* — Structurally divergent and recurrently
  mutated regions of primate genomes. doi:10.1016/j.cell.2024.01.052 (PMC10947866).
- **Yoo 2025** — Yoo et al. 2025, *Nature* 641:401–418 — Complete sequencing of ape
  genomes. doi:10.1038/s41586-025-08816-3 (PMC12058530).
- **Mao 2021** — Mao et al. 2021, *Nature* (s41586-021-03519-x) — bonobo; origin of
  the CASK / DMD / CACNA1C ILS+LoF seed candidates.
- **Annotation** — Ensembl GRCh38 release-110 GTF (gene features).
- **pd-lfs** — Munch group, munch-group.org/pd-lfs (chunked parquet for GitHub).
