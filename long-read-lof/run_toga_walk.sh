#!/usr/bin/env bash
#
# run_toga_walk.sh
# ================
# Genome-wide TOGA walk for a query gene set, writing pd-lfs multi-file parquet
# outputs (git-storable, <50 MB/shard) via intersect_primate_lof.py --toga-dir.
#
# Produces, under the OUT prefix (see "OUT controls..." below):
#   <OUT>.toga_long.parquet/      full tidy table (gene,chrom,lineage,status,assembly,raw_id)
#   <OUT>.intersections.parquet/  one row per (query_gene, lineage)
#   <OUT>.status_matrix.parquet/  gene x lineage -> most-severe status
#   <OUT>.loss_matrix.parquet/    gene x lineage -> 0/1 (UpSet-ready)
# (read any of them back with pd_lfs.parquet.read_parquet, locally or over https.)
#
# Usage:
#   ./run_toga_walk.sh QUERY_GENES [OUT] [extra intersect_primate_lof.py args...]
#
# OUT controls where files land (see intersect_primate_lof.py --out):
#   * a bare label, or omitted (defaults to 'toga')  -> nested in a per-query folder
#       results/<QUERY-stem>/<label>.*   (a curated set like recombination_genes.txt
#       collects in results/recombination_genes/)
#   * a value containing '/'                          -> used verbatim, flat
#       (the genome-wide convention, e.g. results/toga_primates_species)
#
# Examples:
#   # full genome-wide materialisation (all genes), species-collapsed, flat at results/ root
#   ./run_toga_walk.sh data/all_toga_genes.txt results/toga_primates_species
#
#   # a curated gene set, chrX only, strict true-loss, with an alias map
#   #   -> results/recombination_genes/toga.{intersections,loss_matrix,status_matrix}.parquet/
#   ./run_toga_walk.sh recombination_genes.txt toga \
#       --alias data/recombination_alias.tsv --chrom chrX --loss-status L
#
#   # OUT omitted -> default label 'toga', still nested by query stem
#   ./run_toga_walk.sh recombination_genes.txt
#
#   # keep each assembly separate instead of collapsing per species
#   ./run_toga_walk.sh data/all_toga_genes.txt results/toga_primates_perasm \
#       --toga-label assembly
#
# Overridable via environment variables (defaults shown):
#   PYTHON      python interpreter        (../.pixi/envs/default/bin/python, else `python`)
#   TOGA_DIR    downloaded …/Primates/    (genome.senckenberg.de/.../Primates)
#   ANNOTATION  id,gene,chrom table       (data/hg38_toga_genes.tsv)
#   ANNO_COLS   --annotation-cols value   (id,gene,chrom)
#   LABEL       --toga-label value        (species)
#
# Tip: the walk reads 169 gzipped loss_summ files and takes ~25-30 min. To submit
# it to SLURM:  sbatch --mem=8g -c 1 -t 02:00:00 --wrap "./run_toga_walk.sh ... "
#
set -euo pipefail

# Resolve to this script's directory so relative paths work from anywhere.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ $# -lt 1 ]]; then
    awk 'NR>1 && /^#/{print} NR>1 && !/^#/{exit}' "${BASH_SOURCE[0]}"   # header comment as usage
    exit 2
fi

QUERY="$1"; shift
# OUT is optional: take the next arg unless it's an option (starts with '-') or absent.
# A bare label nests under results/<query-stem>/ (intersect_primate_lof.py --out).
OUT="toga"
if [[ $# -ge 1 && "$1" != -* ]]; then OUT="$1"; shift; fi
QSTEM="$(basename "${QUERY%.*}")"   # 'data/all_toga_genes.txt' -> 'all_toga_genes'
# remaining "$@" pass straight through to intersect_primate_lof.py

# Defaults (override by exporting the variable before calling).
PYTHON="${PYTHON:-$SCRIPT_DIR/../.pixi/envs/default/bin/python}"
[[ -x "$PYTHON" ]] || PYTHON="python"
TOGA_DIR="${TOGA_DIR:-genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates}"
ANNOTATION="${ANNOTATION:-data/hg38_toga_genes.tsv}"
ANNO_COLS="${ANNO_COLS:-id,gene,chrom}"
LABEL="${LABEL:-species}"

# Fail early with a clear message if an input is missing.
for f in "$QUERY" "$ANNOTATION" "$TOGA_DIR" intersect_primate_lof.py; do
    [[ -e "$f" ]] || { echo "ERROR: missing required path: $f" >&2; exit 1; }
done
# (intersect_primate_lof.py creates the output folder itself.)

# Best-effort display of the output prefix (assumes the default --results-dir).
if [[ "$OUT" == */* ]]; then DEST="$OUT"; else DEST="results/$QSTEM/$OUT"; fi

echo "=== TOGA genome-wide walk ==========================================="
echo "  python      : $PYTHON"
echo "  query genes : $QUERY"
echo "  toga tree   : $TOGA_DIR"
echo "  annotation  : $ANNOTATION  (cols: $ANNO_COLS)"
echo "  toga-label  : $LABEL"
echo "  out / label : $OUT"
echo "  -> writes   : ${DEST}.*"
echo "  extra args  : $*"
echo "  started     : $(date '+%Y-%m-%d %H:%M:%S')"
echo "====================================================================="

SECONDS=0
"$PYTHON" intersect_primate_lof.py \
    --query "$QUERY" \
    --toga-dir "$TOGA_DIR" \
    --annotation "$ANNOTATION" --annotation-cols "$ANNO_COLS" \
    --toga-label "$LABEL" \
    --out "$OUT" \
    "$@"

echo "====================================================================="
echo "  done in $((SECONDS/60))m $((SECONDS%60))s  ->  ${DEST}.*.parquet/  (see 'wrote:' above)"
