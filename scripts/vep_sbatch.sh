#!/bin/bash
#SBATCH --job-name=vep_annotation
#SBATCH --array=1-24                  # 1-22 + X(23) + Y(24)
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/vep_%a_%j.out
#SBATCH --error=logs/vep_%a_%j.err
#SBATCH --account=xy-drive

# ── map array index to chromosome name ────────────────────────────────────────
declare -a CHROMS=($(printf 'chr%s ' {1..22}) chrX chrY)
CHROM=${CHROMS[$((SLURM_ARRAY_TASK_ID - 1))]}

# ── paths — edit these ────────────────────────────────────────────────────────
VCF_DIR="vep_data/vcf"
OUT_DIR="vep_data/vcf_annotated"
CACHE_DIR="vep_data/.vep"

IN_VCF="${VCF_DIR}/homo_sapiens-${CHROM}.vcf.gz"
OUT_VCF="${OUT_DIR}/homo_sapiens-${CHROM}.vep.vcf.gz"

mkdir -p "$OUT_DIR" logs

# ── skip if output already exists ─────────────────────────────────────────────
if [[ -f "$OUT_VCF" ]]; then
    echo "[$(date)] $CHROM already annotated, skipping."
    exit 0
fi

echo "[$(date)] Starting VEP for $CHROM"

# ── filter out unsupported <.> variants before VEP ────────────────────────────
FILTERED_VCF="/tmp/${CHROM}_filtered_${SLURM_JOB_ID}.vcf.gz"

pixi run bcftools view \
    -e 'ALT="<.>"' \
    "$IN_VCF" \
    -O z -o "$FILTERED_VCF"

pixi run tabix -p vcf "$FILTERED_VCF"

# ── run VEP ───────────────────────────────────────────────────────────────────
pixi run vep \
    --input_file      "$FILTERED_VCF" \
    --output_file     "$OUT_VCF" \
    --vcf \
    --compress_output bgzip \
    --cache \
    --offline \
    --dir_cache       "$CACHE_DIR" \
    --assembly        GRCh38 \
    --fork            "$SLURM_CPUS_PER_TASK" \
    --buffer_size     100000 \
    --no_intergenic \
    --no_stats \
    --symbol \
    --canonical \
    --biotype \
    --hgvs \
    --sift b \  
    --polyphen b \
    --af_gnomade \
    --variant_class \
    --fasta "vep_data/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"


# ── index output ──────────────────────────────────────────────────────────────
pixi run tabix -p vcf "$OUT_VCF"

# ── clean up temp file ────────────────────────────────────────────────────────
rm -f "$FILTERED_VCF" "${FILTERED_VCF}.tbi"

echo "[$(date)] Finished VEP for $CHROM → $OUT_VCF"