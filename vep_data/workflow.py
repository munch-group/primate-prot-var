from pathlib import Path
from gwf import Workflow, AnonymousTarget
from gwf.executors import Pixi
import os
import re

def modify_path(path, **kwargs):
    """
    Utility function for modifying file paths substituting
    the directory (dir), base name (base), or file suffix (suffix).
    """
    was_path = False
    if isinstance(path, Path):
        was_path = True
        path = str(path)

    for key in ['dir', 'base', 'suffix']:
        kwargs.setdefault(key, None)
    assert len(kwargs) == 3

    par, name = os.path.split(path)
    name_no_suffix, suf = os.path.splitext(name)
    if type(kwargs['suffix']) is str:
        suf = kwargs['suffix']
    if kwargs['dir'] is not None:
        par = kwargs['dir']
    if kwargs['base'] is not None:
        name_no_suffix = kwargs['base']

    new_path = os.path.join(par, name_no_suffix + suf)
    if type(kwargs['suffix']) is tuple:
        assert len(kwargs['suffix']) == 2
        new_path, nsubs = re.subn(r'{}$'.format(kwargs['suffix'][0]), kwargs['suffix'][1], new_path)
        assert nsubs == 1, nsubs

    if was_path:
        new_path = Path(new_path)

    return new_path


#   rsync -avP --include='*.vcf.gz' --include='*.vcf.gz.tbi' rsync://ftp.ensembl.org/ensembl/pub/current_variation/vcf/homo_sapiens/  ./vcf/

#   rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz  vep_cache/
#   tar -xzf vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz -C vep_cache/


# ── configuration ─────────────────────────────────────────────────────────────

ACCOUNT       = "xy-brain"           # ← change to your SLURM account
PIXI_ROOT     = Path(".")                # pixi project root (where pixi.toml lives)
ENSEMBL       = 115                      # Ensembl release
CACHE_DIR     = Path("steps/vep_cache")
VCF_DIR       = Path("steps/vcf")
FILT_DIR      = Path("steps/vcf_filtered")
CHUNK_DIR     = Path("steps/vcf_chunks")
ANN_CHUNK_DIR = Path("steps/vcf_annotated_chunks")
ANN_DIR       = Path("steps/vcf_annotated")
PARQUET       = Path("steps/parquet")
FASTA_DIR     = CACHE_DIR
FASTA_GZ      = FASTA_DIR / "_Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
FASTA_BGZ     = FASTA_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
FASTA_RAW     = FASTA_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
FASTA_INDEXED = FASTA_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.fai"
FASTA_INDEX   = FASTA_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.gzi"

CHROMS     = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# Window size (in bp) for splitting VEP runs per chromosome
VEP_WINDOW_SIZE = 10_000_000  # 10 Mb

# Approximate chromosome lengths (GRCh38) for defining windows
CHROM_LENGTHS = {
    "chr1": 249_000_000, "chr2": 243_000_000, "chr3": 198_000_000,
    "chr4": 190_000_000, "chr5": 182_000_000, "chr6": 171_000_000,
    "chr7": 160_000_000, "chr8": 146_000_000, "chr9": 139_000_000,
    "chr10": 134_000_000, "chr11": 135_000_000, "chr12": 134_000_000,
    "chr13": 115_000_000, "chr14": 107_000_000, "chr15": 102_000_000,
    "chr16": 90_000_000,  "chr17": 84_000_000,  "chr18": 81_000_000,
    "chr19": 59_000_000,  "chr20": 65_000_000,  "chr21": 47_000_000,
    "chr22": 51_000_000,  "chrX": 157_000_000,  "chrY": 58_000_000,
}

def chrom_windows(chrom: str, window_size: int = VEP_WINDOW_SIZE):
    """Return a list of (start, end) 1-based windows covering a chromosome."""
    length = CHROM_LENGTHS[chrom]
    windows = []
    for start in range(1, length + 1, window_size):
        end = min(start + window_size - 1, length)
        windows.append((start, end))
    return windows

# ── workflow ───────────────────────────────────────────────────────────────────

gwf = Workflow(
    defaults={
        "account"  : ACCOUNT,
        "walltime" : "02:00:00",
        "memory"   : "8gb",
        "cores"    : 1,
    },
    executor=Pixi(project=str(PIXI_ROOT)),
)

# ── helpers ────────────────────────────────────────────────────────────────────

def ensure_dirs(*dirs):
    return "\n".join(f"mkdir -p {d}" for d in dirs)


# ══════════════════════════════════════════════════════════════════════════════
# ONE-OFF TARGETS: FASTA reference
# ══════════════════════════════════════════════════════════════════════════════

def download_fasta(fasta_gz: Path) -> AnonymousTarget:
    """Download GRCh38 primary assembly FASTA from Ensembl FTP."""
    url = (
        f"https://ftp.ensembl.org/pub/release-{ENSEMBL}"
        f"/fasta/homo_sapiens/dna"
        f"/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    )
    return AnonymousTarget(
        inputs  = [],
        outputs = [str(fasta_gz)],
        options = {"memory": "4gb", "walltime": "18:00:00", "cores": 1},
        spec    = f"""
{ensure_dirs(fasta_gz.parent)}
curl -L -o {fasta_gz} '{url}'
""",
    )


def prepare_fasta(fasta_gz: Path, fasta_fai: Path) -> AnonymousTarget:
    """
    gunzip the downloaded FASTA, recompress with bgzip, and index with samtools.
    VEP requires bgzip (not gzip) for HGVS lookups.
    """
    fasta_raw = modify_path(fasta_gz, suffix=('.gz', ''))
    fasta_bgz = modify_path(fasta_fai, suffix=('.fai', ''))
    return AnonymousTarget(
        inputs  = [str(fasta_gz)],
        outputs = [str(fasta_fai), 
                   str(modify_path(fasta_fai, suffix=('.fai', '.gzi')))
                   ],
        options = {"memory": "42gb", "walltime": "18:00:00", "cores": 1},
        spec    = f"""
pixi run gzip -d -c -f -k {fasta_gz} > {fasta_raw}
pixi run bgzip -f -@ 4 {fasta_raw}
pixi run samtools faidx {fasta_bgz}
""",
    )


gwf.target_from_template(
    name     = "DownloadFasta",
    template = download_fasta(FASTA_GZ),
)

# # bgzipped path: same stem, still .fa.gz but now bgzip-compressed
# FASTA_BGZ = FASTA_RAW.with_suffix(".fa.gz")


fasta_target = gwf.target_from_template(
    name     = "PrepareFasta",
    template = prepare_fasta(FASTA_GZ, FASTA_INDEXED),
)


# ══════════════════════════════════════════════════════════════════════════════
# PER-CHROMOSOME TEMPLATES
# ══════════════════════════════════════════════════════════════════════════════

def filter_vcf(chrom: str, vcf_dir: Path, filt_dir: Path) -> AnonymousTarget:
    """
    Remove variants with ALT='<.>' (unsupported structural variants that
    cause VEP warnings) using bcftools.
    """
    in_vcf   = vcf_dir  / f"homo_sapiens-{chrom}.vcf.gz"
    out_vcf  = filt_dir / f"homo_sapiens-{chrom}.filtered.vcf.gz"
    out_tbi  = Path(str(out_vcf) + ".tbi")
    return AnonymousTarget(
        inputs  = [str(in_vcf)],
        outputs = [str(out_vcf), str(out_tbi)],
        options = {"memory": "8gb", "walltime": "18:00:00", "cores": 1},
        spec    = f"""
{ensure_dirs(filt_dir)}
pixi run bcftools view -e 'ALT="<.>"' {in_vcf} -O z -o {out_vcf}
pixi run tabix -p vcf {out_vcf}
""",
    )


def split_vcf_chunk(
    chrom    : str,
    start    : int,
    end      : int,
    filt_dir : Path,
    chunk_dir: Path,
) -> AnonymousTarget:
    """
    Extract a genomic window from a filtered VCF using bcftools.

    Multiallelic records are split into biallelic ones (norm -m -any):
    multiallelic indels crash VEP 115 (substr error in display_codon,
    e.g. rs758603179), and the forked run then hangs until walltime.
    Records with POS before the window are dropped: bcftools -r also
    returns indels that start before the window and overlap it, which
    would duplicate them in two chunks and break sorting on concat.
    """
    in_vcf  = filt_dir  / f"homo_sapiens-{chrom}.filtered.vcf.gz"
    out_vcf = chunk_dir / f"homo_sapiens-{chrom}_{start}_{end}.vcf.gz"
    out_tbi = Path(str(out_vcf) + ".tbi")
    # bcftools uses chrom names without "chr" prefix in Ensembl VCFs
    bcf_chrom = chrom.replace("chr", "")
    return AnonymousTarget(
        inputs  = [str(in_vcf)],
        outputs = [str(out_vcf), str(out_tbi)],
        options = {"memory": "4gb", "walltime": "02:00:00", "cores": 1},
        spec    = f"""
{ensure_dirs(chunk_dir)}
set -o pipefail
pixi run bcftools view -r {bcf_chrom}:{start}-{end} {in_vcf} -O u \\
    | pixi run bcftools norm -m -any -O u - \\
    | pixi run bcftools view -e 'POS<{start}' -O z -o {out_vcf} -
pixi run tabix -p vcf {out_vcf}
""",
    )


def run_vep_chunk(
    chrom     : str,
    start     : int,
    end       : int,
    chunk_dir : Path,
    ann_dir   : Path,
    fasta_bgz : Path,
    cache_dir : Path,
    cores     : int = 4,
) -> AnonymousTarget:
    """
    Annotate a single chunk VCF with Ensembl VEP using the local cache.
    """
    in_vcf  = chunk_dir / f"homo_sapiens-{chrom}_{start}_{end}.vcf.gz"
    out_vcf = ann_dir   / f"homo_sapiens-{chrom}_{start}_{end}.vep.vcf.gz"
    out_tbi = Path(str(out_vcf) + ".tbi")
    return AnonymousTarget(
        inputs  = [str(in_vcf), str(fasta_bgz)],
        outputs = [str(out_vcf), str(out_tbi)],
        options = {
            "memory"  : "16gb",
            "walltime": "4-00:00:00",
            "cores"   : cores,
        },
        spec = f"""
{ensure_dirs(ann_dir)}

# Count variants (non-header lines). If zero, copy the header-only
# VCF through so downstream concat still works.
N_VARIANTS=$(pixi run bcftools view -H {in_vcf} | wc -l)

if [ "$N_VARIANTS" -eq 0 ]; then
    echo "No variants in chunk — writing header-only VCF"
    pixi run bcftools view -h {in_vcf} | pixi run bgzip -c > {out_vcf}
    pixi run tabix -p vcf {out_vcf}
else
    pixi run vep \\
        --input_file      {in_vcf} \\
        --output_file     {out_vcf} \\
        --force_overwrite \\
        --verbose \\
        --vcf \\
        --compress_output bgzip \\
        --cache \\
        --offline \\
        --dir_cache       {cache_dir} \\
        --assembly        GRCh38 \\
        --fasta           {fasta_bgz} \\
        --fork            {cores} \\
        --buffer_size     5000 \\
        --no_intergenic \\
        --no_stats \\
        --symbol \\
        --canonical \\
        --biotype \\
        --hgvs \\
        --sift            b \\
        --polyphen        b \\
        --af_gnomade \\
        --variant_class

    pixi run tabix -p vcf {out_vcf}
fi
""",
    )


def concat_vep(
    chrom   : str,
    windows : list,
    ann_chunk_dir : Path,
    ann_dir : Path,
) -> AnonymousTarget:
    """
    Concatenate per-chunk VEP outputs back into a single per-chromosome VCF.

    -a -D: chunks produced before the POS filter in SplitVCF contain
    boundary-overlapping indels present in two adjacent chunks; -a keeps
    the output sorted across those overlaps and -D drops the duplicates
    (plain concat leaves the file unsorted there and tabix fails).
    """
    in_vcfs = [
        ann_chunk_dir / f"homo_sapiens-{chrom}_{s}_{e}.vep.vcf.gz"
        for s, e in windows
    ]
    out_vcf = ann_dir / f"homo_sapiens-{chrom}.vep.vcf.gz"
    out_tbi = Path(str(out_vcf) + ".tbi")
    input_list = " ".join(str(v) for v in in_vcfs)
    return AnonymousTarget(
        inputs  = [str(v) for v in in_vcfs],
        outputs = [str(out_vcf), str(out_tbi)],
        options = {"memory": "72gb", "walltime": "23:00:00", "cores": 1},
        spec    = f"""
{ensure_dirs(ann_dir)}
pixi run bcftools concat -a -D {input_list} -O z -o {out_vcf}
pixi run tabix -p vcf {out_vcf}
""",
    )


def vcf_to_parquet(
    chrom      : str,
    ann_dir    : Path,
    parquet_dir: Path,
) -> AnonymousTarget:
    """
    Parse CSQ annotations from the VEP-annotated VCF and write a
    Hive-partitioned Parquet file: parquet/chrom=chrN/part-0.parquet
    The sentinel file (.done) signals completion to gwf.
    """
    in_vcf  = ann_dir    / f"homo_sapiens-{chrom}.vep.vcf.gz"
    out_dir = parquet_dir / f"chrom={chrom}"
    sentinel = out_dir   / ".done"

    # Inline Python script — avoids a separate script file
    script = f"""
import re, sys
import cyvcf2
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path

IN_VCF      = "{in_vcf}"
OUT_DIR     = Path("{out_dir}")
OUT_DIR.mkdir(parents=True, exist_ok=True)
SENTINEL    = Path("{sentinel}")

def get_csq_fields(vcf):
    for line in str(vcf.raw_header).split('\\n'):
        if 'ID=CSQ' in line:
            m = re.search(r'Format: ([^"]+)', line)
            if m:
                return m.group(1).strip().split('|')
    return []

vcf        = cyvcf2.VCF(IN_VCF)
csq_fields = get_csq_fields(vcf)
if not csq_fields:
    print(f"ERROR: No CSQ field found in {{IN_VCF}}", file=sys.stderr)
    sys.exit(1)

rows = []
for v in vcf:
    csq_raw = v.INFO.get('CSQ', '')
    if not csq_raw:
        continue
    for transcript in csq_raw.split(','):
        record = dict(zip(csq_fields, transcript.split('|')))
        if not record.get('SYMBOL'):
            continue
            rows.append({{
                'chrom'            : v.CHROM,
                'pos'              : v.POS,
                'ref'              : v.REF,
                'alt'              : ','.join(v.ALT),
                'variant_id'       : v.ID or '.',
                'gene_symbol'      : record.get('SYMBOL', ''),
                'gene_id'          : record.get('Gene', ''),
                'transcript_id'    : record.get('Feature', ''),
                'consequence_terms': record.get('Consequence', ''),
                'impact'           : record.get('IMPACT', ''),
                'biotype'          : record.get('BIOTYPE', ''),
                'canonical'        : 1 if record.get('CANONICAL') == 'YES' else 0,
                'hgvsc'            : record.get('HGVSc', ''),
                'hgvsp'            : record.get('HGVSp', ''),
                'sift'             : record.get('SIFT', ''),
                'polyphen'         : record.get('PolyPhen', ''),
                'af_gnomad'        : record.get('gnomADe_AF', '') or record.get('AF', ''),
                'clin_sig'         : record.get('CLIN_SIG', ''),
            }})

if rows:
    table = pa.Table.from_pylist(rows)
    # Sort by pos so DuckDB row-group stats enable fast range queries
    import pyarrow.compute as pc
    table = table.sort_by([('pos', 'ascending')])
    pq.write_table(
        table,
        OUT_DIR / 'part-0.parquet',
        compression='zstd',
        use_dictionary=True,
    )
    print(f"Wrote {{len(rows)}} rows to {{OUT_DIR}}")
else:
    print(f"No genic rows found in {{IN_VCF}}")

SENTINEL.touch()
"""

    return AnonymousTarget(
        inputs  = [str(in_vcf)],
        outputs = [str(sentinel)],
        options = {"memory": "16gb", "walltime": "02:00:00", "cores": 1},
        spec    = f"""
{ensure_dirs(out_dir)}
pixi run python - << 'PYEOF'
{script}
PYEOF
""",
    )


# ══════════════════════════════════════════════════════════════════════════════
# WIRE UP PER-CHROMOSOME TARGETS
# ══════════════════════════════════════════════════════════════════════════════

all_parquet_sentinels = []

for chrom in CHROMS:

    # Step 1 — filter
    filter_target = gwf.target_from_template(
        name     = f"FilterVCF_{chrom}",
        template = filter_vcf(chrom, VCF_DIR, FILT_DIR),
    )

    # Step 2 — split, run VEP per chunk, concat
    windows = chrom_windows(chrom)
    for start, end in windows:
        gwf.target_from_template(
            name     = f"SplitVCF_{chrom}_{start}",
            template = split_vcf_chunk(chrom, start, end, FILT_DIR, CHUNK_DIR),
        )
        gwf.target_from_template(
            name     = f"RunVEP_{chrom}_{start}",
            template = run_vep_chunk(chrom, start, end, CHUNK_DIR, ANN_CHUNK_DIR,
                                     FASTA_BGZ, CACHE_DIR, cores=4),
        )
    gwf.target_from_template(
        name     = f"ConcatVEP_{chrom}",
        template = concat_vep(chrom, windows, ANN_CHUNK_DIR, ANN_DIR),
    )

    # Step 3 — Parquet
    parquet_target = gwf.target_from_template(
        name     = f"VCFToParquet_{chrom}",
        template = vcf_to_parquet(chrom, ANN_DIR, PARQUET),
    )

    all_parquet_sentinels.append(str(PARQUET / f"chrom={chrom}" / ".done"))


# ══════════════════════════════════════════════════════════════════════════════
# FINAL ENDPOINT — waits for all chromosomes, prints summary
# ══════════════════════════════════════════════════════════════════════════════

gwf.target(
    name    = "AllDone",
    inputs  = all_parquet_sentinels,
    outputs = [str(PARQUET / "README.md")],
    memory  = "1gb",
    walltime= "00:10:00",
    cores   = 1,
) << f"""
pixi run cat > {PARQUET}/README.md << 'EOF'
# VEP Parquet dataset

Generated by workflow.py using Ensembl release {ENSEMBL} / GRCh38.

## Structure
parquet/
    chrom=chr1/part-0.parquet
    chrom=chr2/part-0.parquet
    ...

## Query example (DuckDB)
import duckdb
duckdb.query(\"\"\"
    SELECT gene_symbol, consequence_terms, impact, hgvsp
    FROM read_parquet('parquet/**/*.parquet', hive_partitioning=true)
    WHERE chrom = 'chr17'
      AND pos BETWEEN 43000000 AND 44000000
      AND canonical = 1
    ORDER BY pos
\"\"\").df()
EOF
echo "All chromosomes complete."
"""
