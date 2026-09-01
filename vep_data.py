"""
vep_data.py
===========
Downloads VEP annotation data from Ensembl and demonstrates fast
random-access query patterns using Parquet+DuckDB and Tabix.

Two data sources:
  1. REST API  – small batches (≤200 variants), good for exploration
  2. FTP bulk  – full pre-computed VEP VCFs, good for production

Dependencies:
    pip install requests pandas pyarrow duckdb pysam tqdm
    (pysam is optional – only needed for tabix queries)
"""

import os
import json
import time
import gzip
import shutil
import subprocess
import requests
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import duckdb
from pathlib import Path
from tqdm import tqdm

# ── directories ────────────────────────────────────────────────────────────────
DATA_DIR    = Path("vep_data")
PARQUET_DIR = DATA_DIR / "parquet"
VCF_DIR     = DATA_DIR / "vcf"
for d in (PARQUET_DIR, VCF_DIR):
    d.mkdir(parents=True, exist_ok=True)

ENSEMBL_REST = "https://rest.ensembl.org"
ENSEMBL_FTP  = "https://ftp.ensembl.org/pub"

# ── 1. REST API – fetch VEP for a list of variants ────────────────────────────

EXAMPLE_VARIANTS = [
    # HGVS notation (works for SNPs, indels, etc.)
    "9:g.22125504G>C",   # CDKN2A
    "7:g.117548628T>A",  # CFTR
    "17:g.43094692G>A",  # BRCA1
    "13:g.32315474A>T",  # BRCA2
    "12:g.25398284C>A",  # KRAS G12V
    "17:g.7674220C>T",   # TP53
    "3:g.178936091A>G",  # PIK3CA
    "10:g.89692905A>G",  # PTEN
    "7:g.140453136A>T",  # BRAF V600E
    "2:g.29443613C>T",   # ALK
]


def fetch_vep_rest(variants: list[str], assembly: str = "GRCh38") -> pd.DataFrame:
    """
    POST up to 200 variants to the Ensembl VEP REST endpoint.
    Returns a flat DataFrame with one row per transcript consequence.
    """
    server = ENSEMBL_REST
    endpoint = "/vep/human/hgvs"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    payload = json.dumps({"hgvs_notations": variants})

    print(f"[REST] Fetching VEP for {len(variants)} variants …")
    r = requests.post(f"{server}{endpoint}", headers=headers, data=payload, timeout=60)
    r.raise_for_status()
    raw = r.json()

    rows = []
    for entry in raw:
        variant_id = entry.get("id", entry.get("input", ""))
        chrom      = entry.get("seq_region_name", "")
        pos        = entry.get("start", None)
        ref        = entry.get("allele_string", "").split("/")[0] if "/" in entry.get("allele_string","") else ""
        alt        = entry.get("allele_string", "").split("/")[-1] if "/" in entry.get("allele_string","") else ""

        # one row per transcript consequence
        for tc in entry.get("transcript_consequences", []):
            rows.append({
                "variant_id"          : variant_id,
                "chrom"               : chrom,
                "pos"                 : pos,
                "ref"                 : ref,
                "alt"                 : alt,
                "gene_id"             : tc.get("gene_id"),
                "gene_symbol"         : tc.get("gene_symbol"),
                "transcript_id"       : tc.get("transcript_id"),
                "biotype"             : tc.get("biotype"),
                "consequence_terms"   : ",".join(tc.get("consequence_terms", [])),
                "impact"              : tc.get("impact"),
                "hgvsc"               : tc.get("hgvsc"),
                "hgvsp"               : tc.get("hgvsp"),
                "sift_score"          : tc.get("sift_score"),
                "sift_prediction"     : tc.get("sift_prediction"),
                "polyphen_score"      : tc.get("polyphen_score"),
                "polyphen_prediction" : tc.get("polyphen_prediction"),
                "cadd_phred"          : entry.get("cadd_phred"),      # top-level
                "af_gnomade"          : tc.get("gnomade_af"),
                "canonical"           : tc.get("canonical", 0),
                "strand"              : tc.get("strand"),
            })

    df = pd.DataFrame(rows)
    df["chrom"] = df["chrom"].astype(str)
    df["pos"]   = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    print(f"[REST] Got {len(df)} transcript consequence rows for {len(raw)} variants.")
    return df


# ── 2. FTP bulk download – pre-computed VEP VCF ───────────────────────────────

def ftp_vcf_url(release: int = 111, species: str = "homo_sapiens") -> str:
    return f"{ENSEMBL_FTP}/release-{release}/variation/vcf/{species}/"


def download_ftp_vcf(
    release: int = 111,
    species: str = "homo_sapiens",
    chrom: str = "1",
    dest_dir: Path = VCF_DIR,
) -> Path:
    """
    Download a single-chromosome VEP VCF + tabix index from Ensembl FTP.
    For chr1 of Homo sapiens this is ~3 GB; use a small chromosome for testing.
    """
    base  = ftp_vcf_url(release, species)
    fname = f"{species.capitalize().replace('_',' ').title().replace(' ','_')}_chr{chrom}.vcf.gz"
    url   = base + fname
    dest  = dest_dir / fname

    if dest.exists():
        print(f"[FTP] {dest.name} already exists, skipping download.")
        return dest

    print(f"[FTP] Downloading {url} …")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with open(dest, "wb") as f, tqdm(total=total, unit="B", unit_scale=True) as bar:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
                bar.update(len(chunk))

    # also grab the .tbi index
    tbi_url  = url + ".tbi"
    tbi_dest = Path(str(dest) + ".tbi")
    print(f"[FTP] Downloading tabix index …")
    r2 = requests.get(tbi_url, timeout=60)
    r2.raise_for_status()
    tbi_dest.write_bytes(r2.content)

    return dest


# ── 3. Save to Parquet ─────────────────────────────────────────────────────────

def save_parquet(df: pd.DataFrame, out_dir: Path = PARQUET_DIR) -> Path:
    """
    Write DataFrame as a Parquet dataset partitioned by chromosome.
    Sorted by (chrom, pos) so DuckDB can exploit min/max statistics.
    """
    df_sorted = df.sort_values(["chrom", "pos"])
    table = pa.Table.from_pandas(df_sorted, preserve_index=False)
    pq.write_to_dataset(
        table,
        root_path=str(out_dir),
        partition_cols=["chrom"],
        compression="zstd",
        use_dictionary=True,
    )
    print(f"[Parquet] Written to {out_dir}/")
    return out_dir


# ── 4. Save to tabix-indexed VCF ──────────────────────────────────────────────

def save_tabix_vcf(df: pd.DataFrame, out_dir: Path = VCF_DIR) -> Path:
    """
    Write a minimal bgzipped VCF and create a tabix index.
    Requires bgzip + tabix (from htslib) on PATH.
    """
    raw_vcf = out_dir / "vep_annotations.vcf"

    with open(raw_vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">\n')
        fh.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence">\n')
        fh.write('##INFO=<ID=IMPACT,Number=1,Type=String,Description="VEP impact">\n')
        fh.write('##INFO=<ID=SIFT,Number=1,Type=Float,Description="SIFT score">\n')
        fh.write('##INFO=<ID=PP2,Number=1,Type=Float,Description="PolyPhen2 score">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        canon = df[df["canonical"] == 1].drop_duplicates(["chrom", "pos", "ref", "alt"])
        for _, row in canon.sort_values(["chrom", "pos"]).iterrows():
            gene   = row["gene_symbol"] or "."
            csq    = row["consequence_terms"] or "."
            impact = row["impact"] or "."
            sift   = row["sift_score"] if pd.notna(row.get("sift_score")) else "."
            pp2    = row["polyphen_score"] if pd.notna(row.get("polyphen_score")) else "."
            info   = f"GENE={gene};CSQ={csq};IMPACT={impact};SIFT={sift};PP2={pp2}"
            fh.write(
                f"{row['chrom']}\t{row['pos']}\t{row['variant_id']}\t"
                f"{row['ref']}\t{row['alt']}\t.\t.\t{info}\n"
            )

    bgz = Path(str(raw_vcf) + ".gz")
    subprocess.run(["bgzip", "-f", str(raw_vcf)], check=True)
    subprocess.run(["tabix", "-p", "vcf", str(bgz)], check=True)
    print(f"[Tabix] Written {bgz} + .tbi")
    return bgz


# ── 5. Query examples ─────────────────────────────────────────────────────────

def demo_duckdb_queries(parquet_dir: Path) -> None:
    """
    Demonstrates common VEP query patterns using DuckDB against Parquet files.
    DuckDB reads only the relevant row-groups — no server required.
    """
    con = duckdb.connect()
    glob = str(parquet_dir / "**" / "*.parquet")

    print("\n" + "═" * 60)
    print("DuckDB / Parquet queries")
    print("═" * 60)

    # ── Q1: positional lookup (single variant) ────────────────────
    print("\n── Q1: All consequences for a single genomic position ──")
    q1 = con.execute(f"""
        SELECT gene_symbol, transcript_id, consequence_terms, impact, hgvsp
        FROM read_parquet('{glob}', hive_partitioning=true)
        WHERE chrom = '17' AND pos = 43094692
        ORDER BY canonical DESC, impact
    """).df()
    print(q1.to_string(index=False))

    # ── Q2: region query ──────────────────────────────────────────
    print("\n── Q2: High-impact variants in a genomic window ──")
    q2 = con.execute(f"""
        SELECT chrom, pos, ref, alt, gene_symbol, consequence_terms, impact
        FROM read_parquet('{glob}', hive_partitioning=true)
        WHERE chrom = '17'
          AND pos BETWEEN 43000000 AND 44000000
          AND impact IN ('HIGH', 'MODERATE')
        ORDER BY pos
    """).df()
    print(q2.to_string(index=False))

    # ── Q3: filter by gene ────────────────────────────────────────
    print("\n── Q3: Canonical transcripts for BRCA1 ──")
    q3 = con.execute(f"""
        SELECT pos, ref, alt, consequence_terms, impact, hgvsc, hgvsp,
               sift_prediction, polyphen_prediction
        FROM read_parquet('{glob}', hive_partitioning=true)
        WHERE gene_symbol = 'BRCA1' AND canonical = 1
        ORDER BY pos
    """).df()
    print(q3.to_string(index=False))

    # ── Q4: filter by consequence type ───────────────────────────
    print("\n── Q4: All stop_gained / frameshift variants ──")
    q4 = con.execute(f"""
        SELECT chrom, pos, gene_symbol, consequence_terms, hgvsp
        FROM read_parquet('{glob}', hive_partitioning=true)
        WHERE consequence_terms LIKE '%stop_gained%'
           OR consequence_terms LIKE '%frameshift%'
    """).df()
    print(q4.to_string(index=False))

    # ── Q5: SIFT / PolyPhen damaging filter ──────────────────────
    print("\n── Q5: Predicted deleterious variants (SIFT + PolyPhen) ──")
    q5 = con.execute(f"""
        SELECT chrom, pos, gene_symbol, consequence_terms,
               sift_score, sift_prediction,
               polyphen_score, polyphen_prediction
        FROM read_parquet('{glob}', hive_partitioning=true)
        WHERE sift_prediction    LIKE '%deleterious%'
          AND polyphen_prediction LIKE '%damaging%'
          AND canonical = 1
        ORDER BY sift_score ASC
    """).df()
    print(q5.to_string(index=False))

    # ── Q6: aggregate – impact counts per gene ────────────────────
    print("\n── Q6: Variant impact counts per gene ──")
    q6 = con.execute(f"""
        SELECT gene_symbol,
               COUNT(DISTINCT pos)                                     AS n_variants,
               SUM(CASE WHEN impact='HIGH'     THEN 1 ELSE 0 END)     AS high,
               SUM(CASE WHEN impact='MODERATE' THEN 1 ELSE 0 END)     AS moderate,
               SUM(CASE WHEN impact='LOW'      THEN 1 ELSE 0 END)     AS low
        FROM read_parquet('{glob}', hive_partitioning=true)
        WHERE canonical = 1
        GROUP BY gene_symbol
        ORDER BY high DESC, moderate DESC
    """).df()
    print(q6.to_string(index=False))

    # ── Q7: export filtered subset ────────────────────────────────
    print("\n── Q7: Export BRCA1 + BRCA2 HIGH/MODERATE variants to CSV ──")
    out_csv = DATA_DIR / "brca_variants.csv"
    con.execute(f"""
        COPY (
            SELECT chrom, pos, ref, alt, gene_symbol,
                   consequence_terms, impact, hgvsc, hgvsp,
                   sift_score, polyphen_score
            FROM read_parquet('{glob}', hive_partitioning=true)
            WHERE gene_symbol IN ('BRCA1', 'BRCA2')
              AND impact IN ('HIGH', 'MODERATE')
              AND canonical = 1
            ORDER BY chrom, pos
        ) TO '{out_csv}' (HEADER, DELIMITER ',')
    """)
    print(f"    → Saved to {out_csv}")


def demo_tabix_queries(bgz_path: Path) -> None:
    """
    Demonstrates tabix (pysam) positional queries — fastest for point/region lookups.
    """
    try:
        import pysam
    except ImportError:
        print("\n[Tabix] pysam not installed — skipping tabix demo. pip install pysam")
        return

    print("\n" + "═" * 60)
    print("Tabix / pysam queries")
    print("═" * 60)

    tbx = pysam.TabixFile(str(bgz_path))

    # ── T1: fetch a specific position ─────────────────────────────
    print("\n── T1: Point lookup chr17:43094692 (BRCA1) ──")
    for row in tbx.fetch("17", 43094691, 43094692):
        cols = row.split("\t")
        print(f"  {cols[0]}:{cols[1]}  {cols[3]}>{cols[4]}  INFO: {cols[7]}")

    # ── T2: region query ──────────────────────────────────────────
    print("\n── T2: All variants in chr7:140,000,000–140,600,000 (BRAF region) ──")
    for row in tbx.fetch("7", 140_000_000, 140_600_000):
        cols = row.split("\t")
        print(f"  {cols[0]}:{cols[1]}  {cols[2]}  {cols[7][:80]}")

    # ── T3: iterate whole chromosome ─────────────────────────────
    print("\n── T3: Count variants on chr9 ──")
    n = sum(1 for _ in tbx.fetch("9"))
    print(f"  chr9 variant count: {n}")

    tbx.close()


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("VEP data download + query demo")
    print("=" * 60)

    # ── Step 1: fetch via REST API ────────────────────────────────
    df = fetch_vep_rest(EXAMPLE_VARIANTS)

    # ── Step 2: persist ───────────────────────────────────────────
    save_parquet(df)

    bgz = None
    try:
        bgz = save_tabix_vcf(df)
    except FileNotFoundError:
        print("[Tabix] bgzip/tabix not found on PATH — skipping VCF output.")
        print("        Install htslib: conda install -c bioconda htslib")

    # ── Step 3: query demos ───────────────────────────────────────
    demo_duckdb_queries(PARQUET_DIR)

    # if bgz and bgz.exists():
    #     demo_tabix_queries(bgz)


    # ── Step 4: show FTP bulk download commands ───────────────────
    print("\n" + "═" * 60)
    print("Bulk FTP download (production use)")
    print("═" * 60)
    print("""
For full genome, use rsync (resumable, much faster than wget):

  # All human pre-computed VEP VCFs + tabix indices
  rsync -avP --include='*.vcf.gz' --include='*.vcf.gz.tbi' \\
    rsync://ftp.ensembl.org/ensembl/pub/current_variation/vcf/homo_sapiens/ \\
    ./vep_data/vcf/

  # Or a single chromosome (e.g. chr22, smallest autosome)
  wget https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/\\
Homo_sapiens_chr22.vcf.gz
  wget https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/\\
Homo_sapiens_chr22.vcf.gz.tbi

  # Convert VCF to Parquet after download (requires cyvcf2):
  #   pip install cyvcf2 pyarrow
  #   python -c "
  #     import cyvcf2, pandas as pd, pyarrow.parquet as pq, pyarrow as pa
  #     vcf = cyvcf2.VCF('Homo_sapiens_chr22.vcf.gz')
  #     rows = [{'chrom': v.CHROM, 'pos': v.POS, 'ref': v.REF,
  #              'alt': ','.join(v.ALT), 'csq': v.INFO.get('CSQ','')}
  #             for v in vcf]
  #     pq.write_table(pa.Table.from_pylist(rows), 'chr22.parquet',
  #                    compression='zstd')
  #   "
""")

    print("\n✓ Done. Files written to:", DATA_DIR.resolve())


if __name__ == "__main__":
    main()

