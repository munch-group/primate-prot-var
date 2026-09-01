#!/usr/bin/env python3
"""
merge_loss_matrices.py
======================

Combine the per-source gene x lineage loss/presence matrices produced by
`intersect_primate_lof.py` (TOGA, Mao, the human 1KG-SV source, ...) into ONE
unified gene x lineage 0/1 matrix that drops straight into the group's geneset.py
UpSet.

Why this exists
---------------
The sources use incompatible lineage vocabularies (TOGA binomial species, Mao
common names, Human_* tiers), so the engine writes each one independently and the
docs compare them per-gene. This script makes that comparison concrete: it
outer-joins every source on the `gene` key and namespaces each source's lineage
columns (`toga:`, `mao:`, `human:`) so provenance is explicit and columns never
collide.

The "absent vs not-lost" subtlety
---------------------------------
A 0 in the merged matrix means "not flagged as lost by that source" — which
conflates "tested and intact" with "gene not in that source's universe at all".
To keep that honest we also write a companion COVERAGE MASK
(`<out-stem>.coverage.tsv`) with one `<label>:_tested` column per source (1 = the
gene is in that source's gene universe). Read the two together.

Inputs may be plain .tsv or pd-lfs parquet directories (e.g. the genome-wide
`toga_primates_species.loss_matrix.parquet/`).

Usage
-----
  python merge_loss_matrices.py \
      --matrix toga=results/recombination_genes/toga_strictL.loss_matrix.tsv \
      --matrix mao=results/recombination_genes/mao.matrix.tsv \
      --matrix human=results/recombination_genes/human.loss_matrix.tsv \
      --out results/recombination_genes/merged.loss_matrix.tsv
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from pathlib import Path

import pandas as pd

# Columns that are gene-level metadata, never a lineage 0/1 column.
META_COLS = {
    "gene", "chrom", "gene_id", "gene_type", "svtype", "consequence",
    "af", "ac", "an", "n_hom_alt", "n_het", "sv_id", "raw_id",
    "query_gene", "lineage", "status", "all_statuses",
}
# metadata we carry into the merged output, in this order (after the gene key)
CARRY_META = ["chrom", "gene_id"]


def read_matrix(path: str) -> pd.DataFrame:
    """Read a source matrix from .tsv or a pd-lfs parquet directory."""
    if path.endswith(".parquet") or os.path.isdir(path):
        try:
            from pd_lfs.parquet import read_parquet
        except Exception as e:                       # pragma: no cover
            raise RuntimeError(f"pd-lfs needed to read parquet {path}: {e}")
        return read_parquet(path)
    return pd.read_csv(path, sep="\t", dtype=str)


def split_columns(df: pd.DataFrame):
    """Return (gene_key, meta_cols_present, lineage_cols)."""
    cols = list(df.columns)
    if "gene" not in cols:
        raise ValueError(
            f"matrix has no 'gene' column (found {cols[:6]}...). Pass a "
            f"loss_matrix/matrix file, not an intersections file.")
    meta = [c for c in cols if c in META_COLS and c != "gene"]
    lineages = [c for c in cols if c != "gene" and c not in META_COLS]
    return "gene", meta, lineages


def to01(series: pd.Series, label: str, col: str) -> pd.Series:
    """Coerce a lineage column to 0/1 ints, warning if it looks non-binary."""
    num = pd.to_numeric(series, errors="coerce")
    if num.isna().any():
        bad = series[num.isna()].dropna().unique()[:5]
        print(f"  [warn] {label}:{col} has non-numeric values {list(bad)} "
              f"-> treated as 0 (is this a status_matrix, not a loss/presence one?)",
              file=sys.stderr)
    return (num.fillna(0) > 0).astype(int)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--matrix", action="append", default=[], metavar="LABEL=PATH",
                   help="A source matrix as LABEL=PATH (repeatable). LABEL becomes "
                        "the lineage-column prefix and a <LABEL>:_tested coverage column.")
    p.add_argument("--out", default="results/merged/loss_matrix.tsv",
                   help="Unified gene x lineage 0/1 matrix (default: results/merged/loss_matrix.tsv).")
    p.add_argument("--coverage-out", default=None,
                   help="Coverage-mask path (default: <out-stem>.coverage.tsv).")
    p.add_argument("--prefix-lineages", dest="prefix", action="store_true", default=True,
                   help="Prefix lineage columns with '<label>:' (default ON).")
    p.add_argument("--no-prefix-lineages", dest="prefix", action="store_false")
    p.add_argument("--universe", action="append", default=[], metavar="LABEL=PATH",
                   help="Override a source's TESTED gene universe for the coverage mask "
                        "with a gene list (one symbol per line, or a TSV with a 'gene' "
                        "column). Needed for positive-only sources (Mao, human) whose "
                        "matrix lists only FLAGGED genes although ALL coding genes were "
                        "scanned: e.g. --universe human=data/hg38_gene_cds_span.tsv so a "
                        "gene with no SV reads as tested-and-clean, not untested.")
    p.add_argument("--include-yoo", action="store_true",
                   help="Include a matrix labelled 'yoo' (excluded by default: it is "
                        "gene GAINS, not loss, and uses non-hg38 contig coordinates).")
    p.add_argument("--fill-missing", type=int, default=0,
                   help="Value for a gene absent from a source (default 0).")
    args = p.parse_args()

    if not args.matrix:
        p.error("provide at least one --matrix LABEL=PATH")

    # optional tested-universe overrides (for positive-only sources)
    universe = {}
    for spec in args.universe:
        label, path = spec.split("=", 1)
        genes = set()
        with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
            first = fh.readline().rstrip("\n").split("\t")
            gi = first.index("gene") if "gene" in first else 0
            if "gene" not in first and first and first[gi].strip():
                genes.add(first[gi].strip())          # headerless: first line is data
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > gi and f[gi].strip():
                    genes.add(f[gi].strip())
        universe[label.strip()] = genes

    lineage_frames = []      # list of (label, DataFrame indexed by gene, 0/1 lineage cols)
    coverage = {}            # label -> set(genes)
    meta_maps = {m: {} for m in CARRY_META}   # gene -> value (first source wins)
    order = []

    for spec in args.matrix:
        if "=" not in spec:
            p.error(f"--matrix must be LABEL=PATH, got {spec!r}")
        label, path = spec.split("=", 1)
        label, path = label.strip(), path.strip()
        if label == "yoo" and not args.include_yoo:
            print(f"  [skip] '{label}' excluded by default (gains, non-hg38). "
                  f"Use --include-yoo to add it.", file=sys.stderr)
            continue
        df = read_matrix(path)
        gkey, meta, lineages = split_columns(df)
        df = df.copy()
        df[gkey] = df[gkey].astype(str)
        df = df.drop_duplicates(subset=[gkey])
        coverage[label] = set(df[gkey])
        # carry metadata (first source to define a gene's value wins)
        for m in CARRY_META:
            if m in df.columns:
                for g, v in zip(df[gkey], df[m].astype(str)):
                    if g not in meta_maps[m] and v and v.lower() != "nan":
                        meta_maps[m][g] = v
        # 0/1 lineage block, indexed by gene, prefixed (build in one shot to
        # avoid DataFrame fragmentation when a source has many lineages, e.g. TOGA)
        data = {c: to01(df[c], label, c).values for c in lineages}
        block = pd.DataFrame(data, index=pd.Index(df[gkey].values, name=gkey))
        if args.prefix:
            block.columns = [f"{label}:{c}" for c in lineages]
        lineage_frames.append((label, block))
        order.append(label)
        print(f"  [read] {label}: {len(df):,} genes, {len(lineages)} lineage cols "
              f"<- {path}", file=sys.stderr)

    if not lineage_frames:
        p.error("no matrices left to merge")

    # union of all genes
    all_genes = sorted(set().union(*[set(b.index) for _, b in lineage_frames]))
    merged = pd.DataFrame(index=pd.Index(all_genes, name="gene"))
    for _, block in lineage_frames:
        merged = merged.join(block, how="left")
    merged = merged.fillna(args.fill_missing).astype(int)

    # leading metadata
    out = merged.reset_index()
    for m in reversed(CARRY_META):                # insert chrom then gene_id after 'gene'
        out.insert(1, m, out["gene"].map(meta_maps[m]).fillna(""))

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    # coverage mask
    cov_path = args.coverage_out or str(Path(args.out).with_suffix("")) + ".coverage.tsv"
    cov = pd.DataFrame({"gene": all_genes})
    for label in order:
        tested = universe.get(label, coverage[label])   # universe override if given
        cov[f"{label}:_tested"] = cov["gene"].map(lambda g, T=tested: int(g in T))
    cov.to_csv(cov_path, sep="\t", index=False)

    n_lin = sum(len(b.columns) for _, b in lineage_frames)
    print(f"[merge] sources: {', '.join(order)}", file=sys.stderr)
    print(f"[merge] unified matrix: {len(all_genes):,} genes x {n_lin} lineage cols "
          f"(+{len(CARRY_META)} meta) -> {args.out}", file=sys.stderr)
    print(f"[merge] coverage mask -> {cov_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
