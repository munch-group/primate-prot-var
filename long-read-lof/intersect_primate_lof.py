#!/usr/bin/env python3
"""
intersect_primate_lof.py
========================

Cross-reference a query gene list (e.g. meiosis_repair_genes.py) against
lineage-specific loss-of-function / structural-variant / gene-status data from
long-read primate genome resources, and emit per-branch intersections plus an
UpSet-ready matrix.

Two input modes
---------------
1. --table  : a single supplementary workbook / TSV (Mao 2024 Data S2, Yoo 2025
              Suppl. Table VIII.34, or any gene+lineage[+status] table).
2. --toga-dir : a downloaded TOGA tree (…/human_hg38_reference/Primates/). The
              walker finds every */loss_summ_data.tsv(.gz), stamps the species /
              assembly from the directory name, and builds a tidy long table.
              This is the genome-wide "all genes, all species" path.

Outputs (status-aware whenever a status is available — always true for --toga-dir):
  <out>.intersections   one row per (gene, lineage): chrom, status, all_statuses
  <out>.status_matrix   gene × lineage -> most-severe status code
  <out>.loss_matrix     gene × lineage -> 0/1 (1 = a "loss" status); UpSet-ready
  <out>.toga_long       (--toga-dir only) the merged long table, for reuse
  (without any status, --table mode falls back to a single <out>.matrix 0/1)

  --meta-cols (--table only) carries extra per-gene source columns through to
  every output as leading columns after `chrom` (e.g. a stable Ensembl id and a
  biotype, so symbol-less GENCODE clone names like AC007993 stay joinable /
  filterable): --meta-cols GENE_ID:gene_id,GENE_TYPE:gene_type.

  --toga-dir writes each output as a pd-lfs multi-file parquet *directory*
  (<out>.<name>.parquet/ of part-*.parquet + _manifest.json) so the genome-wide
  tables stay under GitHub's file-size limit; --table writes plain <name>.tsv.
  Either form is reusable via --table (read_parquet handles the directory).

Output location
---------------
A bare --out label (no '/') nests every output in a per-query folder named after
the query gene list: --query recombination_genes.txt --out toga writes
results/recombination_genes/toga.* — so all cuts of one curated set live together.
A value containing '/' (or an absolute path) is used verbatim and flat, e.g.
--out results/toga_primates_species for the genome-wide walk. --results-dir
changes the base directory (default: results).

No network calls. pandas (+ openpyxl for .xlsx, pd-lfs for parquet). gz inputs
read transparently.

Author: scaffold for Kasper / Munch group
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

# pd-lfs (munch-group.org/pd-lfs): write_parquet shards a DataFrame across
# part-*.parquet files (<50 MB each) with a _manifest.json index, so even the
# genome-wide TOGA tables fit under GitHub's file-size limit. read_parquet
# reverses it (local path or https URL) and restores the original dtypes.
try:
    from pd_lfs.parquet import read_parquet, write_parquet
except ImportError as _pd_lfs_err:  # surfaced only if a parquet path is used
    read_parquet = write_parquet = None
    _PD_LFS_IMPORT_ERROR = _pd_lfs_err

SPLIT_RE = re.compile(r"[,;/|\s]+")
VERSION_RE = re.compile(r"\.\d+$")

# TOGA-style severity ranking for collapsing multiple calls per (gene, lineage).
STATUS_RANK = {
    "I": 0, "INTACT": 0,
    "PI": 1, "PARTIALLY_INTACT": 1, "PARTIAL_INTACT": 1,
    "UL": 2, "UNCERTAIN_LOSS": 2,
    "PG": 2, "PARALOG": 2, "PARALOGOUS_PROJECTION": 2,
    "PM": 3, "PARTIAL_MISSING": 3,
    "L": 4, "LOST": 4, "LOSS": 4, "PSEUDOGENE": 4,
    "M": 5, "MISSING": 5,
    # Human SV-LoF tiers (annotate_sv_lof.py): low- vs high-confidence loss-of-function.
    # Registered only so most_severe() prefers HC over LC; these are NOT in
    # DEFAULT_LOSS_STATUSES — the human source must pass --loss-status explicitly.
    "LC_LOF": 2, "HC_LOF": 4,
}
DEFAULT_LOSS_STATUSES = {
    "UL", "UNCERTAIN_LOSS", "PG", "PARALOG", "PARALOGOUS_PROJECTION",
    "PM", "PARTIAL_MISSING", "L", "LOST", "LOSS", "PSEUDOGENE", "M", "MISSING",
}


# --------------------------------------------------------------------------- #
# Normalisation helpers
# --------------------------------------------------------------------------- #
def norm(sym) -> str:
    if sym is None:
        return ""
    s = str(sym).strip().upper()
    return VERSION_RE.sub("", s)


def norm_status(s) -> str:
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return ""
    return str(s).strip().upper().replace(" ", "_").replace("-", "_")


def norm_chrom(c) -> str:
    if c is None or (isinstance(c, float) and pd.isna(c)):
        return ""
    s = str(c).strip().lower()
    if s.startswith("chr"):
        s = s[3:]
    return s.upper()


def explode_cell(cell) -> list[str]:
    if pd.isna(cell):
        return []
    parts = SPLIT_RE.split(str(cell))
    return [norm(p) for p in parts if norm(p)]


def logical_suffix(path: Path) -> str:
    """Suffix ignoring a trailing .gz (so foo.tsv.gz -> .tsv)."""
    name = path.name[:-3] if path.name.endswith(".gz") else path.name
    return Path(name).suffix.lower()


def most_severe(statuses: set[str]) -> str:
    known = [s for s in statuses if s]
    if not known:
        return ""
    return max(known, key=lambda s: STATUS_RANK.get(s, 1))


def canon(sym, amap: dict[str, str]) -> str:
    n = norm(sym)
    return amap.get(n, n)


# --------------------------------------------------------------------------- #
# Loading (gz-aware; pandas infers compression from the .gz extension)
# --------------------------------------------------------------------------- #
def load_alias_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    amap: dict[str, str] = {}
    df = pd.read_csv(path, sep="\t", header=None, dtype=str).fillna("")
    for alias, canon_ in zip(df[0], df[1]):
        a, c = norm(alias), norm(canon_)
        if a and c:
            amap[a] = c
    return amap


def read_query_genes(path: Path, column: str | None) -> set[str]:
    suf = logical_suffix(path)
    if suf in {".csv", ".tsv", ".txt"}:
        if suf == ".txt" and column is None:
            genes = [ln.strip() for ln in path.read_text().splitlines()]
            return {norm(g) for g in genes if g and not g.startswith("#")}
        sep = "," if suf == ".csv" else "\t"
        df = pd.read_csv(path, sep=sep, dtype=str)
        if column is None:
            column = df.columns[0]
        return {norm(g) for g in df[column].dropna()}
    raise ValueError(f"Unsupported query file type: {suf}")


def read_table(path: Path, sheet) -> dict[str, pd.DataFrame]:
    suf = logical_suffix(path)
    if path.is_dir() or suf == ".parquet":
        # a pd-lfs dataset directory produced by --toga-dir mode
        if read_parquet is None:
            raise RuntimeError(
                "pd-lfs is required to read a parquet --table but could not be "
                f"imported ({_PD_LFS_IMPORT_ERROR}). Install it: pip install pd-lfs."
            )
        return {path.stem: read_parquet(str(path))}
    if suf in {".xlsx", ".xls", ".xlsm"}:
        xl = pd.ExcelFile(path)
        if sheet is None:
            return {s: xl.parse(s, dtype=str) for s in xl.sheet_names}
        return {str(sheet): xl.parse(sheet, dtype=str)}
    sep = "," if suf == ".csv" else "\t"
    return {path.stem: pd.read_csv(path, sep=sep, dtype=str)}


def load_annotation(path: Path | None, cols: str):
    """
    Optional reference annotation mapping the loss_summ identifier -> (symbol, chrom).
    `cols` is a comma list naming the id, symbol, chrom columns, e.g. 'transcript,gene,chrom'.
    Symbol and/or chrom may be omitted with '-' (e.g. 'id,-,chrom').
    Works on TSV/CSV/BED-like files (gz ok); reads with a header.
    """
    if path is None:
        return {}, {}
    names = [c.strip() for c in cols.split(",")]
    while len(names) < 3:
        names.append("-")
    id_c, sym_c, chr_c = names[:3]
    sep = "," if logical_suffix(path) == ".csv" else "\t"
    df = pd.read_csv(path, sep=sep, dtype=str)
    id2sym, id2chrom = {}, {}
    for _, row in df.iterrows():
        rid = norm(row.get(id_c))
        if not rid:
            continue
        if sym_c != "-" and sym_c in df.columns and not pd.isna(row.get(sym_c)):
            id2sym[rid] = norm(row[sym_c])
        if chr_c != "-" and chr_c in df.columns and not pd.isna(row.get(chr_c)):
            id2chrom[rid] = norm_chrom(row[chr_c])
    return id2sym, id2chrom


# --------------------------------------------------------------------------- #
# TOGA tree walker
# --------------------------------------------------------------------------- #
def parse_toga_dirname(name: str):
    """'Pan_troglodytes__chimpanzee__HLpanTroT' -> ('PAN_TROGLODYTES','chimpanzee','HLpanTroT').
    Robust to '__Primates__' or '__-__' placeholders and trinomial species names."""
    parts = name.split("__")
    species = norm(parts[0]) if parts else norm(name)
    assembly = parts[-1] if len(parts) > 1 else ""
    common = parts[1] if len(parts) > 2 else ""
    return species, common, assembly


def read_loss_summ(path: Path, level_keep: str):
    """loss_summ_data.tsv(.gz): headerless 3 cols = entry_type, identifier, status.
    Returns list of (identifier, status) filtered to level_keep ('GENE' default; 'ANY' = no filter)."""
    df = pd.read_csv(path, sep="\t", header=None, dtype=str,
                     names=["level", "id", "status"], usecols=[0, 1, 2])
    out = []
    lk = level_keep.upper()
    for _, r in df.iterrows():
        if lk != "ANY" and str(r["level"]).strip().upper() != lk:
            continue
        rid, st = str(r["id"]).strip(), norm_status(r["status"])
        if rid:
            out.append((rid, st))
    return out


def walk_toga_dir(root: Path, label: str, level_keep: str):
    """Yield (raw_id, lineage, status, assembly_dirname) across the tree.
    label: 'species' (genus_species), 'assembly' (full dir name), 'assembly_id' (last token)."""
    files = sorted(root.glob("*/loss_summ_data.tsv.gz")) + \
            sorted(root.glob("*/loss_summ_data.tsv"))
    if not files:
        # also allow root itself being a single assembly dir
        files = sorted(root.glob("loss_summ_data.tsv*"))
    seen = set()
    for f in files:
        if f in seen:
            continue
        seen.add(f)
        dirname = f.parent.name
        species, _common, asm = parse_toga_dirname(dirname)
        if label == "species":
            lineage = species
        elif label == "assembly_id":
            lineage = asm or dirname
        else:
            lineage = dirname
        for rid, st in read_loss_summ(f, level_keep):
            yield rid, lineage, st, dirname


# --------------------------------------------------------------------------- #
# Core accumulation for --table mode
# --------------------------------------------------------------------------- #
def accumulate_table(df, gene_col, lineage_col, status_col, chrom_col,
                     chrom_keep, amap, store, chrom_of, meta_specs, meta_of):
    for _, row in df.iterrows():
        if chrom_col and chrom_keep is not None:
            if norm_chrom(row.get(chrom_col)) not in chrom_keep:
                continue
        genes = explode_cell(row.get(gene_col))
        if not genes:
            continue
        lineage = (str(row[lineage_col]).strip()
                   if lineage_col and lineage_col in df.columns and not pd.isna(row.get(lineage_col))
                   else "ALL")
        status = norm_status(row.get(status_col)) if status_col else "HIT"
        ch = norm_chrom(row.get(chrom_col)) if chrom_col else ""
        meta_vals = []
        for src, out in meta_specs:
            v = row.get(src)
            meta_vals.append((out, "" if pd.isna(v) else str(v).strip()))
        for g in genes:
            cg = canon(g, amap)
            store.setdefault(cg, {}).setdefault(lineage, set()).add(status)
            if ch:
                chrom_of[cg] = ch
            for out, v in meta_vals:
                if v:
                    meta_of.setdefault(cg, {}).setdefault(out, v)  # first non-empty wins


# --------------------------------------------------------------------------- #
# Output writer (pd-lfs multi-file parquet for --toga-dir, plain TSV otherwise)
# --------------------------------------------------------------------------- #
def write_df(df, out_base: str, name: str, as_parquet: bool) -> str:
    """Write `df` and return the path written.

    as_parquet=True  -> pd-lfs dataset directory `<out_base>.<name>.parquet/`
                        (part-*.parquet shards + _manifest.json; git-friendly).
    as_parquet=False -> plain `<out_base>.<name>.tsv`.
    """
    if as_parquet:
        if write_parquet is None:
            raise RuntimeError(
                "pd-lfs is required to write parquet outputs but could not be "
                f"imported ({_PD_LFS_IMPORT_ERROR}). Install it (pip install pd-lfs) "
                "or run in --table mode (which writes plain .tsv)."
            )
        # pd-lfs infers the Arrow schema from an empty slice (df.iloc[:0]); a
        # plain object column infers as `null` there and then fails to write its
        # real string values. Give string columns a concrete pandas "string"
        # dtype so high-cardinality ones (gene, raw_id) round-trip correctly.
        obj_cols = df.select_dtypes(include="object").columns
        if len(obj_cols):
            df[obj_cols] = df[obj_cols].astype("string")
        path = f"{out_base}.{name}.parquet"
        write_parquet(df, path)
        return path + "/"
    path = f"{out_base}.{name}.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


def parse_meta_cols(spec):
    """Parse a --meta-cols spec into an ordered list of (src, out) pairs.

    'GENE_ID:gene_id,GENE_TYPE' -> [('GENE_ID','gene_id'), ('GENE_TYPE','GENE_TYPE')].
    Whitespace is trimmed; empty items are dropped; 'SRC' (no ':') keeps the name.
    """
    pairs = []
    for item in (spec or "").split(","):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            src, name = item.split(":", 1)
            src, name = src.strip(), name.strip()
            pairs.append((src, name or src))
        else:
            pairs.append((item, item))
    return pairs


def resolve_out_base(out, query_path: Path, results_dir: str, default_label: str) -> str:
    """Compute the output path prefix shared by every <prefix>.<name> file.

    A bare --out label (no '/') nests outputs in a per-query folder:
        <results_dir>/<query-stem>/<label>     e.g. results/recombination_genes/toga
    so all files derived from one curated query gene list collect in one folder.
    A value containing '/' (or an absolute path) is used verbatim, flat — this
    preserves the genome-wide convention (`--out results/toga_primates_species`)
    and any explicit path the caller gives. The parent dir is created if needed.
    """
    label = out if out is not None else default_label
    if out is not None and ("/" in out or Path(out).is_absolute()):
        out_base = Path(out)                                     # explicit path, flat
    else:
        out_base = Path(results_dir) / query_path.stem / label   # per-query folder
    out_base.parent.mkdir(parents=True, exist_ok=True)
    return str(out_base)


# --------------------------------------------------------------------------- #
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--query", required=True, type=Path)
    p.add_argument("--query-col", default=None)
    # --- single-table mode ---
    p.add_argument("--table", default=None, type=Path)
    p.add_argument("--sheet", default=None)
    p.add_argument("--gene-col", default=None)
    p.add_argument("--lineage-col", default=None)
    p.add_argument("--status-col", default=None)
    p.add_argument("--chrom-col", default=None)
    p.add_argument("--meta-cols", default=None,
                   help="(--table only) Comma list of extra source columns to carry "
                        "through to every output as per-gene metadata, placed right "
                        "after `chrom`. Use SRC to keep the name or SRC:out to rename, "
                        "e.g. 'GENE_ID:gene_id,GENE_TYPE:gene_type'. Lets symbol-less "
                        "GENCODE clone names (AC007993) stay joinable on a stable id "
                        "and filterable by biotype. First non-empty value per gene wins.")
    p.add_argument("--list-columns", action="store_true")
    # --- TOGA tree mode ---
    p.add_argument("--toga-dir", default=None, type=Path,
                   help="Path to a downloaded …/Primates/ tree (walks */loss_summ_data.tsv.gz).")
    p.add_argument("--toga-label", default="species",
                   choices=["species", "assembly", "assembly_id"],
                   help="Lineage granularity. 'species' collapses multiple assemblies "
                        "of the same species (most-severe).")
    p.add_argument("--toga-level", default="GENE",
                   help="Which loss_summ entry type to keep: GENE (default), TRANSCRIPT, "
                        "PROJECTION, or ANY.")
    # --- shared ---
    p.add_argument("--annotation", default=None, type=Path,
                   help="Reference table mapping the loss_summ id -> symbol/chrom "
                        "(needed for --chrom and for symbol matching in --toga-dir mode).")
    p.add_argument("--annotation-cols", default="id,gene,chrom",
                   help="Comma list naming id,symbol,chrom columns in --annotation. "
                        "Use '-' to skip one, e.g. 'transcript,gene,chrom' or 'id,-,chrom'.")
    p.add_argument("--chrom", default=None,
                   help="Comma-separated chromosomes to keep (e.g. 'chrX' or 'X,Y').")
    p.add_argument("--loss-status", default=None,
                   help="Comma-separated statuses counted as 'loss' (overrides TOGA default).")
    p.add_argument("--alias", default=None, type=Path)
    p.add_argument("--out", default=None,
                   help="Output name or path prefix. A bare label (no '/') nests "
                        "outputs under <results-dir>/<query-stem>/<label>.* so a "
                        "curated query like recombination_genes.txt lands in "
                        "results/recombination_genes/. A value containing '/' (or an "
                        "absolute path) is used verbatim, flat (the genome-wide "
                        "convention, e.g. 'results/toga_primates_species'). "
                        "Default label: 'toga' for --toga-dir, else the --table stem.")
    p.add_argument("--results-dir", default="results",
                   help="Base directory for the per-query output folder "
                        "(default: results); only used when --out is a bare label.")
    args = p.parse_args()

    if not args.table and not args.toga_dir:
        p.error("provide either --table or --toga-dir.")
    if args.table and args.toga_dir:
        p.error("--table and --toga-dir are mutually exclusive.")

    # --list-columns only meaningful for --table
    if args.list_columns:
        if not args.table:
            p.error("--list-columns applies to --table mode.")
        for name, df in read_table(args.table, args.sheet).items():
            print(f"\n# sheet: {name}  (rows={len(df)})")
            for c in df.columns:
                print(f"    {c!r}")
        return 0

    chrom_keep = {norm_chrom(c) for c in args.chrom.split(",")} if args.chrom else None
    loss_set = ({norm_status(s) for s in args.loss_status.split(",")}
                if args.loss_status else DEFAULT_LOSS_STATUSES)
    amap = load_alias_map(args.alias)
    id2sym, id2chrom = load_annotation(args.annotation, args.annotation_cols)
    query = {canon(g, amap) for g in read_query_genes(args.query, args.query_col)}

    default_label = "toga" if args.toga_dir else (args.table.stem if args.table else "out")
    out_base = resolve_out_base(args.out, args.query, args.results_dir, default_label)

    # --meta-cols is a --table concept (carry extra source columns through); the
    # TOGA walk has a fixed schema, so it ignores it.
    if args.meta_cols and args.toga_dir:
        print("  [note] --meta-cols ignored in --toga-dir mode", file=sys.stderr)
    meta_specs = parse_meta_cols(args.meta_cols) if args.table else []
    meta_names = [out for _, out in meta_specs]

    store: dict[str, dict[str, set[str]]] = {}
    chrom_of: dict[str, str] = {}
    meta_of: dict[str, dict[str, str]] = {}

    if args.toga_dir:
        if args.chrom and not id2chrom:
            p.error("--chrom in --toga-dir mode needs --annotation with a chrom column "
                    "(loss_summ_data.tsv has no chromosome).")
        long_rows = []
        for rid, lineage, status, asm in walk_toga_dir(args.toga_dir, args.toga_label, args.toga_level):
            nid = norm(rid)
            gene = id2sym.get(nid, nid)          # translate id->symbol if annotation given
            gene = canon(gene, amap)
            ch = id2chrom.get(nid, "")
            if chrom_keep is not None and ch not in chrom_keep:
                continue
            store.setdefault(gene, {}).setdefault(lineage, set()).add(status)
            if ch:
                chrom_of[gene] = ch
            long_rows.append({"gene": gene, "chrom": ch, "lineage": lineage,
                              "status": status, "assembly": asm, "raw_id": rid})
        toga_long_df = pd.DataFrame(
            long_rows,
            columns=["gene", "chrom", "lineage", "status", "assembly", "raw_id"])
        toga_long_path = write_df(toga_long_df, out_base, "toga_long", as_parquet=True)
        status_aware = True
    else:
        if args.gene_col is None:
            p.error("--gene-col is required for --table (run --list-columns first).")
        if args.chrom and not args.chrom_col:
            p.error("--chrom requires --chrom-col.")
        for name, df in read_table(args.table, args.sheet).items():
            if args.gene_col not in df.columns:
                print(f"  [skip] sheet {name!r}: no {args.gene_col!r}", file=sys.stderr)
                continue
            accumulate_table(df, args.gene_col, args.lineage_col, args.status_col,
                             args.chrom_col, chrom_keep, amap, store, chrom_of,
                             meta_specs, meta_of)
        status_aware = bool(args.status_col)

    hits = {g: store[g] for g in query if g in store}
    lineages = sorted({l for d in hits.values() for l in d})

    # --toga-dir emits pd-lfs multi-file parquet (git-friendly); --table emits TSV
    as_parquet = bool(args.toga_dir)

    # ---- intersections ----
    rows = []
    for g in sorted(hits):
        meta = meta_of.get(g, {})
        for lin in sorted(hits[g]):
            stset = hits[g][lin]
            rows.append({"query_gene": g, "chrom": chrom_of.get(g, ""),
                         **{m: meta.get(m, "") for m in meta_names}, "lineage": lin,
                         "status": most_severe(stset) if status_aware else "HIT",
                         "all_statuses": ";".join(sorted(s for s in stset if s))})
    cols = ["query_gene", "chrom"] + meta_names + ["lineage", "status", "all_statuses"]
    if not status_aware:
        cols = ["query_gene", "chrom"] + meta_names + ["lineage"]
        rows = [{k: r[k] for k in cols} for r in rows]
    written = []
    if args.toga_dir:
        written.append(toga_long_path)
    written.append(write_df(pd.DataFrame(rows, columns=cols),
                            out_base, "intersections", as_parquet))

    if status_aware:
        smat, lmat = [], []
        for g in sorted(hits):
            base = {"gene": g, "chrom": chrom_of.get(g, ""),
                    **{m: meta_of.get(g, {}).get(m, "") for m in meta_names}}
            sr, lr = dict(base), dict(base)
            for lin in lineages:
                sev = most_severe(hits[g][lin]) if lin in hits[g] else ""
                sr[lin] = sev
                lr[lin] = int(sev in loss_set)
            smat.append(sr); lmat.append(lr)
        written.append(write_df(pd.DataFrame(smat, columns=["gene", "chrom"] + meta_names + lineages),
                                out_base, "status_matrix", as_parquet))
        written.append(write_df(pd.DataFrame(lmat, columns=["gene", "chrom"] + meta_names + lineages),
                                out_base, "loss_matrix", as_parquet))
    else:
        pmat = []
        for g in sorted(hits):
            r = {"gene": g, "chrom": chrom_of.get(g, ""),
                 **{m: meta_of.get(g, {}).get(m, "") for m in meta_names}}
            for lin in lineages:
                r[lin] = int(lin in hits[g])
            pmat.append(r)
        written.append(write_df(pd.DataFrame(pmat, columns=["gene", "chrom"] + meta_names + lineages),
                                out_base, "matrix", as_parquet))

    # ---- summary ----
    print(f"query genes:            {len(query)}")
    if chrom_keep:
        print(f"chrom filter:           {','.join(sorted(chrom_keep))}")
    print(f"lineages found:         {len(lineages)}")
    print(f"genes in data:          {len(store)}")
    print(f"intersection (genes):   {len(hits)}")
    if lineages and lineages != ["ALL"]:
        print("per-lineage hits" + (" / losses:" if status_aware else ":"))
        for lin in lineages:
            n = sum(1 for g in hits if lin in hits[g])
            if status_aware:
                nloss = sum(1 for g in hits if lin in hits[g]
                            and most_severe(hits[g][lin]) in loss_set)
                print(f"    {lin:<28} hits={n:<4} loss={nloss}")
            else:
                print(f"    {lin:<28} {n}")
    print("\nwrote: " + "\n       ".join(written))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# ===========================================================================
# DATA / DOWNLOAD NOTES
# ===========================================================================
# TOGA precomputed primate data (Senckenberg mirror, hg38-referenced):
#   wget -e robots=off -r -np -nH --cut-dirs=4 \
#     -A 'loss_summ_data.tsv.gz' \
#     https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates/
#   (the -A filter avoids pulling the multi-GB codon/protein FASTAs; loss_summ_data.tsv.gz
#    is present in EVERY assembly dir under both file-naming conventions in the tree.)
#
# Per-assembly dirs are named  Genus_species__Common_name__AssemblyID  and there can be
# several assemblies per species. The chromosome is NOT in loss_summ_data; build an hg38
# id->symbol->chrom map once (--annotation) from the human reference annotation that TOGA
# used (TOGAInput/human_hg38 in github.com/hillerlab/TOGA2) or any hg38 gene table, e.g. a
# TSV with columns: id  gene  chrom.
#
# Paper supplements (curated, single workbook -> use --table):
#   Mao 2024 Cell  : Europe PMC  PMC10947866/supplementaryFiles  (Data S2)
#   Yoo 2025 Nature: Europe PMC  PMC12058530/supplementaryFiles  (Suppl. Table VIII.34)
#
# Examples
# --------
#   # genome-wide walk, flat output prefix (the '/' keeps it at results/ root)
#   python intersect_primate_lof.py --query data/all_toga_genes.txt \
#       --toga-dir genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates \
#       --annotation data/hg38_toga_genes.tsv --annotation-cols id,gene,chrom \
#       --out results/toga_primates_species
#
#   # a curated subset: bare --out label -> nested under results/<query-stem>/
#   #   writes results/recombination_genes/toga.{intersections,loss_matrix,...}
#   python intersect_primate_lof.py --query recombination_genes.txt \
#       --toga-dir .../Primates --annotation data/hg38_toga_genes.tsv \
#       --chrom chrX --out toga
#
#   # re-cut that subset from the genome-wide parquet without re-walking
#   python intersect_primate_lof.py --query recombination_genes.txt \
#       --table results/toga_primates_species.toga_long.parquet \
#       --gene-col gene --lineage-col lineage --status-col status --chrom-col chrom \
#       --loss-status L --out toga_strictL        # -> results/recombination_genes/
#
#   # single supplementary workbook
#   python intersect_primate_lof.py --query recombination_genes.txt \
#       --table mao2024_suppl.xlsx --list-columns
#
#   # Mao SV-disruption, carrying the stable Ensembl id + biotype through so
#   # symbol-less GENCODE clone names (AC007993) stay joinable / filterable
#   python intersect_primate_lof.py --query data/all_query_genes.txt \
#       --table data/mao_gene_disrupted.tsv \
#       --gene-col GENE --lineage-col LINEAGE --chrom-col CHR \
#       --meta-cols GENE_ID:gene_id,GENE_TYPE:gene_type \
#       --out results/mao2024/gene_disrupted
# ===========================================================================
