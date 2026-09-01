#!/usr/bin/env python3
"""
annotate_sv_lof.py
==================

Annotate a population structural-variant (SV) VCF for HUMAN gene loss-of-function
and emit an intermediate gene x human-lineage table that `intersect_primate_lof.py
--table` then standardises into the project's matrix format.

Built for the 1000 Genomes long-read SV catalog (Schloissnig et al. 2025,
"1KG_ONT_VIENNA"; hg38; ~167k sequence-resolved SVs over 1,019 samples / 26
populations) but works on any hg38 SV VCF with SVTYPE + END/SVLEN.

Why a separate script (not part of the engine)
----------------------------------------------
The engine is pure pandas / no-network. SV->gene annotation needs a VCF parser
(cyvcf2) and an exon model (GTF). Keeping it here mirrors how Mao/Yoo were
pre-cleaned into data/*.tsv before being fed to the engine via --table.

Pipeline
--------
1. CDS interval model (build once, cached): one streaming pass over the rel-110
   Ensembl GTF -> per-coding-gene CDS exon intervals (+ a `canonical` flag for
   MANE_Select / Ensembl_canonical transcripts) and a per-gene CDS span/length
   table. Restricted to gene_ids present in the id,gene,chrom join table so gene
   symbols are byte-identical to every other source.
2. Overlap (pure python, 100 kb binning; no bedtools/pyranges): each SV -> the set
   of coding genes whose CDS it hits.
3. LoF rules per SVTYPE -> status HC_LOF / LC_LOF (see STATUS rules below).
4. Human lineage tiers + superpopulations -> one row per (gene, lineage), carrying
   max AF and carrier counts.

STATUS rules (Decision: DEL/INV strict HC, DUP/INS broad LC)
------------------------------------------------------------
  HC_LOF : DEL overlapping >=1 canonical CDS exon (incl. whole-gene-CDS deletion);
           INV with a breakpoint inside a canonical CDS exon.
  LC_LOF : DUP overlapping CDS; INS inside a CDS exon; INV that spans (engulfs)
           canonical CDS exons without a breakpoint inside them.
  (dropped: overlaps that touch only non-canonical CDS, or only intron/UTR.)
The default loss set downstream is HC_LOF (pass --loss-status HC_LOF to the engine;
--loss-status HC_LOF,LC_LOF for the broad cut).

Human lineage axis (Decision: 3 tiers + 5 superpopulations)
-----------------------------------------------------------
  Human_anyLoF          a LoF SV is present at all (AC>0)
  Human_commonLoF       global AF >= --min-common (default 0.01)
  Human_homozygousLoF   >=1 sample homozygous-alt for the LoF SV
  AFR_LoF AMR_LoF EAS_LoF EUR_LoF SAS_LoF   LoF allele observed in that superpop
Labels deliberately avoid `HSA` (Yoo gains) and `HOMO_SAPIENS` (TOGA reference).

Output TSV columns (engine --table input)
-----------------------------------------
  GENE GENE_ID CHR LINEAGE STATUS SVTYPE CONSEQUENCE AF AC AN N_HOM_ALT N_HET SV_ID
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

BIN = 100_000                       # binning window for the overlap sweep
SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
TAG_RE = re.compile(r'tag "([^"]+)"')


# --------------------------------------------------------------------------- #
# CDS interval model
# --------------------------------------------------------------------------- #
def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_join(join_tsv: str):
    """data/hg38_toga_genes.tsv (id,gene,chrom) -> {ensg: (symbol, chrom)}."""
    id2 = {}
    with _open(join_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        gi, ge, ch = idx["id"], idx["gene"], idx["chrom"]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(gi, ge, ch):
                continue
            ensg = f[gi].strip()
            if ensg:
                # symbol-less genes (91 in the TOGA set) fall back to the ENSG id,
                # mirroring how TOGA outputs keep bare ENSG — never a blank symbol.
                id2[ensg] = (f[ge].strip() or ensg, norm_chrom(f[ch]))
    return id2


def norm_chrom(c) -> str:
    s = str(c).strip()
    if s[:3].lower() == "chr":
        s = s[3:]
    return s.upper()


def build_cds_cache(gtf: str, join_tsv: str, out_intervals: str, out_span: str):
    """Stream the GTF once; write a deduped CDS-interval cache + per-gene span table.

    An interval is `canonical` if any MANE_Select / Ensembl_canonical transcript
    contributes that exact CDS exon. Genes are restricted to the join table and
    named/located from it (authoritative, matches every other source).
    """
    id2 = load_join(join_tsv)
    # (gene_id, start, end) -> canonical(bool); chrom kept per gene_id
    exons: dict[tuple, bool] = {}
    gchrom: dict[str, str] = {}
    n_cds = 0
    with _open(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "CDS":
                continue
            m = GENE_ID_RE.search(f[8])
            if not m:
                continue
            gid = m.group(1)
            if gid not in id2:
                continue
            n_cds += 1
            chrom = norm_chrom(f[0])
            start, end = int(f[3]), int(f[4])
            tags = set(TAG_RE.findall(f[8]))
            canon = ("MANE_Select" in tags) or ("Ensembl_canonical" in tags)
            key = (gid, start, end)
            exons[key] = exons.get(key, False) or canon
            gchrom.setdefault(gid, chrom)

    # per-gene span / canonical length
    span = {}  # gid -> [cds_min, cds_max, canon_len, n_canon, symbol, chrom]
    for (gid, start, end), canon in exons.items():
        sym, jchrom = id2[gid]
        chrom = gchrom.get(gid, jchrom)
        s = span.get(gid)
        if s is None:
            span[gid] = [start, end, 0, 0, sym, chrom]
            s = span[gid]
        else:
            if start < s[0]:
                s[0] = start
            if end > s[1]:
                s[1] = end
        if canon:
            s[2] += (end - start + 1)
            s[3] += 1

    # write interval cache, sorted by chrom,start
    rows = []
    for (gid, start, end), canon in exons.items():
        sym, jchrom = id2[gid]
        chrom = gchrom.get(gid, jchrom)
        rows.append((chrom, start, end, gid, sym, 1 if canon else 0))
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    op = gzip.open(out_intervals, "wt") if str(out_intervals).endswith(".gz") else open(out_intervals, "w")
    with op as out:
        out.write("chrom\tstart\tend\tgene_id\tgene\tcanonical\n")
        for r in rows:
            out.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\n")
    with open(out_span, "w") as out:
        out.write("gene_id\tgene\tchrom\tcds_min\tcds_max\tcds_len_canonical\tn_cds_exons_canonical\n")
        for gid, s in sorted(span.items(), key=lambda kv: (kv[1][5], kv[1][0])):
            out.write(f"{gid}\t{s[4]}\t{s[5]}\t{s[0]}\t{s[1]}\t{s[2]}\t{s[3]}\n")
    print(f"[cds-cache] CDS lines kept: {n_cds:,} | distinct exons: {len(exons):,} | "
          f"coding genes: {len(span):,}", file=sys.stderr)
    print(f"[cds-cache] wrote {out_intervals} and {out_span}", file=sys.stderr)


def load_cds_cache(intervals: str, span: str):
    """Load the cache into a binned interval index + per-gene span dict."""
    bins: dict[str, dict[int, list]] = defaultdict(lambda: defaultdict(list))
    with _open(intervals) as fh:
        fh.readline()
        for line in fh:
            chrom, start, end, gid, gene, canon = line.rstrip("\n").split("\t")
            start, end, canon = int(start), int(end), canon == "1"
            iv = (start, end, gid, canon)
            for b in range(start // BIN, end // BIN + 1):
                bins[chrom][b].append(iv)
    spand = {}
    with _open(span) as fh:
        fh.readline()
        for line in fh:
            gid, gene, chrom, cmin, cmax, clen, nex = line.rstrip("\n").split("\t")
            spand[gid] = (int(cmin), int(cmax), int(clen), gene, chrom)
    return bins, spand


# --------------------------------------------------------------------------- #
# Overlap + LoF rules
# --------------------------------------------------------------------------- #
def overlapping_genes(chrom, s, e, bins):
    """Return {gene_id: dict(canon, noncanon, n_canon, canon_bp, span_intervals)}.

    s,e are 1-based inclusive SV bounds. Aggregates per gene over its CDS exons.
    """
    cb = bins.get(chrom)
    if not cb:
        return {}
    seen = set()
    per = {}
    for b in range(s // BIN, e // BIN + 1):
        for iv in cb.get(b, ()):  # iv = (start, end, gid, canon)
            if iv in seen:
                continue
            seen.add(iv)
            ivs, ive, gid, canon = iv
            if ivs > e or ive < s:           # inclusive overlap test
                continue
            g = per.get(gid)
            if g is None:
                g = per[gid] = {"canon": False, "noncanon": False, "n_canon": 0,
                                "canon_del_len": 0, "canon_bp": False, "canon_pos": False}
            if canon:
                g["canon"] = True
                g["n_canon"] += 1
                g["canon_del_len"] += min(e, ive) - max(s, ivs) + 1
                # breakpoint of SV falls inside this canonical exon?
                if ivs <= s <= ive or ivs <= e <= ive:
                    g["canon_bp"] = True
                # SV start position inside this canonical exon (for INS)?
                if ivs <= s <= ive:
                    g["canon_pos"] = True
            else:
                g["noncanon"] = True
    return per


def classify(svtype, s, e, g, gspan, net_del=0):
    """Return (status, consequence) for one (gene, SV), or (None, None) to drop.

    g = per-gene overlap dict from overlapping_genes(); gspan = span tuple;
    net_del = bp of reference lost (len(REF)-len(ALT)), for COMPLEX events.
    """
    svtype = (svtype or "").upper()
    if not g["canon"]:
        return None, None                     # only non-canonical / no CDS -> drop
    cds_min, cds_max = gspan[0], gspan[1]
    canon_len = gspan[2] or 1
    if svtype == "DEL":
        if s <= cds_min and e >= cds_max:
            return "HC_LOF", "whole_gene_CDS_deletion"
        frac = g["canon_del_len"] / canon_len
        if frac >= 0.10:
            return "HC_LOF", f"CDS_deletion_{frac:.0%}_of_canonical"
        return "HC_LOF", "canonical_CDS_exon_deletion"
    if svtype == "INV":
        if g["canon_bp"]:
            return "HC_LOF", "inversion_breakpoint_in_CDS"
        return "LC_LOF", "inversion_spanning_CDS"
    if svtype == "INS":
        if g["canon_pos"]:
            return "LC_LOF", "insertion_in_CDS_exon"
        return None, None                     # INS not inside a canonical exon -> drop
    if svtype == "COMPLEX":
        if net_del >= 50:                     # net loss of coding sequence
            return "HC_LOF", "complex_deletion_in_CDS"
        return "LC_LOF", "complex_over_CDS"
    if svtype in ("DUP", "CNV"):
        return "LC_LOF", f"{svtype.lower()}_over_CDS"
    if svtype in ("BND", "TRA"):
        if g["canon_bp"]:
            return "HC_LOF", "breakend_in_CDS"
        return None, None
    # unknown SVTYPE: conservative LC if it touches canonical CDS
    return "LC_LOF", f"{svtype.lower() or 'sv'}_over_CDS"


# --------------------------------------------------------------------------- #
# VCF parsing helpers
# --------------------------------------------------------------------------- #
SVID_TYPE_RE = re.compile(r'^[A-Za-z]+')


def sv_typespan(v, s, min_sv_len):
    """Return (svtype, start, end, ref_len, alt_len) or None to skip.

    Handles BOTH symbolic SVs (INFO SVTYPE + END/SVLEN) and sequence-resolved
    biallelic records (REF/ALT are sequences; type from the ID token or REF/ALT
    lengths; reference span = len(REF)). Sub-SV-size records are skipped.
    """
    alt = (v.ALT[0] if v.ALT else "") or ""
    ref_len = len(v.REF or "")
    svt = v.INFO.get("SVTYPE")
    if svt:                                       # symbolic SV
        end = v.INFO.get("END")
        if end is not None:
            end = int(end)
        else:
            svlen = v.INFO.get("SVLEN")
            if svlen is not None:
                svlen = svlen[0] if isinstance(svlen, (list, tuple)) else svlen
                end = s + abs(int(svlen))
            else:
                end = s + max(ref_len, 1) - 1
        alt_len = 0 if alt.startswith("<") else len(alt)
        return str(svt).upper(), s, end, ref_len, alt_len
    # sequence-resolved record
    alt_len = len(alt)
    if max(ref_len, alt_len) < min_sv_len:        # indel/SNV below SV size -> skip
        return None
    end = s + max(ref_len, 1) - 1
    t = None
    if v.ID:
        parts = v.ID.split("-")
        if len(parts) >= 3:
            m = SVID_TYPE_RE.match(parts[2])
            if m:
                t = m.group(0).upper()
    if t not in ("DEL", "INS", "DUP", "INV", "COMPLEX"):
        t = "DEL" if ref_len > alt_len else ("INS" if alt_len > ref_len else "COMPLEX")
    return t, s, end, ref_len, alt_len


def info_float(v, field):
    x = v.INFO.get(field)
    if x is None:
        return None
    if isinstance(x, (list, tuple)):
        x = x[0]
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def load_superpop_map(path):
    """sample<TAB>superpop panel -> {sample: SUPERPOP}. Tolerates IGSR panel
    headers (picks a column named superpop/super_pop/pop-style if present)."""
    smap = {}
    with _open(path) as fh:
        first = fh.readline().rstrip("\n").split("\t")
        cols = [c.strip().lower() for c in first]
        si = sp = None
        if "sample" in cols:
            si = cols.index("sample")
        for cand in ("superpop", "super_pop", "super_population", "population_group"):
            if cand in cols:
                sp = cols.index(cand)
                break
        if si is not None and sp is not None:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > max(si, sp):
                    smap[f[si].strip()] = f[sp].strip().upper()
        else:  # headerless: assume col0=sample, col-1=superpop
            for line in [first] + [l.rstrip("\n").split("\t") for l in fh]:
                f = line if isinstance(line, list) else line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    smap[f[0].strip()] = f[-1].strip().upper()
    return smap


# --------------------------------------------------------------------------- #
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtf", default="tmp/Homo_sapiens.GRCh38.110.gtf.gz",
                   help="Ensembl GTF (same release as --join) for the CDS model.")
    p.add_argument("--join", default="data/hg38_toga_genes.tsv",
                   help="id,gene,chrom table; restricts + names the gene universe.")
    p.add_argument("--cds-cache", default="data/hg38_cds_intervals.tsv.gz",
                   help="CDS interval cache (built on first run).")
    p.add_argument("--span-cache", default="data/hg38_gene_cds_span.tsv",
                   help="Per-gene CDS span/length cache.")
    p.add_argument("--build-cds-cache", action="store_true",
                   help="(Re)build the CDS caches from --gtf/--join, then exit.")
    p.add_argument("--vcf", default=None, help="SV VCF(.gz) to annotate.")
    p.add_argument("--out", default="data/human_1kgp_sv_lof.tsv",
                   help="Output intermediate human-source TSV.")
    p.add_argument("--af-field", default="AF")
    p.add_argument("--ac-field", default="AC")
    p.add_argument("--an-field", default="AN")
    p.add_argument("--pop-af-prefix", default="AF_",
                   help="INFO per-population AF prefix (e.g. AF_AFR). Used if present.")
    p.add_argument("--min-common", type=float, default=0.01,
                   help="Global AF threshold for Human_commonLoF (default 0.01).")
    p.add_argument("--min-af", type=float, default=0.0,
                   help="Skip SVs below this global AF before annotation.")
    p.add_argument("--min-sv-len", type=int, default=50,
                   help="Min REF/ALT length for a sequence-resolved record to count "
                        "as an SV (default 50; drops SNVs/small indels). Ignored for "
                        "symbolic SVs that carry SVTYPE.")
    p.add_argument("--superpop-map", default=None,
                   help="sample->superpop panel; enables per-superpop AF from genotypes "
                        "when the VCF lacks AF_<POP> INFO fields.")
    p.add_argument("--chrom", default=None, help="Restrict to one chromosome (e.g. 21).")
    args = p.parse_args()

    need_build = args.build_cds_cache or not (Path(args.cds_cache).exists()
                                              and Path(args.span_cache).exists())
    if need_build:
        build_cds_cache(args.gtf, args.join, args.cds_cache, args.span_cache)
        if args.build_cds_cache:
            return 0

    if not args.vcf:
        p.error("--vcf is required (unless only --build-cds-cache).")

    from cyvcf2 import VCF
    bins, span = load_cds_cache(args.cds_cache, args.span_cache)

    vcf = VCF(args.vcf)
    samples = list(vcf.samples)
    smap = load_superpop_map(args.superpop_map) if args.superpop_map else {}
    # sample-index lists per superpop (for genotype-based per-pop AF)
    pop_idx = {pop: [i for i, s in enumerate(samples) if smap.get(s) == pop]
               for pop in SUPERPOPS} if smap else {}

    chrom_keep = norm_chrom(args.chrom) if args.chrom else None

    # aggregate per (gene, lineage): best status, max AF, max hom, representative SV
    agg: dict[tuple, dict] = {}
    n_sv = n_hit = 0
    for v in vcf:
        chrom = norm_chrom(v.CHROM)
        if chrom_keep and chrom != chrom_keep:
            continue
        n_sv += 1
        s = v.start + 1                       # cyvcf2 .start is 0-based
        ts = sv_typespan(v, s, args.min_sv_len)
        if ts is None:                        # sub-SV-size record (SNV/small indel)
            continue
        svtype, s, e, ref_len, alt_len = ts
        net_del = max(0, ref_len - alt_len)
        if e < s:
            s, e = e, s
        af = info_float(v, args.af_field)
        ac = v.INFO.get(args.ac_field)
        an = v.INFO.get(args.an_field)
        ac = ac[0] if isinstance(ac, (list, tuple)) else ac
        # homozygous/het carriers: prefer INFO AC_Hom/AC_Het (allele counts in
        # hom/het genotypes), else fall back to decoded genotypes.
        ac_hom = v.INFO.get("AC_Hom")
        ac_het = v.INFO.get("AC_Het")
        if ac_hom is not None:
            ac_hom = ac_hom[0] if isinstance(ac_hom, (list, tuple)) else ac_hom
            n_hom_alt = int(ac_hom) // 2      # alt-allele count in homs -> #samples
        else:
            try:
                n_hom_alt = int(v.num_hom_alt)
            except Exception:
                n_hom_alt = None
        if ac_het is not None:
            ac_het = ac_het[0] if isinstance(ac_het, (list, tuple)) else ac_het
            n_het = int(ac_het)
        else:
            try:
                n_het = int(v.num_het)
            except Exception:
                n_het = None
        if af is None and ac and an:
            try:
                af = float(ac) / float(an)
            except (TypeError, ZeroDivisionError):
                af = None
        if af is not None:
            af = round(af, 6)                  # cyvcf2 returns float32; tidy it
        if af is not None and af < args.min_af:
            continue

        genes = overlapping_genes(chrom, s, e, bins)
        if not genes:
            continue

        # per-superpop AF (INFO field first, else genotype counts)
        pop_af = {}
        for pop in SUPERPOPS:
            x = info_float(v, f"{args.pop_af_prefix}{pop}")
            if x is not None:
                pop_af[pop] = x
        if not pop_af and pop_idx:
            gts = v.gt_types  # cyvcf2: 0 HOM_REF, 1 HET, 3 HOM_ALT, 2 UNKNOWN
            for pop, idxs in pop_idx.items():
                if not idxs:
                    continue
                ac_p = an_p = 0
                for i in idxs:
                    t = gts[i]
                    if t == 1:
                        ac_p += 1; an_p += 2
                    elif t == 3:
                        ac_p += 2; an_p += 2
                    elif t == 0:
                        an_p += 2
                pop_af[pop] = (ac_p / an_p) if an_p else 0.0

        hit_any = False
        for gid, g in genes.items():
            status, conseq = classify(svtype, s, e, g, span[gid], net_del)
            if status is None:
                continue
            hit_any = True
            sym, gchrom = span[gid][3], span[gid][4]
            # which human lineages does this SV populate?
            lineages = ["Human_anyLoF"]
            if af is not None and af >= args.min_common:
                lineages.append("Human_commonLoF")
            if n_hom_alt:
                lineages.append("Human_homozygousLoF")
            for pop in SUPERPOPS:
                if pop_af.get(pop, 0.0) > 0.0:
                    lineages.append(f"{pop}_LoF")
            rank = {"LC_LOF": 0, "HC_LOF": 1}
            for lin in lineages:
                key = (gid, lin)
                cur = agg.get(key)
                cand = {"gene": sym, "gene_id": gid, "chrom": gchrom, "status": status,
                        "svtype": svtype or "", "consequence": conseq,
                        "af": af if af is not None else "",
                        "ac": ac if ac is not None else "",
                        "an": an if an is not None else "",
                        "n_hom_alt": n_hom_alt if n_hom_alt is not None else "",
                        "n_het": n_het if n_het is not None else "",
                        "sv_id": v.ID or f"{chrom}-{s}-{svtype}-{e-s}"}
                if cur is None:
                    agg[key] = cand
                    continue
                # keep the single strongest SV (status tier, then AF) as the
                # coherent representative for this (gene, lineage) cell, so the
                # carried svtype/consequence/af all describe one real SV.
                cur_af = cur["af"] if isinstance(cur["af"], float) else -1.0
                cand_af = af if isinstance(af, float) else -1.0
                if (rank[status], cand_af) > (rank[cur["status"]], cur_af):
                    agg[key] = cand
        if hit_any:
            n_hit += 1

    # write output
    cols = ["GENE", "GENE_ID", "CHR", "LINEAGE", "STATUS", "SVTYPE", "CONSEQUENCE",
            "AF", "AC", "AN", "N_HOM_ALT", "N_HET", "SV_ID"]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    # Order rows so each gene's STRONGEST SV row comes first: the engine carries
    # per-gene meta as "first non-empty wins", so this makes the gene-level
    # svtype/consequence/af describe that gene's most-severe SV.
    srank = {"LC_LOF": 0, "HC_LOF": 1}
    def _key(kv):
        (gid, lin), r = kv
        af = r["af"] if isinstance(r["af"], float) else -1.0
        return (r["gene"], -srank[r["status"]], -af, lin)
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for (gid, lin), r in sorted(agg.items(), key=_key):
            out.write("\t".join(str(x) for x in [
                r["gene"], r["gene_id"], r["chrom"], lin, r["status"], r["svtype"],
                r["consequence"], r["af"], r["ac"], r["an"], r["n_hom_alt"],
                r["n_het"], r["sv_id"]]) + "\n")

    n_genes = len({k[0] for k in agg})
    print(f"[annotate] SVs scanned: {n_sv:,} | SVs hitting a CDS gene: {n_hit:,}", file=sys.stderr)
    print(f"[annotate] (gene,lineage) rows: {len(agg):,} | distinct genes: {n_genes:,}", file=sys.stderr)
    print(f"[annotate] wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
