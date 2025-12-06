#!/usr/bin/env python3
"""
Compile amino acid variants with >10% population allele frequency
in human, chimpanzee, and gorilla genes.

Data sources:
- Human: Ensembl Variation (gnomAD, 1000 Genomes), or direct gnomAD download
- Chimpanzee: Ensembl Variation (Great Ape Genome Project data)
- Gorilla: Ensembl Variation (Great Ape Genome Project data)

Note: Population variant data for great apes is much more limited than for humans.

Usage:
    python primate_aa_variants.py --method gnomad    # Best for human
    python primate_aa_variants.py --method ensembl   # Works for all species
    python primate_aa_variants.py --method biomart   # Bulk download

Author: Generated for comparative genomics analysis
"""

import requests
import pandas as pd
import time
from typing import Optional
from dataclasses import dataclass, asdict
from collections import defaultdict
import json
from pathlib import Path
from io import StringIO
import pickle
import sys

# Ensembl REST API
ENSEMBL_REST = "https://rest.ensembl.org"

# Species configuration
SPECIES_CONFIG = {
    "human": {
        "ensembl_name": "homo_sapiens",
        "taxon_id": "9606",
        "assembly": "GRCh38",
    },
    "chimpanzee": {
        "ensembl_name": "pan_troglodytes",
        "taxon_id": "9598",
        "assembly": "Pan_tro_3.0",
    },
    "gorilla": {
        "ensembl_name": "gorilla_gorilla",
        "taxon_id": "9593",
        "assembly": "gorGor4",
    }
}

MAF_THRESHOLD = 0.10


@dataclass
class AAVariant:
    """Amino acid variant with population frequency."""
    species: str
    gene_id: str
    gene_name: str
    transcript_id: str
    variant_id: str
    chromosome: str
    genomic_position: int
    ref_allele: str
    alt_allele: str
    aa_position: int
    ref_aa: str
    alt_aa: str
    allele_frequency: float
    consequence: str
    source: str


class EnsemblClient:
    """Client for Ensembl REST API with rate limiting."""

    def __init__(self, verbose: bool = True):
        self.session = requests.Session()
        self.session.headers.update({
            "Content-Type": "application/json",
            "Accept": "application/json"
        })
        self.verbose = verbose

    def get(self, endpoint: str, params: dict = None, allow_404: bool = False) -> dict:
        """GET request with retry logic.

        Args:
            endpoint: API endpoint
            params: Query parameters
            allow_404: If True, return {} for 404 errors. If False, raise exception.

        Raises:
            requests.exceptions.HTTPError: For non-404 HTTP errors
            requests.exceptions.RequestException: For network errors after retries
        """
        url = f"{ENSEMBL_REST}{endpoint}"

        for attempt in range(3):
            try:
                response = self.session.get(url, params=params, timeout=60)

                if response.status_code == 429:
                    retry_after = int(response.headers.get("Retry-After", 5))
                    if self.verbose:
                        print(f"\n  Rate limited, waiting {retry_after}s...", flush=True)
                    time.sleep(retry_after)
                    continue

                if response.status_code == 404:
                    if allow_404:
                        return {}
                    else:
                        raise requests.exceptions.HTTPError(
                            f"404 Not Found: {url}", response=response
                        )

                if response.status_code == 400:
                    error_msg = response.text
                    raise requests.exceptions.HTTPError(
                        f"400 Bad Request: {url}\n{error_msg}", response=response
                    )

                response.raise_for_status()
                return response.json()

            except requests.exceptions.RequestException as e:
                if attempt == 2:
                    raise Exception(f"Failed after 3 attempts: {url}\n{str(e)}")
                if self.verbose:
                    print(f"\n  Retry {attempt + 1}/3: {str(e)[:100]}...", flush=True)
                time.sleep(2 ** attempt)

        raise Exception(f"Failed to fetch: {url}")

    def post(self, endpoint: str, data: dict) -> dict:
        """POST request with retry logic."""
        url = f"{ENSEMBL_REST}{endpoint}"

        for attempt in range(3):
            try:
                response = self.session.post(url, json=data, timeout=120)

                if response.status_code == 429:
                    retry_after = int(response.headers.get("Retry-After", 5))
                    if self.verbose:
                        print(f"\n  Rate limited, waiting {retry_after}s...", flush=True)
                    time.sleep(retry_after)
                    continue

                response.raise_for_status()
                return response.json()

            except requests.exceptions.RequestException as e:
                if attempt == 2:
                    raise Exception(f"Failed after 3 attempts: {url}\n{str(e)}")
                if self.verbose:
                    print(f"\n  Retry {attempt + 1}/3: {str(e)[:100]}...", flush=True)
                time.sleep(2 ** attempt)

        raise Exception(f"Failed to post: {url}")


# ============================================================================
# Method 1: Ensembl REST API (works for all species, but slow)
# ============================================================================

def get_chromosome_info(client: EnsemblClient, species: str) -> dict:
    """Get chromosome names and lengths for a species.

    Returns:
        dict: {chrom_name: length}
    """
    ensembl_species = SPECIES_CONFIG[species]["ensembl_name"]
    assembly = SPECIES_CONFIG[species]["assembly"]

    print(f"Fetching chromosome info for {species} ({assembly})...")

    # Get assembly info
    result = client.get(f"/info/assembly/{ensembl_species}")

    if not result or "top_level_region" not in result:
        raise Exception(f"Failed to get assembly info for {species}")

    # Extract chromosome lengths (only main chromosomes, not patches/scaffolds)
    chrom_lengths = {}
    for region in result["top_level_region"]:
        name = region["name"]
        length = region["length"]

        # For humans: chromosomes 1-22, X, Y
        # For great apes: similar, but may vary
        if species == "human":
            if name in [str(i) for i in range(1, 23)] + ["X", "Y"]:
                chrom_lengths[name] = length
        else:
            # For other species, take chromosomes (not scaffolds/contigs)
            if region.get("coord_system_level") == "chromosome":
                chrom_lengths[name] = length

    print(f"  Found {len(chrom_lengths)} chromosomes")
    return chrom_lengths


def get_protein_coding_genes(
    client: EnsemblClient,
    species: str,
    cache_file: Optional[str] = None
) -> list[dict]:
    """Fetch protein-coding genes for a species.

    Args:
        client: Ensembl API client
        species: Species name (human, chimpanzee, gorilla)
        cache_file: Optional cache file to save/load gene list

    Returns:
        List of gene dictionaries
    """
    # Try loading from cache
    if cache_file and Path(cache_file).exists():
        print(f"Loading genes from cache: {cache_file}")
        with open(cache_file, "rb") as f:
            genes = pickle.load(f)
        print(f"  Loaded {len(genes)} genes from cache")
        return genes

    ensembl_species = SPECIES_CONFIG[species]["ensembl_name"]
    print(f"Fetching protein-coding genes for {species}...")

    # Get chromosome lengths
    chrom_lengths = get_chromosome_info(client, species)

    genes = []
    seen_gene_ids = set()

    # Ensembl API has 5Mb max region size
    CHUNK_SIZE = 5000000

    for chrom, chrom_length in chrom_lengths.items():
        print(f"  Chr {chrom} ({chrom_length:,} bp)...", end=" ", flush=True)

        chrom_genes = []
        chunks_queried = 0
        chunks_failed = 0

        # Query entire chromosome in chunks
        for start in range(1, chrom_length + 1, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE - 1, chrom_length)

            endpoint = f"/overlap/region/{ensembl_species}/{chrom}:{start}-{end}"
            params = {"feature": "gene", "biotype": "protein_coding"}

            try:
                result = client.get(endpoint, params, allow_404=True)
                chunks_queried += 1

                if isinstance(result, list):
                    for g in result:
                        if g.get("biotype") == "protein_coding":
                            gene_id = g["gene_id"]
                            # Avoid duplicates from overlapping chunks
                            if gene_id not in seen_gene_ids:
                                seen_gene_ids.add(gene_id)
                                chrom_genes.append({
                                    "gene_id": gene_id,
                                    "gene_name": g.get("external_name", gene_id),
                                    "chromosome": chrom,
                                    "start": g.get("start"),
                                    "end": g.get("end")
                                })

                time.sleep(0.05)  # Rate limiting

            except Exception as e:
                chunks_failed += 1
                print(f"\n  ERROR in chunk {start}-{end}: {e}", file=sys.stderr)
                # Continue with next chunk instead of failing completely
                continue

        genes.extend(chrom_genes)
        print(f"{len(chrom_genes)} genes ({chunks_queried} chunks, {chunks_failed} failed)")

    print(f"Total: {len(genes)} genes across {len(chrom_lengths)} chromosomes")

    # Save to cache
    if cache_file:
        cache_path = Path(cache_file)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_file, "wb") as f:
            pickle.dump(genes, f)
        print(f"Saved gene list to cache: {cache_file}")

    return genes


def get_gene_variants(
    client: EnsemblClient,
    gene: dict,
    species: str,
    maf_threshold: float = MAF_THRESHOLD
) -> list[AAVariant]:
    """Get missense variants with MAF >= threshold for a gene."""
    ensembl_species = SPECIES_CONFIG[species]["ensembl_name"]
    variants = []
    
    # Get gene details
    gene_info = client.get(f"/lookup/id/{gene['gene_id']}", {"expand": "1"}, allow_404=True)
    if not gene_info:
        return variants

    gene_name = gene_info.get("display_name", gene["gene_id"])

    # Find canonical transcript
    canonical_transcript = None
    for t in gene_info.get("Transcript", []):
        if t.get("is_canonical"):
            canonical_transcript = t["id"]
            break

    if not canonical_transcript:
        return variants

    chrom = gene["chromosome"]
    start = gene["start"]
    end = gene["end"]

    # Get variants in gene region
    endpoint = f"/overlap/region/{ensembl_species}/{chrom}:{start}-{end}"
    params = {"feature": "variation"}
    var_result = client.get(endpoint, params, allow_404=True)

    if not isinstance(var_result, list):
        return variants
    
    # Process each variant
    var_ids = [v.get("id") for v in var_result if v.get("id")]
    
    # Batch lookup (up to 200 at a time)
    for i in range(0, len(var_ids), 200):
        batch = var_ids[i:i+200]
        batch_result = client.post(f"/variation/{ensembl_species}", {"ids": batch})
        
        if not batch_result:
            continue
        
        for var_id, var_detail in batch_result.items():
            if not isinstance(var_detail, dict):
                continue
            
            # Get max population frequency
            max_freq = 0.0
            freq_source = ""
            
            for pop in var_detail.get("populations", []):
                freq = pop.get("frequency")
                if freq and freq > max_freq:
                    max_freq = freq
                    freq_source = pop.get("population", "unknown")
            
            # Check MAF field
            maf = var_detail.get("MAF")
            if maf and maf > max_freq:
                max_freq = maf
                freq_source = "MAF"
            
            if max_freq < maf_threshold:
                continue
            
            # Get mappings and consequences
            for mapping in var_detail.get("mappings", []):
                allele_string = mapping.get("allele_string", "")
                if "/" not in allele_string:
                    continue
                
                ref, *alts = allele_string.split("/")
                genomic_pos = mapping.get("start", 0)
            
            # Need VEP to get protein consequences
            vep_result = client.get(f"/vep/{ensembl_species}/id/{var_id}", allow_404=True)

            if not vep_result or not isinstance(vep_result, list):
                continue
            
            for vep_entry in vep_result:
                for tc in vep_entry.get("transcript_consequences", []):
                    if tc.get("transcript_id") != canonical_transcript:
                        continue
                    
                    consequences = tc.get("consequence_terms", [])
                    if "missense_variant" not in consequences:
                        continue
                    
                    aa_change = tc.get("amino_acids", "")
                    if "/" not in aa_change:
                        continue
                    
                    ref_aa, alt_aa = aa_change.split("/")[:2]
                    protein_pos = tc.get("protein_start")
                    
                    if not protein_pos:
                        continue
                    
                    variants.append(AAVariant(
                        species=species,
                        gene_id=gene["gene_id"],
                        gene_name=gene_name,
                        transcript_id=canonical_transcript,
                        variant_id=var_id,
                        chromosome=chrom,
                        genomic_position=genomic_pos,
                        ref_allele=ref,
                        alt_allele=alts[0] if alts else "",
                        aa_position=protein_pos,
                        ref_aa=ref_aa,
                        alt_aa=alt_aa,
                        allele_frequency=max_freq,
                        consequence="missense_variant",
                        source=freq_source
                    ))
        
        time.sleep(0.1)
    
    return variants


def compile_via_ensembl_api(
    species_list: list[str],
    gene_limit: int = None,
    maf_threshold: float = MAF_THRESHOLD,
    output_prefix: str = "aa_variants",
    cache_dir: str = ".cache"
) -> pd.DataFrame:
    """Compile variants using Ensembl REST API."""
    client = EnsemblClient()
    all_variants = []

    for species in species_list:
        print(f"\n{'='*60}")
        print(f"Processing {species.upper()}")
        print(f"{'='*60}")

        cache_file = f"{cache_dir}/{species}_genes.pkl"
        genes = get_protein_coding_genes(client, species, cache_file=cache_file)
        
        if gene_limit:
            genes = genes[:gene_limit]
            print(f"Limited to {gene_limit} genes")
        
        species_variants = []
        
        for i, gene in enumerate(genes):
            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{len(genes)} genes, {len(species_variants)} variants found")
            
            variants = get_gene_variants(client, gene, species, maf_threshold)
            species_variants.extend(variants)
            time.sleep(0.05)
        
        print(f"Found {len(species_variants)} variants for {species}")
        all_variants.extend(species_variants)
        
        # Save intermediate
        if species_variants:
            df = pd.DataFrame([asdict(v) for v in species_variants])
            df.to_csv(f"{output_prefix}_{species}.csv", index=False)
    
    if all_variants:
        combined_df = pd.DataFrame([asdict(v) for v in all_variants])
        combined_df.to_csv(f"{output_prefix}_combined.csv", index=False)
        return combined_df
    
    return pd.DataFrame()


# ============================================================================
# Method 2: BioMart (faster bulk download)
# ============================================================================

def fetch_via_biomart(
    species: str,
    maf_threshold: float = MAF_THRESHOLD,
    output_file: str = None
) -> pd.DataFrame:
    """Fetch variants using Ensembl BioMart."""
    
    # Map species to BioMart dataset names
    dataset_map = {
        "human": "hsapiens_snp",
        "chimpanzee": "ptroglodytes_snp",
        "gorilla": "ggorilla_snp"
    }
    
    dataset = dataset_map.get(species)
    if not dataset:
        print(f"Unknown species: {species}")
        return pd.DataFrame()
    
    biomart_url = "https://www.ensembl.org/biomart/martservice"
    
    # Query for missense variants with frequency data
    xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
    <Dataset name="{dataset}" interface="default">
        <Filter name="consequence_type_tv" value="missense_variant"/>
        <Attribute name="refsnp_id"/>
        <Attribute name="chr_name"/>
        <Attribute name="chrom_start"/>
        <Attribute name="consequence_type_tv"/>
        <Attribute name="ensembl_gene_stable_id"/>
        <Attribute name="associated_gene"/>
        <Attribute name="ensembl_transcript_stable_id"/>
        <Attribute name="minor_allele"/>
        <Attribute name="minor_allele_freq"/>
        <Attribute name="peptide_location"/>
        <Attribute name="amino_acid_variation"/>
    </Dataset>
</Query>"""
    
    print(f"Querying BioMart for {species} missense variants...")
    print("This may take several minutes for genome-wide queries...")
    
    try:
        response = requests.post(
            biomart_url,
            data={"query": xml_query},
            timeout=600
        )
        response.raise_for_status()
        
        df = pd.read_csv(StringIO(response.text), sep="\t")
        
        # Filter by MAF
        if "Minor allele frequency" in df.columns:
            df = df[df["Minor allele frequency"] >= maf_threshold]
        
        print(f"Retrieved {len(df)} variants")
        
        if output_file:
            df.to_csv(output_file, index=False)
            print(f"Saved to {output_file}")
        
        return df
        
    except Exception as e:
        print(f"BioMart query failed: {e}")
        return pd.DataFrame()


# ============================================================================
# Method 3: gnomAD (best for human, comprehensive population data)
# ============================================================================

def parse_hgvsp(hgvsp: str) -> Optional[tuple[str, int, str]]:
    """Parse HGVS protein notation to extract amino acid change.

    Args:
        hgvsp: HGVS protein notation (e.g., "p.Arg175His", "p.G12V")

    Returns:
        Tuple of (ref_aa, position, alt_aa) or None if not a simple substitution

    Examples:
        >>> parse_hgvsp("p.Arg175His")
        ('Arg', 175, 'His')
        >>> parse_hgvsp("p.G12V")
        ('G', 12, 'V')
    """
    if not hgvsp or not hgvsp.startswith("p."):
        return None

    # Remove "p." prefix
    change = hgvsp[2:]

    # Match pattern: amino_acid + position + amino_acid
    # Handles both 3-letter (Arg) and 1-letter (R) codes
    import re
    # Pattern: (letters)(digits)(letters)
    match = re.match(r"^([A-Za-z]+)(\d+)([A-Za-z]+)$", change)

    if match:
        ref_aa = match.group(1)
        position = int(match.group(2))
        alt_aa = match.group(3)
        return (ref_aa, position, alt_aa)

    return None


def fetch_gnomad_via_api(
    genes: list[str],
    maf_threshold: float = MAF_THRESHOLD
) -> pd.DataFrame:
    """Fetch variants from gnomAD GraphQL API (gene by gene).

    Retrieves comprehensive variant annotations including:
    - Population frequencies (global and per-population)
    - In silico predictors (CADD, REVEL, SpliceAI, AlphaMissense, etc.)
    - Functional predictions (PolyPhen, SIFT)
    - Loss-of-function annotations
    """

    gnomad_url = "https://gnomad.broadinstitute.org/api/"

    query = """
    query GeneVariants($geneSymbol: String!, $dataset: DatasetId!) {
        gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
            gene_id
            symbol
            variants(dataset: $dataset) {
                variant_id
                pos
                ref
                alt
                rsids
                flags
                in_silico_predictors {
                    id
                    value
                }
                exome {
                    af
                    ac
                    an
                    ac_hom
                    ac_hemi
                    filters
                    populations {
                        id
                        ac
                        an
                        ac_hom
                    }
                    faf95 {
                        popmax
                        popmax_population
                    }
                    fafmax {
                        faf95_max
                        faf95_max_gen_anc
                        faf99_max
                        faf99_max_gen_anc
                    }
                }
                genome {
                    af
                    ac
                    an
                    ac_hom
                    ac_hemi
                    filters
                    populations {
                        id
                        ac
                        an
                        ac_hom
                    }
                    fafmax {
                        faf95_max
                        faf95_max_gen_anc
                    }
                }
                transcript_consequence {
                    gene_symbol
                    transcript_id
                    consequence_terms
                    major_consequence
                    hgvsp
                    hgvsc
                    is_canonical
                    polyphen_prediction
                    sift_prediction
                    lof
                    lof_filter
                    lof_flags
                }
            }
        }
    }
    """
    
    all_variants = []
    
    for gene_symbol in genes:
        print(f"  Querying gnomAD for {gene_symbol}...", end=" ", flush=True)
        
        try:
            response = requests.post(
                gnomad_url,
                json={
                    "query": query,
                    "variables": {
                        "geneSymbol": gene_symbol,
                        "dataset": "gnomad_r4"
                    }
                },
                timeout=60
            )
            response.raise_for_status()
            data = response.json()
            
            gene_data = data.get("data", {}).get("gene")
            if not gene_data:
                print("not found")
                continue
            
            gene_id = gene_data.get("gene_id")
            variants = gene_data.get("variants", [])
            
            count = 0
            for var in variants:
                # Get exome and genome data
                exome = var.get("exome") or {}
                genome = var.get("genome") or {}

                # Get max AF between exome and genome
                exome_af = exome.get("af") or 0
                genome_af = genome.get("af") or 0
                max_af = max(exome_af, genome_af)

                if max_af < maf_threshold:
                    continue

                # Check for missense in canonical transcript
                # Note: transcript_consequence is singular in gnomAD API
                tc = var.get("transcript_consequence")
                if not tc:
                    continue

                if not tc.get("is_canonical"):
                    continue

                # Check for missense variant
                consequence_terms = tc.get("consequence_terms", [])
                major_consequence = tc.get("major_consequence", "")

                if "missense_variant" not in consequence_terms and major_consequence != "missense_variant":
                    continue

                # Parse HGVS protein notation to get amino acid change
                hgvsp = tc.get("hgvsp")
                if not hgvsp:
                    continue

                aa_change = parse_hgvsp(hgvsp)
                if not aa_change:
                    continue

                ref_aa, aa_position, alt_aa = aa_change

                # Extract in silico predictor scores
                predictors = {}
                for pred in var.get("in_silico_predictors") or []:
                    pred_id = pred.get("id", "").lower()
                    pred_value = pred.get("value")
                    if pred_id and pred_value:
                        predictors[pred_id] = pred_value

                # Extract population frequencies from exome or genome (prefer exome)
                # Compute AF from ac/an since gnomAD v4 doesn't provide af per population
                pop_source = exome if exome.get("af") else genome
                pop_freqs = {}
                for pop in pop_source.get("populations") or []:
                    pop_id = pop.get("id", "").lower()
                    if pop_id:
                        ac = pop.get("ac")
                        an = pop.get("an")
                        # Calculate AF from ac/an
                        if ac is not None and an and an > 0:
                            pop_freqs[f"af_{pop_id}"] = ac / an
                        else:
                            pop_freqs[f"af_{pop_id}"] = None
                        pop_freqs[f"ac_{pop_id}"] = ac
                        pop_freqs[f"an_{pop_id}"] = an
                        pop_freqs[f"ac_hom_{pop_id}"] = pop.get("ac_hom")

                # Get faf95 and fafmax (filtering allele frequencies)
                faf95_data = exome.get("faf95") or genome.get("faf95") or {}
                fafmax_data = exome.get("fafmax") or genome.get("fafmax") or {}

                # Get rsids
                rsids = var.get("rsids") or []
                rsid = rsids[0] if rsids else None

                # Get filters
                exome_filters = exome.get("filters") or []
                genome_filters = genome.get("filters") or []

                all_variants.append({
                    # Basic variant info
                    "species": "human",
                    "gene_id": gene_id,
                    "gene_name": gene_symbol,
                    "transcript_id": tc.get("transcript_id"),
                    "variant_id": var.get("variant_id"),
                    "rsid": rsid,
                    "chromosome": var.get("variant_id", "").split("-")[0],
                    "genomic_position": var.get("pos"),
                    "ref_allele": var.get("ref"),
                    "alt_allele": var.get("alt"),
                    # Amino acid change
                    "aa_position": aa_position,
                    "ref_aa": ref_aa,
                    "alt_aa": alt_aa,
                    "hgvsp": hgvsp,
                    "hgvsc": tc.get("hgvsc"),
                    "consequence": major_consequence or "missense_variant",
                    # Global allele frequencies
                    "allele_frequency": max_af,
                    "exome_af": exome_af if exome_af else None,
                    "genome_af": genome_af if genome_af else None,
                    "exome_ac": exome.get("ac"),
                    "exome_an": exome.get("an"),
                    "exome_ac_hom": exome.get("ac_hom"),
                    "exome_ac_hemi": exome.get("ac_hemi"),
                    "genome_ac": genome.get("ac"),
                    "genome_an": genome.get("an"),
                    "genome_ac_hom": genome.get("ac_hom"),
                    "genome_ac_hemi": genome.get("ac_hemi"),
                    # Fafmax (filtering allele frequency max)
                    "faf95_max": fafmax_data.get("faf95_max"),
                    "faf95_max_population": fafmax_data.get("faf95_max_gen_anc"),
                    "faf99_max": fafmax_data.get("faf99_max"),
                    "faf99_max_population": fafmax_data.get("faf99_max_gen_anc"),
                    # Faf95 popmax
                    "faf95_popmax": faf95_data.get("popmax"),
                    "faf95_popmax_population": faf95_data.get("popmax_population"),
                    # Population-specific frequencies (computed from ac/an)
                    "af_afr": pop_freqs.get("af_afr"),
                    "af_amr": pop_freqs.get("af_amr"),
                    "af_asj": pop_freqs.get("af_asj"),
                    "af_eas": pop_freqs.get("af_eas"),
                    "af_fin": pop_freqs.get("af_fin"),
                    "af_nfe": pop_freqs.get("af_nfe"),
                    "af_sas": pop_freqs.get("af_sas"),
                    "af_mid": pop_freqs.get("af_mid"),
                    "af_ami": pop_freqs.get("af_ami"),
                    # In silico predictors / pathogenicity scores
                    "cadd_phred": predictors.get("cadd"),
                    "revel_score": predictors.get("revel_max") or predictors.get("revel"),
                    "spliceai_score": predictors.get("spliceai_ds_max") or predictors.get("spliceai"),
                    "alphamissense_score": predictors.get("alphamissense"),
                    "pangolin_score": predictors.get("pangolin_largest_ds") or predictors.get("pangolin"),
                    "phylop_score": predictors.get("phylop"),
                    "sift_max": predictors.get("sift_max") or predictors.get("sift"),
                    "polyphen_max": predictors.get("polyphen_max") or predictors.get("polyphen"),
                    # Transcript-level predictions
                    "polyphen_prediction": tc.get("polyphen_prediction"),
                    "sift_prediction": tc.get("sift_prediction"),
                    # LoF annotations
                    "lof": tc.get("lof"),
                    "lof_filter": tc.get("lof_filter"),
                    "lof_flags": tc.get("lof_flags"),
                    # Quality flags and filters
                    "flags": ",".join(var.get("flags") or []) if var.get("flags") else None,
                    "exome_filters": ",".join(exome_filters) if exome_filters else None,
                    "genome_filters": ",".join(genome_filters) if genome_filters else None,
                    # Source
                    "source": "gnomAD_v4"
                })
                count += 1

            print(f"{count} variants")
            time.sleep(0.5)  # Rate limiting
            
        except Exception as e:
            print(f"error: {e}")
    
    return pd.DataFrame(all_variants)


def download_gnomad_constraint_file() -> pd.DataFrame:
    """Download gnomAD gene constraint metrics (smaller file with gene-level summaries)."""
    
    url = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/constraint/gnomad.v4.0.constraint_metrics.tsv"
    
    print("Downloading gnomAD constraint metrics...")
    print("(This contains gene-level constraint scores, not individual variants)")
    
    try:
        df = pd.read_csv(url, sep="\t")
        print(f"Downloaded {len(df)} gene records")
        return df
    except Exception as e:
        print(f"Download failed: {e}")
        return pd.DataFrame()


# ============================================================================
# Method 4: Great Ape Genome Project (for chimp/gorilla population data)
# ============================================================================

def get_great_ape_variants_info():
    """Provide information about Great Ape Genome Project data."""
    
    info = """
Great Ape Genome Project Data Sources
=====================================

For chimpanzee and gorilla population-level variants, the main source is:

1. Great Ape Genome Project (GAGP)
   - Paper: Prado-Martinez et al., Nature 2013
   - VCF files: ftp://ftp.ebi.ac.uk/pub/databases/eva/PRJEB15086/
   - Includes: 79 great apes (bonobos, chimps, gorillas, orangutans)

2. Ensembl Variation
   - Some GAGP variants are loaded into Ensembl
   - Query via: /variation/{species}/{variant_id}
   - Population frequencies may be limited

3. NCBI dbSNP
   - Some great ape variants deposited
   - Search: https://www.ncbi.nlm.nih.gov/snp/

4. Recent studies:
   - Chimp: de Manuel et al., Science 2016 (population structure)
   - Gorilla: Xue et al., Nature 2015 (population genetics)

Download VCF files and process with:
    bcftools view -i 'AF>=0.10' file.vcf.gz | \\
    bcftools csq -f reference.fa -g annotation.gff3 - | \\
    grep 'missense'
"""
    print(info)


# ============================================================================
# Main
# ============================================================================

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Compile amino acid variants with >10% population frequency"
    )
    parser.add_argument(
        "--species", nargs="+",
        default=["human", "chimpanzee", "gorilla"],
        help="Species to analyze"
    )
    parser.add_argument(
        "--maf", type=float, default=0.10,
        help="Minimum allele frequency threshold"
    )
    parser.add_argument(
        "--gene-limit", type=int, default=None,
        help="Limit number of genes (for testing)"
    )
    parser.add_argument(
        "--genes", nargs="+", default=None,
        help="Specific gene symbols to query (for gnomAD method)"
    )
    parser.add_argument(
        "--output", default="aa_variants",
        help="Output file prefix"
    )
    parser.add_argument(
        "--method",
        choices=["ensembl", "biomart", "gnomad", "info"],
        default="ensembl",
        help="Data retrieval method"
    )
    
    args = parser.parse_args()
    
    print(f"="*60)
    print(f"Amino Acid Variants with MAF >= {args.maf}")
    print(f"Species: {', '.join(args.species)}")
    print(f"Method: {args.method}")
    print(f"="*60)
    
    if args.method == "info":
        get_great_ape_variants_info()
        return
    
    if args.method == "ensembl":
        df = compile_via_ensembl_api(
            species_list=args.species,
            gene_limit=args.gene_limit,
            maf_threshold=args.maf,
            output_prefix=args.output
        )
        
    elif args.method == "biomart":
        for species in args.species:
            df = fetch_via_biomart(
                species=species,
                maf_threshold=args.maf,
                output_file=f"{args.output}_{species}_biomart.csv"
            )
            
    elif args.method == "gnomad":
        if "human" not in args.species:
            print("gnomAD method only works for human")
            return

        if args.genes:
            genes = args.genes
        else:
            # Get all human genes first
            print("Fetching human gene list...")
            client = EnsemblClient()
            cache_file = ".cache/human_genes.pkl"
            gene_list = get_protein_coding_genes(client, "human", cache_file=cache_file)
            genes = [g["gene_name"] for g in gene_list if g.get("gene_name")]

            if args.gene_limit:
                genes = genes[:args.gene_limit]

        df = fetch_gnomad_via_api(genes, args.maf)

        if not df.empty:
            df.to_csv(f"{args.output}_human_gnomad.csv", index=False)
            print(f"Saved {len(df)} variants to {args.output}_human_gnomad.csv")
    
    print("\nDone!")


if __name__ == "__main__":
    main()
