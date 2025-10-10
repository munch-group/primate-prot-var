"""
Multiple approaches to retrieve loss-of-function variants and their population-specific 
allele frequencies from gnomAD database
"""

# ==============================================================================
# APPROACH 1: Using gnomAD GraphQL API directly
# ==============================================================================

import requests
import json
import pandas as pd
import time

def query_gnomad_graphql(gene_name, dataset="gnomad_r3"):
    """
    Query gnomAD GraphQL API for loss-of-function variants in a specific gene
    
    Parameters:
    -----------
    gene_name : str
        Gene symbol (e.g., 'BRCA1', 'TP53')
    dataset : str
        gnomAD dataset version ('gnomad_r2_1' for v2.1.1 or 'gnomad_r3' for v3)
    """
    
    # GraphQL endpoint
    url = "https://gnomad.broadinstitute.org/api"
    
    # GraphQL query for loss-of-function variants
    query = """
    query GeneVariants($geneSymbol: String!, $datasetId: DatasetId!) {
      gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
        variants(dataset: $datasetId) {
          variant_id
          pos
          ref
          alt
          consequence
          lof
          lof_filter
          exome {
            ac
            an
            af
            ac_popmax
            af_popmax
            popmax
            ac_afr
            af_afr
            ac_amr
            af_amr
            ac_asj
            af_asj
            ac_eas
            af_eas
            ac_fin
            af_fin
            ac_nfe
            af_nfe
            ac_sas
            af_sas
          }
          genome {
            ac
            an
            af
            ac_popmax
            af_popmax
            popmax
            ac_afr
            af_afr
            ac_amr
            af_amr
            ac_asj
            af_asj
            ac_eas
            af_eas
            ac_fin
            af_fin
            ac_nfe
            af_nfe
            ac_sas
            af_sas
          }
        }
      }
    }
    """
    
    variables = {
        "geneSymbol": gene_name,
        "datasetId": dataset
    }
    
    response = requests.post(
        url,
        json={'query': query, 'variables': variables},
        headers={'Content-Type': 'application/json'}
    )
    
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        print(f"Error: {response.status_code}")
        return None

def extract_lof_variants(gnomad_data):
    """
    Extract loss-of-function variants from gnomAD API response
    
    High-confidence LoF variants are those with:
    - consequence containing: 'stop_gained', 'frameshift', 'splice_donor', 'splice_acceptor'
    - lof = 'HC' (high confidence)
    """
    
    if not gnomad_data or 'data' not in gnomad_data:
        return pd.DataFrame()
    
    variants = gnomad_data['data']['gene']['variants']
    lof_variants = []
    
    lof_consequences = ['stop_gained', 'frameshift_variant', 
                       'splice_donor_variant', 'splice_acceptor_variant']
    
    for variant in variants:
        # Check if it's a loss-of-function variant
        if variant.get('lof') == 'HC' or any(cons in variant.get('consequence', '') 
                                             for cons in lof_consequences):
            
            # Extract population frequencies from exome or genome data
            freq_source = variant.get('exome') or variant.get('genome') or {}
            
            lof_variants.append({
                'variant_id': variant.get('variant_id'),
                'position': variant.get('pos'),
                'ref': variant.get('ref'),
                'alt': variant.get('alt'),
                'consequence': variant.get('consequence'),
                'lof_confidence': variant.get('lof'),
                'lof_filter': variant.get('lof_filter'),
                'af_total': freq_source.get('af'),
                'af_afr': freq_source.get('af_afr'),
                'af_amr': freq_source.get('af_amr'),
                'af_asj': freq_source.get('af_asj'),
                'af_eas': freq_source.get('af_eas'),
                'af_fin': freq_source.get('af_fin'),
                'af_nfe': freq_source.get('af_nfe'),
                'af_sas': freq_source.get('af_sas'),
                'ac_total': freq_source.get('ac'),
                'an_total': freq_source.get('an'),
                'popmax': freq_source.get('popmax'),
                'af_popmax': freq_source.get('af_popmax')
            })
    
    return pd.DataFrame(lof_variants)

# ==============================================================================
# APPROACH 2: Using pynoma package (requires installation)
# ==============================================================================

def get_lof_variants_pynoma(gene_name, gnomad_version=3):
    """
    Use pynoma package to retrieve variants (requires: pip install pynoma)
    
    Note: You need to clone and install pynoma first:
    git clone https://github.com/bioinfo-hcpa/pynoma.git
    pip install -e pynoma
    """
    try:
        from pynoma import GeneSearch
        
        gs = GeneSearch(gnomad_version, gene_name)
        df, clinical_df = gs.get_data()
        
        # Filter for LoF variants based on consequence
        lof_consequences = ['stop_gained', 'frameshift_variant', 
                           'splice_donor_variant', 'splice_acceptor_variant']
        
        lof_df = df[df['Consequence'].str.contains('|'.join(lof_consequences), na=False)]
        
        return lof_df
        
    except ImportError:
        print("pynoma not installed. Please install from GitHub:")
        print("git clone https://github.com/bioinfo-hcpa/pynoma.git")
        print("pip install -e pynoma")
        return None

# ==============================================================================
# APPROACH 3: Using gnomad-db SQLite database (requires download)
# ==============================================================================

def get_lof_variants_sqlite(variants_list, db_path=None):
    """
    Query gnomAD SQLite database for loss-of-function variants
    
    Requires downloading pre-built database from:
    https://zenodo.org/record/6818606/files/gnomad_db_v3.1.2.sqlite3.gz
    
    Or install gnomad-db: pip install gnomad-db
    """
    try:
        from gnomad_db.database import gnomAD_DB
        import sqlite3
        
        if not db_path:
            # Download database if not provided
            download_link = "https://zenodo.org/record/6818606/files/gnomad_db_v3.1.2.sqlite3.gz?download=1"
            output_dir = "./gnomad_db"
            gnomAD_DB.download_and_unzip(download_link, output_dir)
            db_path = f"{output_dir}/gnomad_db_v3.1.2.sqlite3"
        
        # Initialize database
        db = gnomAD_DB(db_path, gnomad_version="v3.1.2")
        
        # Query variants
        result = db.get_info_from_df(pd.DataFrame({'variant': variants_list}))
        
        # Filter for LoF based on consequence
        lof_result = result[result['most_severe_consequence'].str.contains(
            'stop_gained|frameshift|splice_donor|splice_acceptor', na=False)]
        
        return lof_result
        
    except ImportError:
        print("gnomad-db not installed. Install with: pip install gnomad-db")
        return None

# ==============================================================================
# APPROACH 4: Using tabix to query VCF files directly (streaming)
# ==============================================================================

def query_gnomad_tabix(chrom, start, end, build="hg38"):
    """
    Query gnomAD VCF files directly using tabix (requires pysam)
    
    This approach doesn't require downloading the entire dataset
    """
    try:
        import pysam
        
        # gnomAD public URLs for streaming
        if build == "hg38":
            url = "https://storage.googleapis.com/gnomad-public/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr{}.vcf.bgz".format(chrom)
        else:  # hg19/GRCh37
            url = "https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.{}.vcf.bgz".format(chrom)
        
        # Open remote VCF file
        vcf = pysam.TabixFile(url)
        
        lof_variants = []
        
        # Query region
        for row in vcf.fetch(str(chrom), start, end):
            fields = row.split('\t')
            info = fields[7]
            
            # Parse INFO field for LoF annotations and frequencies
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
            
            # Check if it's a LoF variant (look for VEP annotations)
            if 'vep' in info_dict and ('stop_gained' in info_dict['vep'] or 
                                       'frameshift' in info_dict['vep'] or
                                       'splice_donor' in info_dict['vep'] or
                                       'splice_acceptor' in info_dict['vep']):
                
                variant = {
                    'chrom': fields[0],
                    'pos': int(fields[1]),
                    'ref': fields[3],
                    'alt': fields[4],
                    'af': info_dict.get('AF', 'NA'),
                    'af_afr': info_dict.get('AF_afr', 'NA'),
                    'af_amr': info_dict.get('AF_amr', 'NA'),
                    'af_asj': info_dict.get('AF_asj', 'NA'),
                    'af_eas': info_dict.get('AF_eas', 'NA'),
                    'af_fin': info_dict.get('AF_fin', 'NA'),
                    'af_nfe': info_dict.get('AF_nfe', 'NA'),
                    'af_sas': info_dict.get('AF_sas', 'NA'),
                    'ac': info_dict.get('AC', 'NA'),
                    'an': info_dict.get('AN', 'NA')
                }
                lof_variants.append(variant)
        
        return pd.DataFrame(lof_variants)
        
    except ImportError:
        print("pysam not installed. Install with: pip install pysam")
        return None

# ==============================================================================
# MAIN FUNCTION: Batch processing multiple genes
# ==============================================================================

def get_all_lof_variants(gene_list, method='graphql', output_file='lof_variants.csv'):
    """
    Retrieve all loss-of-function variants for a list of genes
    
    Parameters:
    -----------
    gene_list : list
        List of gene symbols
    method : str
        Method to use ('graphql', 'pynoma', 'sqlite', 'tabix')
    output_file : str
        Output CSV filename
    """
    
    all_variants = []
    
    for gene in gene_list:
        print(f"Processing {gene}...")
        
        if method == 'graphql':
            # Use GraphQL API
            data = query_gnomad_graphql(gene, dataset="gnomad_r3")
            if data:
                df = extract_lof_variants(data)
                df['gene'] = gene
                all_variants.append(df)
                time.sleep(0.5)  # Rate limiting
                
        elif method == 'pynoma':
            # Use pynoma package
            df = get_lof_variants_pynoma(gene)
            if df is not None:
                df['gene'] = gene
                all_variants.append(df)
                
        # Add other methods as needed
    
    # Combine all results
    if all_variants:
        final_df = pd.concat(all_variants, ignore_index=True)
        
        # Save to CSV
        final_df.to_csv(output_file, index=False)
        print(f"Saved {len(final_df)} LoF variants to {output_file}")
        
        return final_df
    else:
        print("No variants found")
        return pd.DataFrame()

# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

if __name__ == "__main__":
    
    # Example 1: Query single gene using GraphQL
    print("Example 1: Query BRCA1 using GraphQL API")
    data = query_gnomad_graphql("BRCA1", dataset="gnomad_r3")
    if data:
        lof_df = extract_lof_variants(data)
        print(f"Found {len(lof_df)} LoF variants in BRCA1")
        print("\nFirst 5 variants:")
        print(lof_df.head())
        print("\nPopulation frequencies summary:")
        print(lof_df[['variant_id', 'consequence', 'af_total', 'af_afr', 
                      'af_amr', 'af_eas', 'af_nfe']].head())
    
    # Example 2: Batch processing multiple genes
    print("\n" + "="*50)
    print("Example 2: Batch processing multiple genes")
    genes = ["BRCA1", "BRCA2", "TP53", "PTEN", "APC"]
    results = get_all_lof_variants(genes, method='graphql', 
                                   output_file='cancer_genes_lof.csv')
    
    # Summary statistics
    if not results.empty:
        print(f"\nTotal LoF variants: {len(results)}")
        print(f"Genes analyzed: {results['gene'].nunique()}")
        print("\nVariants per gene:")
        print(results['gene'].value_counts())
        
        # Population with highest frequency
        pop_cols = [col for col in results.columns if col.startswith('af_') and 
                   col != 'af_total' and col != 'af_popmax']
        if pop_cols:
            print("\nAverage AF by population:")
            for pop in pop_cols:
                mean_af = results[pop].replace('NA', None).dropna().astype(float).mean()
                print(f"{pop}: {mean_af:.6f}")