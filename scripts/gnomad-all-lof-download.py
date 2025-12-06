"""
Script to download and process ALL loss-of-function variants from gnomAD
with population-specific allele frequencies
"""

import pandas as pd
import requests
import gzip
import io
from urllib.request import urlretrieve
import os

# ==============================================================================
# OPTION 1: Download pre-processed LoF constraint files
# ==============================================================================

def download_gnomad_constraint_data():
    """
    Download gnomAD gene constraint data which includes LoF metrics
    This is a smaller, pre-processed dataset
    """
    
    # URLs for gnomAD constraint files
    urls = {
        'v2.1.1': 'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz',
        'v3.1.2': 'https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/constraint/gnomad.v3.1.2.constraint_metrics.tsv'
    }
    
    constraint_data = {}
    
    for version, url in urls.items():
        print(f"Downloading gnomAD {version} constraint data...")
        
        if url.endswith('.bgz'):
            # Handle bgzipped file
            response = requests.get(url, stream=True)
            with gzip.open(io.BytesIO(response.content), 'rt') as f:
                df = pd.read_csv(f, sep='\t')
        else:
            # Handle regular TSV
            df = pd.read_csv(url, sep='\t')
        
        constraint_data[version] = df
        print(f"Downloaded {len(df)} transcript records for {version}")
    
    return constraint_data

# ==============================================================================
# OPTION 2: Download and process full VCF files (WARNING: Very large files!)
# ==============================================================================

def download_gnomad_vcf_header(chrom='1', version='3.1.2'):
    """
    Download just the header of gnomAD VCF to understand the format
    """
    
    if version.startswith('3'):
        url = f"https://storage.googleapis.com/gcp-public-data--gnomad/release/{version}/vcf/genomes/gnomad.genomes.v{version}.sites.chr{chrom}.vcf.bgz"
    else:
        url = f"https://storage.googleapis.com/gcp-public-data--gnomad/release/{version}/vcf/genomes/gnomad.genomes.r{version}.sites.{chrom}.vcf.bgz"
    
    print(f"Fetching VCF header from chromosome {chrom}...")
    
    # Download first few KB to get header
    headers = {'Range': 'bytes=0-50000'}
    response = requests.get(url, headers=headers)
    
    # Decompress and extract header lines
    content = gzip.decompress(response.content).decode('utf-8', errors='ignore')
    header_lines = [line for line in content.split('\n') if line.startswith('#')]
    
    return '\n'.join(header_lines)

def stream_process_gnomad_lof(chrom='21', version='3.1.2', max_variants=10000):
    """
    Stream process gnomAD VCF file to extract LoF variants
    
    WARNING: Full chromosome files are 10-100GB each!
    This example uses chromosome 21 (smallest) and limits to first N variants
    """
    
    try:
        import pysam
    except ImportError:
        print("pysam required for VCF streaming. Install with: pip install pysam")
        return None
    
    if version.startswith('3'):
        url = f"https://storage.googleapis.com/gcp-public-data--gnomad/release/{version}/vcf/genomes/gnomad.genomes.v{version}.sites.chr{chrom}.vcf.bgz"
    else:
        url = f"https://storage.googleapis.com/gcp-public-data--gnomad/release/{version}/vcf/exomes/gnomad.exomes.r{version}.sites.{chrom}.vcf.bgz"
    
    print(f"Streaming chromosome {chrom} from gnomAD {version}...")
    print("This will process only the first {:,} variants as an example".format(max_variants))
    
    vcf = pysam.VariantFile(url)
    
    lof_variants = []
    total_processed = 0
    
    # LoF consequences to look for in CSQ field
    lof_consequences = [
        'stop_gained',
        'frameshift_variant', 
        'splice_acceptor_variant',
        'splice_donor_variant',
        'stop_lost',
        'start_lost'
    ]
    
    for record in vcf:
        if total_processed >= max_variants:
            break
            
        total_processed += 1
        
        # Check VEP CSQ annotation for LoF
        if 'CSQ' in record.info:
            csq_string = str(record.info.get('CSQ', [''])[0])
            
            # Check if any LoF consequence is present
            is_lof = any(cons in csq_string for cons in lof_consequences)
            
            if is_lof:
                # Extract population frequencies
                variant_data = {
                    'chrom': record.chrom,
                    'pos': record.pos,
                    'id': record.id,
                    'ref': record.ref,
                    'alt': ','.join([str(a) for a in record.alts]) if record.alts else '',
                    'filter': ','.join(record.filter.keys()) if record.filter else 'PASS',
                    'af': record.info.get('AF', [None])[0],
                    'ac': record.info.get('AC', [None])[0],
                    'an': record.info.get('AN', None),
                    'af_afr': record.info.get('AF_afr', [None])[0],
                    'af_amr': record.info.get('AF_amr', [None])[0],
                    'af_asj': record.info.get('AF_asj', [None])[0],
                    'af_eas': record.info.get('AF_eas', [None])[0],
                    'af_fin': record.info.get('AF_fin', [None])[0],
                    'af_nfe': record.info.get('AF_nfe', [None])[0],
                    'af_oth': record.info.get('AF_oth', [None])[0],
                    'af_sas': record.info.get('AF_sas', [None])[0],
                    'csq': csq_string[:200]  # First 200 chars of consequence
                }
                
                lof_variants.append(variant_data)
        
        if total_processed % 1000 == 0:
            print(f"Processed {total_processed:,} variants, found {len(lof_variants)} LoF...")
    
    vcf.close()
    
    df = pd.DataFrame(lof_variants)
    print(f"\nFound {len(df)} LoF variants in first {total_processed:,} variants of chr{chrom}")
    
    return df

# ==============================================================================
# OPTION 3: Use pre-built SQLite databases (Recommended for full dataset)
# ==============================================================================

def setup_gnomad_sqlite_db(version='v3.1.2'):
    """
    Download and setup gnomAD SQLite database for efficient querying
    This is the recommended approach for working with the full dataset
    """
    
    print("Setting up gnomAD SQLite database...")
    print("This will download a ~46GB compressed file (98GB uncompressed)")
    print("Make sure you have enough disk space!")
    
    db_urls = {
        'v4.0_wgs': 'https://zenodo.org/records/10066323/files/gnomad_db_wgs_v4.0.sqlite3.gz',
        'v4.0_wes': 'https://zenodo.org/records/10066310/files/gnomad_db_wes_v4.0.sqlite3.gz',
        'v3.1.2': 'https://zenodo.org/record/6818606/files/gnomad_db_v3.1.2.sqlite3.gz',
        'v2.1.1': 'https://zenodo.org/record/5770384/files/gnomad_db_v2.1.1.sqlite3.gz'
    }
    
    if version not in db_urls:
        print(f"Version {version} not available. Choose from: {list(db_urls.keys())}")
        return None
    
    url = db_urls[version]
    filename = url.split('/')[-1].replace('?download=1', '')
    
    # Check if already downloaded
    if os.path.exists(filename.replace('.gz', '')):
        print(f"Database already exists: {filename.replace('.gz', '')}")
        return filename.replace('.gz', '')
    
    print(f"Downloading from {url}...")
    urlretrieve(url, filename)
    
    print("Decompressing...")
    import gzip
    import shutil
    
    with gzip.open(filename, 'rb') as f_in:
        with open(filename.replace('.gz', ''), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    # Remove compressed file to save space
    os.remove(filename)
    
    print(f"Database ready: {filename.replace('.gz', '')}")
    return filename.replace('.gz', '')

def query_all_lof_from_sqlite(db_path, output_file='all_lof_variants.csv', 
                              batch_size=100000):
    """
    Query all LoF variants from SQLite database with population frequencies
    """
    
    import sqlite3
    
    print(f"Connecting to database: {db_path}")
    conn = sqlite3.connect(db_path)
    
    # Query for all LoF variants
    # Note: Column names may vary by version, adjust as needed
    query = """
    SELECT 
        chrom,
        pos,
        ref,
        alt,
        variant_id,
        most_severe_consequence,
        AF,
        AF_afr,
        AF_amr,
        AF_asj,
        AF_eas,
        AF_fin,
        AF_nfe,
        AF_oth,
        AF_sas,
        AC,
        AN,
        popmax,
        AF_popmax
    FROM variants
    WHERE most_severe_consequence LIKE '%stop_gained%'
       OR most_severe_consequence LIKE '%frameshift%'
       OR most_severe_consequence LIKE '%splice_acceptor%'
       OR most_severe_consequence LIKE '%splice_donor%'
       OR most_severe_consequence LIKE '%start_lost%'
       OR most_severe_consequence LIKE '%stop_lost%'
    """
    
    print("Querying all LoF variants...")
    print("This may take several minutes for the full dataset...")
    
    # Read in chunks to manage memory
    chunks = []
    for chunk in pd.read_sql_query(query, conn, chunksize=batch_size):
        chunks.append(chunk)
        print(f"Processed {len(chunks) * batch_size:,} variants...")
    
    conn.close()
    
    # Combine all chunks
    df = pd.concat(chunks, ignore_index=True)
    
    print(f"\nTotal LoF variants found: {len(df):,}")
    
    # Save to CSV
    print(f"Saving to {output_file}...")
    df.to_csv(output_file, index=False)
    
    return df

# ==============================================================================
# OPTION 4: Download pre-filtered LoF files (if available)
# ==============================================================================

def download_preprocessed_lof():
    """
    Download pre-processed files containing high-confidence LoF variants
    """
    
    # Download pLoF (predicted loss-of-function) variants file
    # These are variants that passed LOFTEE filters
    
    urls = {
        'pLoF_v2': 'https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/gnomad.v2.1.1.all_lofs.txt.bgz',
        'homozygous_lof': 'https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/gnomad.v2.1.1.lof_homozygous.txt.bgz'
    }
    
    dfs = {}
    
    for name, url in urls.items():
        print(f"Downloading {name}...")
        response = requests.get(url)
        
        with gzip.open(io.BytesIO(response.content), 'rt') as f:
            df = pd.read_csv(f, sep='\t')
            dfs[name] = df
            print(f"Downloaded {len(df)} variants")
    
    return dfs

# ==============================================================================
# MAIN EXECUTION EXAMPLE
# ==============================================================================

if __name__ == "__main__":
    
    print("gnomAD Loss-of-Function Variants Extraction")
    print("=" * 60)
    
    # Example 1: Download constraint data (quick, small files)
    print("\n1. Downloading constraint data...")
    constraint_data = download_gnomad_constraint_data()
    
    # Example 2: Stream process a small sample from VCF
    print("\n2. Streaming sample LoF variants from chromosome 21...")
    sample_lof = stream_process_gnomad_lof(chrom='21', max_variants=5000)
    
    if sample_lof is not None and not sample_lof.empty:
        print("\nSample of LoF variants with population frequencies:")
        print(sample_lof[['chrom', 'pos', 'ref', 'alt', 'af', 
                         'af_afr', 'af_eas', 'af_nfe']].head(10))
        
        # Save sample
        sample_lof.to_csv('sample_lof_chr21.csv', index=False)
        print("\nSaved sample to sample_lof_chr21.csv")
    
    # Example 3: For full dataset, use SQLite (recommended)
    print("\n3. For full dataset analysis, use SQLite database:")
    print("   Uncomment the lines below to download and query the full dataset")
    print("   WARNING: This requires ~100GB of disk space!")
    
    # Uncomment to download and use full dataset:
    # db_path = setup_gnomad_sqlite_db(version='v3.1.2')
    # if db_path:
    #     all_lof = query_all_lof_from_sqlite(db_path, 'all_gnomad_lof.csv')
    #     print(f"\nSummary by population:")
    #     pop_cols = [col for col in all_lof.columns if col.startswith('AF_')]
    #     for col in pop_cols:
    #         mean_af = all_lof[col].mean()
    #         print(f"{col}: {mean_af:.6f}")
    
    print("\n" + "=" * 60)
    print("Processing complete!")
    print("\nNotes:")
    print("- For small queries: Use GraphQL API or streaming VCF")
    print("- For full dataset: Download SQLite database or VCF files")
    print("- Population codes: AFR=African, AMR=American, ASJ=Ashkenazi Jewish,")
    print("                   EAS=East Asian, FIN=Finnish, NFE=Non-Finnish European,")
    print("                   SAS=South Asian, OTH=Other")
