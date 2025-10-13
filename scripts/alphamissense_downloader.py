import pandas as pd
import requests
import gzip
import io
from pathlib import Path


# # Download all files to default directory
# python alphamissense_downloader.py

# # Download to specific directory  
# python alphamissense_downloader.py --dir /path/to/data

# # Only amino acid substitutions (most comprehensive)
# python alphamissense_downloader.py --files AlphaMissense_aa_substitutions.tsv.gz

# # Only genomic variants for hg38
# python alphamissense_downloader.py --files AlphaMissense_hg38.tsv.gz

# # Use zenodo_get tool (good for large datasets)
# python alphamissense_downloader.py --method zenodo_get


# from alphamissense_downloader import AlphaMissenseDownloader

# # Download everything
# downloader = AlphaMissenseDownloader("my_data_dir")
# downloader.download_all()

# # Or use convenience functions
# download_amino_acid_substitutions_only("my_data_dir")
# download_genomic_variants_only("my_data_dir", genome_build="hg38")  
# quick_download_gene_averages("my_data_dir")

def download_ttll10_alphamissense_data():
    """
    Download and extract AlphaMissense pathogenicity scores for TTLL10 (Q6ZVM7)
    """
    
    # TTLL10 UniProt ID
    ttll10_uniprot = "Q6ZVM7"
    
    # Method 1: Download from Zenodo (Main dataset)
    print("Downloading AlphaMissense amino acid substitutions dataset...")
    
    # Zenodo direct download URL for amino acid substitutions file
    zenodo_url = "https://zenodo.org/record/8208688/files/AlphaMissense_aa_substitutions.tsv.gz"
    
    # Download the compressed file
    response = requests.get(zenodo_url, stream=True)
    response.raise_for_status()
    
    print("Processing data...")
    
    # Read compressed file directly into pandas
    with gzip.open(io.BytesIO(response.content), 'rt') as f:
        # Read header first to understand structure
        header = f.readline().strip().split('\t')
        print(f"Columns: {header}")
        
        # Read data in chunks to filter for TTLL10
        ttll10_data = []
        chunk_size = 10000
        
        while True:
            chunk_lines = []
            for _ in range(chunk_size):
                line = f.readline()
                if not line:
                    break
                chunk_lines.append(line.strip().split('\t'))
            
            if not chunk_lines:
                break
                
            # Convert to DataFrame chunk
            chunk_df = pd.DataFrame(chunk_lines, columns=header)
            
            # Filter for TTLL10
            ttll10_chunk = chunk_df[chunk_df['uniprot_id'] == ttll10_uniprot]
            
            if not ttll10_chunk.empty:
                ttll10_data.append(ttll10_chunk)
                print(f"Found {len(ttll10_chunk)} TTLL10 variants in this chunk")
    
    # Combine all TTLL10 data
    if ttll10_data:
        ttll10_df = pd.concat(ttll10_data, ignore_index=True)
        print(f"Total TTLL10 variants found: {len(ttll10_df)}")
        
        # Convert pathogenicity score to float
        ttll10_df['alphamissense_pathogenicity'] = ttll10_df['alphamissense_pathogenicity'].astype(float)
        
        # Save to CSV
        output_file = "TTLL10_AlphaMissense_pathogenicity.csv"
        ttll10_df.to_csv(output_file, index=False)
        print(f"Data saved to: {output_file}")
        
        return ttll10_df
    else:
        print("No TTLL10 data found!")
        return None

def download_using_zenodo_get():
    """
    Alternative method using zenodo_get tool
    """
    import subprocess
    import os
    
    # Install zenodo_get if not available
    try:
        import zenodo_get
    except ImportError:
        print("Installing zenodo_get...")
        subprocess.check_call(["pip", "install", "zenodo_get"])
    
    # Download specific file
    zenodo_record = "8208688"
    filename = "AlphaMissense_aa_substitutions.tsv.gz"
    
    cmd = f"zenodo_get {zenodo_record} -g {filename}"
    subprocess.run(cmd, shell=True)
    
    print(f"Downloaded {filename}")

def extract_ttll10_from_local_file(filename="AlphaMissense_aa_substitutions.tsv.gz"):
    """
    Extract TTLL10 data from locally downloaded file
    """
    ttll10_uniprot = "Q6ZVM7"
    
    print(f"Extracting TTLL10 data from {filename}...")
    
    # Read and filter the compressed file
    ttll10_variants = []
    
    with gzip.open(filename, 'rt') as f:
        header = f.readline().strip().split('\t')
        
        for line in f:
            fields = line.strip().split('\t')
            if fields[0] == ttll10_uniprot:  # uniprot_id is first column
                ttll10_variants.append(fields)
    
    # Create DataFrame
    ttll10_df = pd.DataFrame(ttll10_variants, columns=header)
    ttll10_df['alphamissense_pathogenicity'] = ttll10_df['alphamissense_pathogenicity'].astype(float)
    
    # Save results
    output_file = "TTLL10_AlphaMissense_pathogenicity.csv"
    ttll10_df.to_csv(output_file, index=False)
    
    print(f"Extracted {len(ttll10_df)} TTLL10 variants")
    print(f"Data saved to: {output_file}")
    
    return ttll10_df

def download_from_hegelab_api():
    """
    Alternative: Download from HegedLab web interface API
    """
    # This uses the web interface at alphamissense.hegelab.org
    api_url = f"https://alphamissense.hegelab.org/api/protein/Q6ZVM7"
    
    try:
        response = requests.get(api_url)
        if response.status_code == 200:
            data = response.json()
            
            # Convert to DataFrame
            df = pd.DataFrame(data)
            df.to_csv("TTLL10_AlphaMissense_hegelab.csv", index=False)
            print("Data downloaded from HegedLab API")
            return df
        else:
            print(f"API request failed with status {response.status_code}")
    except Exception as e:
        print(f"API request failed: {e}")
    
    return None

# Comprehensive analysis function
def analyze_ttll10_pathogenicity(df):
    """
    Analyze the pathogenicity scores for TTLL10
    """
    print("\n=== TTLL10 AlphaMissense Analysis ===")
    print(f"Total variants: {len(df)}")
    
    # Categorize by pathogenicity
    likely_benign = df[df['alphamissense_pathogenicity'] < 0.34]
    ambiguous = df[(df['alphamissense_pathogenicity'] >= 0.34) & 
                   (df['alphamissense_pathogenicity'] <= 0.564)]
    likely_pathogenic = df[df['alphamissense_pathogenicity'] > 0.564]
    
    print(f"Likely benign: {len(likely_benign)} ({len(likely_benign)/len(df)*100:.1f}%)")
    print(f"Ambiguous: {len(ambiguous)} ({len(ambiguous)/len(df)*100:.1f}%)")
    print(f"Likely pathogenic: {len(likely_pathogenic)} ({len(likely_pathogenic)/len(df)*100:.1f}%)")
    
    # Find high-risk positions
    print(f"\nMean pathogenicity score: {df['alphamissense_pathogenicity'].mean():.3f}")
    print(f"Highest risk variants (score > 0.8):")
    
    high_risk = df[df['alphamissense_pathogenicity'] > 0.8].sort_values('alphamissense_pathogenicity', ascending=False)
    if not high_risk.empty:
        for _, row in high_risk.head(10).iterrows():
            print(f"  Position {row['position']}: {row['ref_aa']} → {row['alt_aa']} (score: {row['alphamissense_pathogenicity']:.3f})")
    
    return {
        'likely_benign': likely_benign,
        'ambiguous': ambiguous, 
        'likely_pathogenic': likely_pathogenic
    }

# Main execution
if __name__ == "__main__":
    print("Downloading TTLL10 AlphaMissense data...")
    
    # Try different methods
    try:
        # Method 1: Direct download and filter
        ttll10_data = download_ttll10_alphamissense_data()
        
        if ttll10_data is not None:
            # Analyze the data
            analysis = analyze_ttll10_pathogenicity(ttll10_data)
            
            print("\nFirst few variants:")
            print(ttll10_data.head())
            
        else:
            print("Trying alternative download method...")
            # Method 2: Use zenodo_get
            download_using_zenodo_get()
            ttll10_data = extract_ttll10_from_local_file()
            
    except Exception as e:
        print(f"Download failed: {e}")
        print("You may need to manually download from:")
        print("https://zenodo.org/record/8208688/files/AlphaMissense_aa_substitutions.tsv.gz")