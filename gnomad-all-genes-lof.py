"""
Scalable pipeline to retrieve loss-of-function variants for ALL human protein-coding genes
with gene names, CDS positions, mutation types, and population frequencies
"""

import pandas as pd
import numpy as np
import sqlite3
import requests
import gzip
import io
import os
import time
from multiprocessing import Pool, cpu_count
from typing import List, Dict, Tuple
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ==============================================================================
# APPROACH 1: Download and process gnomAD annotated files (RECOMMENDED)
# ==============================================================================

class GnomADLoFProcessor:
    """
    Main class for processing all LoF variants across all protein-coding genes
    """
    
    def __init__(self, version='v2.1.1', genome_build='GRCh37'):
        self.version = version
        self.genome_build = genome_build
        self.lof_categories = {
            'stop_gained': 'in-frame stop',
            'stop_lost': 'stop lost',
            'start_lost': 'start lost',
            'frameshift_variant': 'frameshift',
            'splice_acceptor_variant': 'splice acceptor',
            'splice_donor_variant': 'splice donor',
            'transcript_ablation': 'transcript ablation',
            'transcript_amplification': 'transcript amplification'
        }
        
    def download_annotated_vcf(self, chrom):
        """
        Download annotated VCF with VEP annotations including gene names and CDS positions
        """
        if self.version == 'v2.1.1':
            base_url = "https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/"
            file_name = f"gnomad.exomes.r2.1.1.sites.{chrom}.vcf.bgz"
        elif self.version == 'v3.1.2':
            base_url = "https://storage.googleapis.com/gnomad-public/release/3.1.2/vcf/genomes/"
            file_name = f"gnomad.genomes.v3.1.2.sites.chr{chrom}.vcf.bgz"
        else:
            raise ValueError(f"Version {self.version} not supported")
        
        url = base_url + file_name
        logger.info(f"Downloading chromosome {chrom} from {url}")
        
        # Download to local file for processing
        local_file = f"gnomad_chr{chrom}.vcf.bgz"
        if not os.path.exists(local_file):
            os.system(f"wget -q {url} -O {local_file}")
        
        return local_file
    
    def parse_vep_annotation(self, csq_string, csq_header):
        """
        Parse VEP CSQ annotation to extract gene, position, and consequence
        """
        if not csq_string:
            return []
        
        # Parse CSQ header to get field positions
        fields = csq_header.replace('Consequence annotations from Ensembl VEP. Format: ', '').split('|')
        field_idx = {field: i for i, field in enumerate(fields)}
        
        annotations = []
        for transcript_annotation in csq_string.split(','):
            values = transcript_annotation.split('|')
            
            if len(values) != len(fields):
                continue
            
            # Extract relevant fields
            try:
                annotation = {
                    'gene_symbol': values[field_idx.get('SYMBOL', -1)] if 'SYMBOL' in field_idx else '',
                    'gene_id': values[field_idx.get('Gene', -1)] if 'Gene' in field_idx else '',
                    'transcript_id': values[field_idx.get('Feature', -1)] if 'Feature' in field_idx else '',
                    'consequence': values[field_idx.get('Consequence', -1)] if 'Consequence' in field_idx else '',
                    'cds_position': values[field_idx.get('CDS_position', -1)] if 'CDS_position' in field_idx else '',
                    'protein_position': values[field_idx.get('Protein_position', -1)] if 'Protein_position' in field_idx else '',
                    'amino_acids': values[field_idx.get('Amino_acids', -1)] if 'Amino_acids' in field_idx else '',
                    'codons': values[field_idx.get('Codons', -1)] if 'Codons' in field_idx else '',
                    'impact': values[field_idx.get('IMPACT', -1)] if 'IMPACT' in field_idx else '',
                    'lof': values[field_idx.get('LoF', -1)] if 'LoF' in field_idx else '',
                    'lof_filter': values[field_idx.get('LoF_filter', -1)] if 'LoF_filter' in field_idx else '',
                    'lof_flags': values[field_idx.get('LoF_flags', -1)] if 'LoF_flags' in field_idx else ''
                }
                
                # Only keep if it's a LoF variant
                if any(lof_type in annotation['consequence'] for lof_type in self.lof_categories.keys()):
                    annotations.append(annotation)
                    
            except IndexError:
                continue
        
        return annotations

# ==============================================================================
# APPROACH 2: Use pre-built Hail tables (Most efficient for full dataset)
# ==============================================================================

def setup_hail_environment():
    """
    Set up Hail for processing gnomAD data efficiently
    Hail is the framework used by gnomAD internally
    """
    try:
        import hail as hl
        
        # Initialize Hail
        hl.init(default_reference='GRCh38', log='hail.log')
        
        return hl
    except ImportError:
        logger.error("Hail not installed. Install with: pip install hail")
        return None

def process_with_hail(output_file='all_genes_lof.csv'):
    """
    Process all LoF variants using Hail (most efficient for full dataset)
    """
    hl = setup_hail_environment()
    if not hl:
        return None
    
    # Load gnomAD Hail table
    # You can download these from: https://gnomad.broadinstitute.org/downloads
    if os.path.exists('gnomad.exomes.r2.1.1.sites.ht'):
        ht = hl.read_table('gnomad.exomes.r2.1.1.sites.ht')
    else:
        # Download from gnomAD
        logger.info("Downloading gnomAD Hail table...")
        ht = hl.read_table('gs://gnomad-public-requester-pays/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht')
    
    # Filter for LoF variants
    lof_consequences = [
        'stop_gained',
        'frameshift_variant',
        'splice_acceptor_variant',
        'splice_donor_variant',
        'start_lost',
        'stop_lost'
    ]
    
    # Filter variants
    ht_lof = ht.filter(
        hl.any(lambda csq: hl.any(
            lambda term: term.lower().contains(csq),
            ht.vep.most_severe_consequence.lower()
        ), lof_consequences)
    )
    
    # Select relevant fields
    ht_lof = ht_lof.select(
        gene_symbol=ht_lof.vep.worst_csq_by_gene.gene_symbol,
        gene_id=ht_lof.vep.worst_csq_by_gene.gene_id,
        transcript_id=ht_lof.vep.worst_csq_by_gene.transcript_id,
        consequence=ht_lof.vep.most_severe_consequence,
        cds_position=ht_lof.vep.worst_csq_by_gene.cds_position,
        protein_position=ht_lof.vep.worst_csq_by_gene.protein_position,
        lof=ht_lof.vep.worst_csq_by_gene.lof,
        lof_filter=ht_lof.vep.worst_csq_by_gene.lof_filter,
        af=ht_lof.freq[0].AF,
        af_afr=ht_lof.freq[1].AF,  # Adjust indices based on population order
        af_amr=ht_lof.freq[2].AF,
        af_asj=ht_lof.freq[3].AF,
        af_eas=ht_lof.freq[4].AF,
        af_fin=ht_lof.freq[5].AF,
        af_nfe=ht_lof.freq[6].AF,
        af_sas=ht_lof.freq[7].AF,
        ac=ht_lof.freq[0].AC,
        an=ht_lof.freq[0].AN
    )
    
    # Export to pandas
    df = ht_lof.to_pandas()
    
    # Save to file
    df.to_csv(output_file, index=False)
    logger.info(f"Saved {len(df)} LoF variants to {output_file}")
    
    return df

# ==============================================================================
# APPROACH 3: Download and process pre-annotated files
# ==============================================================================

def download_all_genes_list():
    """
    Download list of all human protein-coding genes from Ensembl
    """
    # Using Ensembl BioMart REST API
    url = "http://www.ensembl.org/biomart/martservice"
    
    query = '''<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0">
        <Dataset name="hsapiens_gene_ensembl" interface="default">
            <Attribute name="ensembl_gene_id"/>
            <Attribute name="hgnc_symbol"/>
            <Attribute name="chromosome_name"/>
            <Attribute name="start_position"/>
            <Attribute name="end_position"/>
            <Attribute name="gene_biotype"/>
        </Dataset>
    </Query>'''
    
    response = requests.post(url, data={'query': query})
    
    # Parse response
    df_genes = pd.read_csv(io.StringIO(response.text), sep='\t')
    
    # Filter for protein-coding genes on standard chromosomes
    standard_chr = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
    df_genes = df_genes[
        (df_genes['Gene type'] == 'protein_coding') &
        (df_genes['Chromosome/scaffold name'].isin(standard_chr))
    ]
    
    logger.info(f"Found {len(df_genes)} protein-coding genes")
    
    return df_genes

def process_gnomad_annotated_files(chromosomes=None):
    """
    Process pre-annotated gnomAD files that include gene annotations
    """
    if chromosomes is None:
        chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
    
    all_lof_variants = []
    
    for chrom in chromosomes:
        logger.info(f"Processing chromosome {chrom}")
        
        # Download annotated file with consequences
        url = f"https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.{chrom}.vcf.bgz"
        
        try:
            import pysam
            vcf = pysam.VariantFile(url)
            
            # Get CSQ header
            csq_header = vcf.header.info['CSQ'].description
            
            lof_variants = []
            count = 0
            
            for record in vcf:
                count += 1
                if count % 10000 == 0:
                    logger.info(f"Chr{chrom}: Processed {count:,} variants, found {len(lof_variants)} LoF")
                
                # Parse VEP annotations
                if 'CSQ' in record.info:
                    csq_string = record.info['CSQ'][0] if record.info['CSQ'] else ''
                    
                    # Check for LoF consequences
                    lof_types = ['stop_gained', 'frameshift', 'splice_acceptor', 'splice_donor', 'start_lost']
                    
                    if any(lof_type in csq_string for lof_type in lof_types):
                        # Parse annotations
                        processor = GnomADLoFProcessor()
                        annotations = processor.parse_vep_annotation(csq_string, csq_header)
                        
                        for ann in annotations:
                            if ann['gene_symbol']:  # Only keep if gene symbol exists
                                variant_data = {
                                    'gene_name': ann['gene_symbol'],
                                    'gene_id': ann['gene_id'],
                                    'transcript_id': ann['transcript_id'],
                                    'chrom': record.chrom,
                                    'pos': record.pos,
                                    'ref': record.ref,
                                    'alt': ','.join([str(a) for a in record.alts]) if record.alts else '',
                                    'variant_id': f"{record.chrom}-{record.pos}-{record.ref}-{record.alts[0]}" if record.alts else '',
                                    'mutation_type': processor.lof_categories.get(
                                        ann['consequence'].split('&')[0], 
                                        ann['consequence']
                                    ),
                                    'consequence': ann['consequence'],
                                    'cds_position': ann['cds_position'],
                                    'protein_position': ann['protein_position'],
                                    'codons': ann['codons'],
                                    'amino_acids': ann['amino_acids'],
                                    'lof_confidence': ann['lof'],
                                    'lof_filter': ann['lof_filter'],
                                    'af': record.info.get('AF', [None])[0],
                                    'af_afr': record.info.get('AF_afr', [None])[0],
                                    'af_amr': record.info.get('AF_amr', [None])[0],
                                    'af_asj': record.info.get('AF_asj', [None])[0],
                                    'af_eas': record.info.get('AF_eas', [None])[0],
                                    'af_fin': record.info.get('AF_fin', [None])[0],
                                    'af_nfe': record.info.get('AF_nfe', [None])[0],
                                    'af_sas': record.info.get('AF_sas', [None])[0],
                                    'ac': record.info.get('AC', [None])[0],
                                    'an': record.info.get('AN', None)
                                }
                                lof_variants.append(variant_data)
            
            vcf.close()
            all_lof_variants.extend(lof_variants)
            logger.info(f"Chromosome {chrom} complete: {len(lof_variants)} LoF variants")
            
        except Exception as e:
            logger.error(f"Error processing chromosome {chrom}: {e}")
            continue
    
    # Create DataFrame
    df = pd.DataFrame(all_lof_variants)
    
    return df

# ==============================================================================
# APPROACH 4: Parallel processing using multiple cores
# ==============================================================================

def process_chromosome_parallel(args):
    """
    Process a single chromosome for parallel execution
    """
    chrom, version = args
    processor = GnomADLoFProcessor(version=version)
    
    try:
        # Process chromosome
        df = process_single_chromosome(chrom, processor)
        
        # Save intermediate result
        output_file = f"lof_chr{chrom}.parquet"
        df.to_parquet(output_file, index=False)
        
        return f"Chromosome {chrom}: {len(df)} LoF variants"
    except Exception as e:
        return f"Error on chromosome {chrom}: {e}"

def process_single_chromosome(chrom, processor):
    """
    Process a single chromosome and extract LoF variants
    """
    import pysam
    
    # Determine URL based on version
    if processor.version == 'v2.1.1':
        url = f"https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.{chrom}.vcf.bgz"
    else:
        url = f"https://storage.googleapis.com/gnomad-public/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr{chrom}.vcf.bgz"
    
    vcf = pysam.VariantFile(url)
    
    # Get CSQ header
    csq_header = vcf.header.info.get('CSQ', None)
    if csq_header:
        csq_header = csq_header.description
    
    lof_variants = []
    
    for record in vcf:
        if 'CSQ' in record.info:
            csq_string = str(record.info.get('CSQ', [''])[0])
            
            # Quick check for LoF consequences
            if any(lof in csq_string for lof in processor.lof_categories.keys()):
                annotations = processor.parse_vep_annotation(csq_string, csq_header)
                
                for ann in annotations:
                    if ann['gene_symbol']:
                        variant_data = {
                            'gene_name': ann['gene_symbol'],
                            'chrom': record.chrom,
                            'pos': record.pos,
                            'ref': record.ref,
                            'alt': ','.join([str(a) for a in record.alts]) if record.alts else '',
                            'mutation_type': processor.lof_categories.get(
                                ann['consequence'].split('&')[0],
                                ann['consequence']
                            ),
                            'cds_position': ann['cds_position'],
                            'protein_position': ann['protein_position'],
                            'consequence': ann['consequence'],
                            'af': record.info.get('AF', [None])[0],
                            'af_afr': record.info.get('AF_afr', [None])[0],
                            'af_amr': record.info.get('AF_amr', [None])[0],
                            'af_asj': record.info.get('AF_asj', [None])[0],
                            'af_eas': record.info.get('AF_eas', [None])[0],
                            'af_fin': record.info.get('AF_fin', [None])[0],
                            'af_nfe': record.info.get('AF_nfe', [None])[0],
                            'af_sas': record.info.get('AF_sas', [None])[0]
                        }
                        lof_variants.append(variant_data)
    
    vcf.close()
    
    return pd.DataFrame(lof_variants)

def run_parallel_processing(version='v2.1.1', n_cores=None):
    """
    Run parallel processing across all chromosomes
    """
    chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
    
    if n_cores is None:
        n_cores = min(cpu_count() - 1, len(chromosomes))
    
    logger.info(f"Processing {len(chromosomes)} chromosomes using {n_cores} cores")
    
    # Prepare arguments
    args = [(chrom, version) for chrom in chromosomes]
    
    # Run parallel processing
    with Pool(n_cores) as pool:
        results = pool.map(process_chromosome_parallel, args)
    
    # Print results
    for result in results:
        logger.info(result)
    
    # Combine all parquet files
    logger.info("Combining results...")
    dfs = []
    for chrom in chromosomes:
        file_path = f"lof_chr{chrom}.parquet"
        if os.path.exists(file_path):
            df = pd.read_parquet(file_path)
            dfs.append(df)
            os.remove(file_path)  # Clean up intermediate file
    
    final_df = pd.concat(dfs, ignore_index=True)
    
    return final_df

# ==============================================================================
# MAIN EXECUTION PIPELINE
# ==============================================================================

def create_final_summary_table(df):
    """
    Create a summary table with key statistics per gene
    """
    summary = df.groupby('gene_name').agg({
        'variant_id': 'count',
        'mutation_type': lambda x: x.value_counts().to_dict(),
        'af': 'mean',
        'af_afr': 'mean',
        'af_amr': 'mean',
        'af_asj': 'mean',
        'af_eas': 'mean',
        'af_fin': 'mean',
        'af_nfe': 'mean',
        'af_sas': 'mean'
    }).reset_index()
    
    summary.columns = [
        'gene_name',
        'total_lof_variants',
        'mutation_types',
        'mean_af',
        'mean_af_afr',
        'mean_af_amr',
        'mean_af_asj',
        'mean_af_eas',
        'mean_af_fin',
        'mean_af_nfe',
        'mean_af_sas'
    ]
    
    return summary

def main():
    """
    Main pipeline to process all human protein-coding genes
    """
    
    logger.info("="*60)
    logger.info("gnomAD Loss-of-Function Analysis Pipeline")
    logger.info("Processing ALL human protein-coding genes")
    logger.info("="*60)
    
    # Choose processing method based on what's available
    method = 'parallel'  # Options: 'hail', 'parallel', 'sequential'
    
    if method == 'hail':
        # Method 1: Use Hail (most efficient for full dataset)
        logger.info("Using Hail for processing...")
        df = process_with_hail('all_genes_lof.csv')
        
    elif method == 'parallel':
        # Method 2: Parallel processing of VCF files
        logger.info("Using parallel processing...")
        df = run_parallel_processing(version='v2.1.1', n_cores=4)
        
    else:
        # Method 3: Sequential processing (slower but simpler)
        logger.info("Using sequential processing...")
        df = process_gnomad_annotated_files()
    
    if df is not None and not df.empty:
        # Save full results
        output_file = f"gnomad_all_genes_lof_{datetime.now().strftime('%Y%m%d')}.csv"
        df.to_csv(output_file, index=False)
        logger.info(f"Saved {len(df)} LoF variants to {output_file}")
        
        # Create summary table
        summary = create_final_summary_table(df)
        summary_file = f"gnomad_lof_summary_by_gene_{datetime.now().strftime('%Y%m%d')}.csv"
        summary.to_csv(summary_file, index=False)
        logger.info(f"Saved gene summary to {summary_file}")
        
        # Print statistics
        logger.info("\n" + "="*60)
        logger.info("ANALYSIS COMPLETE")
        logger.info("="*60)
        logger.info(f"Total LoF variants: {len(df):,}")
        logger.info(f"Total genes with LoF: {df['gene_name'].nunique():,}")
        logger.info(f"Average LoF per gene: {len(df) / df['gene_name'].nunique():.1f}")
        
        # Mutation type breakdown
        logger.info("\nMutation type breakdown:")
        for mut_type, count in df['mutation_type'].value_counts().head(10).items():
            logger.info(f"  {mut_type}: {count:,} ({count/len(df)*100:.1f}%)")
        
        # Top genes by LoF count
        logger.info("\nTop 10 genes by LoF variant count:")
        top_genes = df['gene_name'].value_counts().head(10)
        for gene, count in top_genes.items():
            logger.info(f"  {gene}: {count:,} variants")
    
    else:
        logger.error("No data processed. Check logs for errors.")
    
    logger.info("\nPipeline complete!")

if __name__ == "__main__":
    main()
