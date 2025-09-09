import sys, os
import sgkit as sg
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import gffutils
import pandas as pd


def parse_gtf(gtf_file_name):

    # Parse GTF file (works with both GTF and GFF3 formats)
    db = gffutils.create_db(
        gtf_file_name, 
        dbfn=":memory:", 
        force=True,     
        keep_order=True,
        merge_strategy="merge", 
        sort_attribute_values=True, 
        disable_infer_transcripts=True,
        )

    genes_processed = {}

    for transcript in db.features_of_type('transcript'):
        # Get gene info
        gene_name = transcript.attributes.get('gene_name', [''])[0]
        gene_id = transcript.attributes.get('gene_id', [''])[0]
        if not gene_name and not gene_id:
            continue
            
        gene_key = gene_name if gene_name else gene_id
        
        # Get CDS features for this transcript
        cds_features = list(db.children(transcript, featuretype='CDS'))
        
        if not cds_features:
            continue
        
        # Sort CDS features by start position
        cds_features.sort(key=lambda x: x.start)
        cds_coords = [(cds.start-1, cds.end) for cds in cds_features]
        cds_length = sum([e-s for s, e in cds_coords])

        genes_processed = {}
        # get only longest (canonical) transcript
        if gene_key not in genes_processed or (cds_length > genes_processed[gene_key]['length'] and cds_length < 6):

            genes_processed[gene_key] = {
                'transcript_id': transcript.attributes.get('transcript_id', [transcript.id])[0],
                'chromosome': transcript.chrom,
                'gene_name': gene_name,
                'cds_coords': cds_coords,
                'strand': transcript.strand,
            }

    return genes_processed


def get_cds_haplotypes(genome, chrom, cds_coords, strand, region_ds=None, sample_idx=None):

    chrom_seq = genome[chrom]

    # make the reference cds looks ok:
    ref_cds = ''.join([str(chrom_seq[s:e]) for s, e in cds_coords])
    s, e = zip(*cds_coords)
    cds_start, cds_end = min(s), max(e)
    if strand == '-':
        ref_cds = chrom_seq[cds_start-3:cds_start] + ref_cds # include stop codon at beginning
        ref_cds = str(Seq(ref_cds).reverse_complement().seq)
    else:
        ref_cds = ref_cds + chrom_seq[cds_end:cds_end+3] # include stop codon at end
    assert ref_cds[:3] == 'ATG'
    assert len(ref_cds) % 3 == 0
    stop_codons = set(['TAA', 'TAG', 'TGA'])
    assert ref_cds[-3:] in stop_codons
    ref_cds = ref_cds[:-3] # remove stop codon again    
    codons = set([ref_cds[i:i+3] for i in range(0, len(ref_cds)-3+1, 3)])
    assert not codons.intersection(stop_codons)

    if region_ds is None and sample_idx is None:
        return ref_cds

    # make reference haplotypes
    hap1 = list(chrom_seq)
    hap2 = list(chrom_seq)

    cds_haps = []
    # add variants to haplotypes
    for i in range(len(region_ds.variants)):
        alleles = region_ds.variant_allele.values[i]
        ref = alleles[0]
        alts = [a for a in alleles[1:] if a != '' and a != '.']

        # Get genotype for this sample at this variant if available
        if 'call_genotype' in region_ds:
            genotypes = region_ds.call_genotype.values[i, sample_idx]  # Get genotype for variant i, sample sample_idx
            
            pos = int(region_ds.variant_position.values[i])
            
            if genotypes[0] >= 0 and genotypes[0] < len(alleles):
                allele1 = alleles[genotypes[0]]
                if allele1 and allele1 != '.':
                    if len(ref) == 1 and len(allele1) == 1:
                        hap1[pos] = allele1
                    # TODO: Handle indels properly
                    
            if genotypes[1] >= 0 and genotypes[1] < len(alleles):
                allele2 = alleles[genotypes[1]]
                if allele2 and allele2 != '.':
                    if len(ref) == 1 and len(allele2) == 1:
                        hap2[pos] = allele2
                    # TODO: Handle indels properly

    cds1_parts = []
    cds2_parts = []
    for s, e in cds_coords:
        cds1_parts.append(''.join(hap1[s:e]))
        cds2_parts.append(''.join(hap2[s:e]))
    
    return ''.join(cds1_parts), ''.join(cds2_parts)

chromosomes = {
    'NC_044976.1': 'chr1',
    'NC_044977.1': 'chr2',
    'NC_044978.1': 'chr3',
    'NC_044979.1': 'chr4',
    'NC_044980.1': 'chr5',
    'NC_044981.1': 'chr6',
    'NC_044982.1': 'chr7',
    'NC_044983.1': 'chr8',
    'NC_044984.1': 'chr9',
    'NC_044985.1': 'chr10',
    'NC_044986.1': 'chr11',
    'NC_044987.1': 'chr12',
    'NC_044988.1': 'chr13',
    'NC_044989.1': 'chr14',
    'NC_044990.1': 'chr15',
    'NC_044991.1': 'chr16',
    'NC_044992.1': 'chr17',
    'NC_044993.1': 'chr18',
    'NC_044994.1': 'chr19',
    'NC_044995.1': 'chr20',
    'NC_044996.1': 'chrX',
    'NC_044997.1': 'chrY',
}

gene_list = ['TTLL10', 'DYNLT3', 'CFAP47', 'DIAPH2', 'PRICKLE3']

# Load genome sequence
genome = SeqIO.to_dict(SeqIO.parse("papAnu4.fa", "fasta"))
for key in genome:
    genome[key] = str(genome[key].seq)

species_list = ['Papio_hamadryas', 'Papio_anubis', 'Papio_papio', 
                'Papio_kindae', 'Papio_cynocephalus', 'Papio_ursinus']

transcript_info = parse_gtf('papAnu4.ncbiRefSeq.gtf')
#transcript_info = parse_gtf('TTLL10.gtf')

records = []



for species in species_list:

    for chrom_id in chromosomes:

        zarr_path = os.path.expanduser(f"~/primatediversity/people/erik/Diversity_Selection_Primates/zarr_data/{species}_ssp/{chrom_id}")

        if not os.path.exists(zarr_path):
            print(f"Zarr path {zarr_path} does not exist, skipping...", file=sys.stderr)
            continue

        ds = sg.load_dataset(zarr_path)

        for gene_name in transcript_info:
            print(gene_name)
            if transcript_info[gene_name]['chromosome'] != chromosomes[chrom_id]:
                continue
            if gene_list and gene_name not in gene_list:
                continue

            transcript = transcript_info[gene_name]
            cds_coords = transcript['cds_coords']
            
            ref_cds = get_cds_haplotypes(genome, chromosomes[chrom_id], cds_coords, transcript['strand'])
            prot_ref = str(Seq(ref_cds).translate())

            start_region = min(start for start, end in cds_coords)
            end_region = max(end for start, end in cds_coords)
            position_mask = (ds.variant_position >= start_region) & (ds.variant_position <= end_region)
            
            # Compute mask if it's a dask array
            if hasattr(position_mask, 'compute'):
                position_mask_computed = position_mask.compute()
            else:
                position_mask_computed = position_mask
            
            # Filter dataset to gene region
            region_ds = ds.sel(variants=position_mask_computed)
            
            # iterate over samples 
            for j, sample_id in enumerate(region_ds.sample_id.values):            
                sample_idx = np.where(region_ds.sample_id.values == sample_id)[0][0]

                # Extract CDS sequences for both haplotypes of this sample
                for h, cds_hap in enumerate(get_cds_haplotypes(genome, chromosomes[chrom_id], cds_coords, transcript['strand'], region_ds, sample_idx)):
                    prot_hap = str(Seq(cds_hap).translate())
                    inframe_stop = None
                    nsyn = []
                    for i, (ref, hap) in enumerate(zip(prot_ref, prot_hap)):
                        if ref != hap:
                            nsyn.append(f'{ref}{i}{hap}')
                            if hap == '*':
                                inframe_stop = i
                    for s in nsyn:
                        records.append([transcript['transcript_id'], gene_name, sample_id, h+1, len(prot_hap), s, int(inframe_stop)])
                    nsyn_str = ','.join(nsyn) if nsyn else None
                    if inframe_stop is not None:
                        prot_hap = prot_hap[:inframe_stop] + prot_hap[inframe_stop:].lower()
                    print(f">{transcript['transcript_id']} gene={gene_name} sample={sample_id} haplotype={h+1} species={species} length={len(prot_hap)} IFS={inframe_stop} NSYN={nsyn_str}")
                    print(prot_hap)


pd.DataFrame.from_records(records, columns=['transcript_id', 'gene_name', 'sample_id', 'haplotype', 'protein_length', 'change', 'inframe_stop']).to_csv('prot_variants.csv', index=False)

# Number of variants discovered in this study                 
# RefSeqID    Chromosome  Human Ortholog (23 pairs)   Length (bp) - assembly v1.0 SNVs identified Indels identified
# NC_044976.1 Chr1    1   218,172,882 1,253,765   280,563
# NC_044977.1 Chr2    3   193,660,750 1,106,523   234,198
# NC_044978.1 Chr3    7/21    184,919,515 1,185,115   240,921
# NC_044979.1 Chr4    6   182,120,902 1,103,802   244,805
# NC_044980.1 Chr5    4   173,900,761 1,040,353   220,224
# NC_044981.1 Chr6    5   167,138,247 1,001,262   212,306
# NC_044982.1 Chr7    14/15   161,768,468 887,977 204,077
# NC_044983.1 Chr8    8   140,274,886 864,665 176,801
# NC_044984.1 Chr9    10  127,591,819 779,449 174,978
# NC_044985.1 Chr10   20/22   126,462,689 750,767 157,493
# NC_044986.1 Chr11   12  125,913,696 802,993 169,432
# NC_044987.1 Chr12   2q  123,343,450 715,469 155,406
# NC_044988.1 Chr13   2p  106,849,001 626,461 141,214
# NC_044989.1 Chr14   11  106,654,974 648,270 133,342
# NC_044990.1 Chr15   9   91,985,775  630,624 130,269
# NC_044991.1 Chr16   17  91,184,193  600,654 140,801
# NC_044992.1 Chr17   13  74,525,926  465,873 120,520
# NC_044993.1 Chr18   18  72,894,408  482,937 105,971
# NC_044994.1 Chr19   19  72,123,344  462,891 94,878
# NC_044995.1 Chr20   16  50,021,108  384,351 120,698
# NC_044996.1 ChrX        142,711,496 356,041 126,813
# NC_044997.1 ChrY        8,309,886   28,099  12,105
# Unplaced    11,122 scaffolds        127,276,471     
# total           2,869,804,647   17,033,371  3,838,917
