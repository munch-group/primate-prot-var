import requests
import subprocess
import pandas as pd
import sys
import tempfile
import os
from tqdm import tqdm
 

def hgnc_to_ensembl(hgnc_symbol):
    """Convert HGNC symbol to Ensembl ID"""
    url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{hgnc_symbol}"
    response = requests.get(url, headers={"Content-Type": "application/json"})
    if response.status_code == 200:
        data = response.json()
        for entry in data:
            if entry['type'] == 'gene':
                return entry['id']

def ensembl_to_uniprot(ensembl_id):
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/id/{ensembl_id}?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()    
    decoded = r.json()
    for entry in decoded:
        if entry['dbname'] == 'Uniprot_gn':
            return entry['primary_id']

def get_human_cds_seqs(ensembl_id, sequence):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{ensembl_id}?type={sequence};multiple_sequences=1"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    records = {entry['id']: entry['seq'] for entry in decoded}
    return records

def get_orthologs(ensembl_id, taxids, sequence, aligned=False):
    server = "https://rest.ensembl.org"
    ext = f"/homology/id/human/{ensembl_id}?compara=vertebrates;type=orthologues;sequence={sequence};aligned={int(aligned)};{';'.join([f'target_taxon={t}' for t in taxids])}"
    print(ext)
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    orthologs_file = f"{ensembl_id}_{sequence}_{'aln' if aligned else 'seq'}.fa"
    with open(orthologs_file, "w") as f:
        for entry in decoded['data']:
            seq_key = 'align_seq' if aligned else 'seq'
            for homology in entry['homologies']:
                print(f">{homology['target']['species']}\n{homology['target'][seq_key]}", file=f)
            break
        else:
            raise ValueError(f"No homologies found for {ensembl_id}")
    return orthologs_file

def get_taxid(taxon_name):
    server = "https://rest.ensembl.org"
    ext = f"/taxonomy/id/{taxon_name}?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    return decoded['id']

import argparse

taxa = [
    'primates',        # all primates    
    'simiiformes',     # all monkeys / apes
    'catarrhini',      # old world monkeys / apes
    'cercopithecidae', # old world monkeys (african/asian with tails)
    'hominoidea',      # all apes (great apes and gibbons)
    'hominidae',       # great apes (chimpanzee, bonobo, gorilla, orangutan)
]
parser = argparse.ArgumentParser(description="Get orthologs for a gene")
parser.add_argument('-t', "--taxon", choices=taxa, dest='target_taxon', action='append', help="Target taxon")
parser.add_argument("gene_symbol", help="Gene symbol")
args = parser.parse_args()

# id conversions
ensembl_id = hgnc_to_ensembl(args.gene_symbol)
uniprot_id = ensembl_to_uniprot(ensembl_id)

# translate taxon names to taxon ids
target_taxon_ids = [get_taxid(taxon) for taxon in args.target_taxon]

# get orthologs cds sequences
orthologs_nt_seq = get_orthologs(ensembl_id, target_taxon_ids, sequence="cdna", aligned=False)

# add human cds variants
with open(orthologs_nt_seq, "a") as f:
    for name, cds in get_human_cds_seqs(ensembl_id, "cdna").items():
        f.write(f">{name}\n{cds}\n")

# file prefix
prefix = f'{args.gene_symbol}_{ensembl_id}_{uniprot_id}_{'-'.join([x.replace(' ', '_') for x in args.target_taxon])}'

# align cds with macse
temp_dir = tempfile.mkdtemp(prefix="pre_",suffix="_suf")
tmp_out_aa_aln = f"{temp_dir}/aligned_AA.fa"
tmp_out_nt_aln = f"{temp_dir}/aligned_NT.fa"
cmd = f"macse -prog alignSequences -seq {orthologs_nt_seq} -out_NT {tmp_out_nt_aln} -out_AA {tmp_out_aa_aln} -gc_def 1 -local_realign_init 1 -local_realign_dec 1"
print(cmd)
subprocess.run(cmd.split())

# enrich alignment with extra sequences (if any)
# extra_cds_seqs = 'results/TTLL10.fa'
# cmd = f'cmd = f"macse -prog enrichAlignment -align {out_nt_aln} -seq {extra_cds_seqs}'
# print(cmd)
# subprocess.run(cmd.split())

# export stats and final alignments                                             
out_aa_aln = f"{prefix}_protein.fa"
out_nt_aln = f"{prefix}_cds.fa"
out_stats_aln = f"{prefix}_stats.csv"
cmd = f"macse -prog exportAlignment -align {tmp_out_nt_aln} -out_NT {out_nt_aln} -out_AA {out_aa_aln} -out_stat_per_seq {out_stats_aln}"
# cmd = f"macse -prog exportAlignment -align {tmp_out_nt_aln} -codonForInternalStop NNN -codonForInternalFS --- -charForRemainingFS --- -out_NT {out_nt_aln} -out_AA {out_aa_aln} -out_stat_per_seq VERSION2_{out_stats_aln}"
print(cmd)
subprocess.run(cmd.split())


