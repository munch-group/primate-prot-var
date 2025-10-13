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

def get_cdna_seqs(ensembl_id, taxids):
    server = "https://rest.ensembl.org"
    ext = f"/homology/id/human/{ensembl_id}?compara=vertebrates;type=orthologues;sequence=cdna;aligned=0;target_taxon={taxid};target_taxon={'simiiformes'}"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    orthologs_nt_seq = f"{gene_symbol}_nt_seq.fasta"
    with open(orthologs_nt_seq, "w") as f:
        for entry in decoded['data']:
            for homology in entry['homologies']:
                print(f">{homology['target']['species']}\n{homology['target']['seq']}", file=f)
            break
        else:
            raise ValueError(f"No homologies found for {gene_symbol}")
    return orthologs_nt_seq

def get_aa_alignment(ensembl_id, taxids):
    taxon_options = ';'.join(['target_taxon={t}' for t in taxids])
    server = "https://rest.ensembl.org"
    ext = f"/homology/id/human/{ensembl_id}?compara=vertebrates;type=orthologues;sequence=protein;aligned=1;target_taxon={taxid};target_taxon={'simiiformes'}"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    orthologs_aa_aln = f"{gene_symbol}_aa_align.fasta"
    with open(orthologs_aa_aln, "w") as f:
        for entry in decoded['data']:
            for homology in entry['homologies']:                
                print(f">{homology['target']['species']}\n{homology['target']['align_seq']}", file=f)
            break
        else:
            raise ValueError(f"No homologies found for {gene_symbol}")
    return orthologs_aa_aln

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

    ext = f"/taxonomy/name/{taxon_name}?"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    decoded = r.json()
    for entry in decoded:
        if entry['scientific_name'].lower() == taxon_name.lower():
            return entry['id']
    else:
        raise ValueError(f"Taxon {taxon_name} not found")


import argparse

parser = argparse.ArgumentParser(description="Get orthologs for a gene")
parser.add_argument('-t', "--taxon", dest='target_taxon', action='append', help="Target taxon")
parser.add_argument("gene_symbol", help="Gene symbol")
args = parser.parse_args()


# primates:     
# simiiformes:     all monkeys / apes
# catarrhini:      old world monkeys / apes
# cercopithecidae: old world monkeys
# hominoidea:      all apes
# hominidae:       great apes


#target_taxon, gene_symbol = 'catarrhini',  'TTLL10'
#_, target_taxon, gene_symbol = sys.argv

ensembl_id = hgnc_to_ensembl(args.gene_symbol)
uniprot_id = ensembl_to_uniprot(ensembl_id)

target_taxon_ids = [get_taxid(taxon) for taxon in args.target_taxon]

orthologs_nt_seq = get_orthologs(ensembl_id, target_taxon_ids, sequence="cdna", aligned=False)
orthologs_aa_aln = get_orthologs(ensembl_id, target_taxon_ids, sequence="protein", aligned=True)

# orthologs_nt_seq = get_cdna_seqs(ensembl_id, target_taxon_ids)
# orthologs_aa_aln = get_aa_alignment(ensembl_id, target_taxon_ids)

temp_dir = tempfile.mkdtemp(prefix="pre_",suffix="_suf")

tmp_out_aa_aln = f"{temp_dir}/aligned_AA.fa"
tmp_out_nt_aln = f"{temp_dir}/aligned_NT.fa"

out_aa_aln = f"{args.gene_symbol}_{ensembl_id}_{uniprot_id}_{'-'.join(args.target_taxon)}_aln_aa.fa"
out_nt_aln = f"{args.gene_symbol}_{ensembl_id}_{uniprot_id}_{'-'.join(args.target_taxon)}_cds_seqs.fa"

cmd = f"macse -prog alignSequences -seq {orthologs_nt_seq} -out_NT {tmp_out_nt_aln} -out_AA {tmp_out_aa_aln} -gc_def 1 -local_realign_init 1 -local_realign_dec 1"
print(cmd)
subprocess.run(cmd.split())  # Be polite to the server

cmd = f"macse -prog exportAlignment -align {tmp_out_nt_aln} -codonForInternalStop NNN -codonForInternalFS --- -charForRemainingFS --- -out_NT {out_nt_aln} -out_AA {out_aa_aln}"
print(cmd)
subprocess.run(cmd.split())


