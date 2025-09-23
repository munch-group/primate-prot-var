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
        if entry['dbname'] == 'Uniprot_gn' and entry['primary_id'].startswith('A'):
            return entry['primary_id']


subprocess.check_call("wget https://zenodo.org/records/8208688/files/AlphaMissense_aa_substitutions.tsv.gz".split())
alpha_missense_file = 'alpha_missense_hg38.h5'

df = pd.read_csv('AlphaMissense_aa_substitutions.tsv.gz', sep="\t", comment='#')

store = pd.HDFStore(alpha_missense_file, 'a')  
keys = set(store.keys())
groups = df.groupby('uniprot_id')
for name, group in tqdm(groups):
    print(name)
    if name in keys:
        continue
    try:
        store.put(name, group)  
    except Exception as e:
        print(f"Error storing {name}: {e}")
store.get(name)  
print(len(store.keys()))
store.close()  


#pd.read_hdf('alpha_missense_hg38.h5', 'df' ,where=f'uniprot_id>{uniprot_id}')

sys.exit()



taxid = {
    'primates':    9443,
    'simiiformes': 314293,  # monkeys / apes
    'catarrhini':  9526,    # old world monkeys / apes
    'hominoidea':  314295,  # all apes
    'hominidae':   9604,    # great apes
}

target_taxon, gene_symbol = 'catarrhini',  'TTLL10'

ensembl_id = hgnc_to_ensembl(gene_symbol)
uniprot_id = ensembl_to_uniprot(ensembl_id)

server = "https://rest.ensembl.org"

ext = f"/homology/id/human/{ensembl_id}?compara=vertebrates;type=orthologues;sequence=cdna;aligned=0;target_taxon={taxid[target_taxon]}"
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
if not r.ok:
    r.raise_for_status()
    sys.exit()
decoded = r.json()
orthologs_nt_seq = f"{gene_symbol}_nt_seq.fasta"
with open(orthologs_nt_file, "w") as f:
    for entry in decoded['data']:
        for homology in entry['homologies']:
            print(f">{homology['target']['species']}\n{homology['target']['seq']}", file=f)
        break

ext = f"/homology/id/human/{ensembl_id}?compara=vertebrates;type=orthologues;sequence=protein;aligned=1;target_taxon={taxid[target_taxon]}"
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


temp_dir = tempfile.mkdtemp(prefix="pre_",suffix="_suf")

tmp_out_aa_aln = f"{temp_dir}/aligned_AA.fa"
tmp_out_nt_aln = f"{temp_dir}/aligned_NT.fa"

out_aa_aln = f"{gene_symbol}_{ensembl_id}_{uniprot_id}_aligned_AA.fa"
out_nt_aln = f"{gene_symbol}_{ensembl_id}_{uniprot_id}_aligned_NT.fa"

# align
cmd = f"macse -prog alignSequences -seq {orthologs_nt_seq} -out_NT {tmp_out_nt_aln} -out_AA {tmp_out_aa_aln} -gc_def 1 -local_realign_init 1 -local_realign_dec 1"
print(cmd)
subprocess.run(cmd.split())  # Be polite to the server

# export

cmd = f"macse -prog exportAlignment -align {tmp_out_nt_aln} -codonForInternalStop NNN -codonForInternalFS --- -charForRemainingFS --- -out {out_nt_aln} -out_AA {out_aa_aln}"
print(cmd)
subprocess.run(cmd.split())



#subprocess.run(f"macse -prog reportGapsAA2NT -align_AA {orthologs_aa_file} -seq {orthologs_nt_file} -out_NT output_NT.fasta".split())




# if not os.path.exists(alpha_missense_file):
#     cmds = f"""
#     wget https://zenodo.org/records/8208688/files/AlphaMissense_aa_substitutions.tsv.gz?download=1
#     mv AlphaMissense_aa_substitutions.tsv.gz?download=1 AlphaMissense_aa_substitutions.tsv.gz
#     gzip -d --stdout AlphaMissense_aa_substitutions.tsv.gz | grep -v '#' > alpha_missense_hg38.tsv
#     python -c 'import pandas ; pandas.read_csv("alpha_missense_hg38.tsv", sep="\t").to_hdf("{alpha_missense_file}", key="df", format="table", data_columns=["uniprot_id"]);'
#     """.split('\n')
#     for cmd in cmds:
#         if cmd.strip():
#             subprocess.check_call(cmd.split())

