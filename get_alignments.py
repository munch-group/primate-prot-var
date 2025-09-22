from tqdm import tqdm
import requests, sys
import subprocess


def hgnc_to_ensembl(hgnc_symbol):
    """Convert HGNC symbol to Ensembl ID"""
    url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{hgnc_symbol}"
    response = requests.get(url, headers={"Content-Type": "application/json"})
    
    if response.status_code == 200:
        data = response.json()
        for entry in data:
            if entry['type'] == 'gene':
                return entry['id']
    return None

# Example usage
for gene_symbol in ["TTLL10"]:
    ensembl_id = hgnc_to_ensembl(gene_symbol)

    print(gene_symbol)

    server = "https://rest.ensembl.org"
    # ext = f"/homology/id/human/{ensembl_id}?compara=vertebrates;type=orthologues;sequence=cdna;aligned=1;target_taxon=9443"
    ext = f"/homology/id/human/{ensembl_id}?compara=vertebrates;type=orthologues;sequence=cdna;aligned=0;target_taxon=9443"

    #ext = "/homology/id/human/ENSG00000157764?sequence=cdna;target_taxon=9443;type=orthologues"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    
    decoded = r.json()
    orthologs_file = f"{gene_symbol}_orthologs.fasta"
    with open(orthologs_file, "w") as f:
        for entry in decoded['data']:
            for homology in entry['homologies']:
                print(f">{homology['target']['species']}\n{homology['target']['seq']}", file=f)
            break

    subprocess.run(f"java -jar macse_v2.07.jar -prog alignSequences -seq {orthologs_file} -out_NT {gene_symbol}_aligned_NT.fasta -out_AA {gene_symbol}_aligned_AA.fasta -gc_def 1".split())  # Be polite to the server


    