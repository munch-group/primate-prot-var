import subprocess
from tqdm import tqdm
import pandas as pd
import numpy as np

#subprocess.check_call("wget https://zenodo.org/records/8208688/files/AlphaMissense_aa_substitutions.tsv.gz".split())
alpha_missense_file = 'alpha_missense_hg38.h5'

df = pd.read_csv('AlphaMissense_aa_substitutions.tsv.gz', sep="\t", comment='#')
store = pd.HDFStore('alpha_missense_hg38.h5', 'w')  
keys = set(store.keys())
groups = df.groupby('uniprot_id')
for name, group in tqdm(groups):
    gr = group.copy()
    # gr['ref'] = [ord(x[0]) for x in group.protein_variant]
    # gr['pos'] = [int(x[1:-1]) for x in group.protein_variant]
    # gr['alt'] = [ord(x[-1]) for x in group.protein_variant]
    # gr['ref'] = gr['ref'].astype(np.uint8)
    # gr['alt'] = gr['alt'].astype(np.uint8)
    # gr['pos'] = gr['pos'].astype(np.uint16)
    # gr.drop(columns=['protein_variant'], inplace=True)

    if name in keys:
        continue
    try:
        store.put(f'{name[0].upper()}/{name}', gr, index=False)  
    except Exception as e:
#        print(f"Error storing {name}: {e}")
        raise e
store.close()  
