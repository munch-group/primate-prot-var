import subprocess
from tqdm import tqdm
import pandas as pd

subprocess.check_call("wget https://zenodo.org/records/8208688/files/AlphaMissense_aa_substitutions.tsv.gz".split())
alpha_missense_file = 'alpha_missense_hg38.h5'

df = pd.read_csv('AlphaMissense_aa_substitutions.tsv.gz', sep="\t", comment='#')
store = pd.HDFStore('alpha_missense_hg38.h5', 'a')  
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
