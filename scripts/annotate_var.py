
import pandas as pd
import sys
import numpy as np

_, gene_id, alignment_file_name, hdf_file_name = sys.argv
store = pd.HDFStore(hdf_file_name, 'r')  
key = f'{gene_id[0].upper()}/{gene_id}'
df = store.get(key)
# df['ref'] = [chr(x) for x in df['ref']]
# df['alt'] = [chr(x) for x in df['alt']]

print(df.reset_index(drop=True).head())

