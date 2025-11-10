#%%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from functions import seq_retrieve
from functions import df_tiling
import regex as re
import Bio
import pickle
import sys

"""This section takes Ryder's catalog and turns it to a dictionary of AH variants"""
variants_df=pd.read_excel(r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\Tomas_Apes_Catalog_parsing_5"
                        r"\archaic_derived_MPRA_oligos SF.xlsx")
variants_df['AH_variant']=variants_df.apply(lambda x: x['archaic_seq'][134],axis=1)
var_dict=dict(zip(tuple(zip(variants_df.chr,variants_df.variant_pos)),variants_df.AH_variant))
# the keys in this dictionary use integers as chromosome id and not string (eg. "chr11")
with open(r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\Tomas_Apes_Catalog_parsing_5"
                        r"\archaic_derived_variants_dict.pkl", 'wb') as f:
    pickle.dump(var_dict, f)
#%%

#%%


with open(r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\Tomas_Apes_Catalog_parsing_5"
          r"\archaic_derived_variants_dict.pkl", 'rb') as f:
    variant_dictionary = pickle.load(f)
#%%
# here we create a dictionary with strings instead of integers
new_dict={("chr"+str(k[0]),k[1]): v for k,v in variant_dictionary.items()}
with open(r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\Tomas_Apes_Catalog_parsing_5"
                        r"\archaic_derived_variants_dict_str.pkl", 'wb') as f:
    pickle.dump(new_dict, f)
#%%
def archaic_modify(row, d):
    seq = row['sequence_hg19']
    chrom = re.findall("\d+|X|Y", row['chr'])[0]
    if chrom.isnumeric():
        chrom=int(chrom)
    start = row['start']
    end = row['end']
    new_d = {k: v for k, v in d.items() if k[0] == chrom and end >= k[1] >= start}
    print(len(new_d))
    if len(new_d) == 0:
        return seq
    archaic_seq_list = list(seq)
    for k in new_d.keys():
        var_nuc = new_d[k]
        oligo_pos = k[1] - start
        archaic_seq_list[oligo_pos] = var_nuc
    archaic_seq = "".join(archaic_seq_list)
    return archaic_seq


print(sys.argv)
path = sys.argv[1]
df = pd.read_csv(path, index_col=0)
df['sequence_AH'] = df.apply(archaic_modify, axis=1, args=(variant_dictionary,))

#%%
# QA
chosen_oligos=variants_df.iloc[0:2].copy()
#%%
test_df=pd.DataFrame(data={"chr":chosen_oligos['chr'].apply(lambda x: str(x)),"start":chosen_oligos['seq_start'],
                           "end":chosen_oligos['seq_end'],"sequence_hg19":chosen_oligos['modernhg19_seq']})
#%%
arc_seq=test_df.apply(archaic_modify,axis=1,args=(var_dict,))

