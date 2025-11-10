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
#%%
variants_df = pd.read_excel(r"C:\Users\omerro\Dropbox (Weizmann Institute)\Gokhman lab general info"
                             r"\USEFUL DATASETS\Expression\MPRAs\Modern human derived MPRA\elife-63713-supp1-v1.xlsx",
                             header=2, sheet_name=1)
#%%
variants_df['AH_variant']=variants_df.apply(lambda x: x['Archaic sequence sequence'][99],axis=1)
variants_df['MH_variant']=variants_df.apply(lambda x: x['Modern sequence sequence'][99],axis=1)

#%%
var_dict1=dict(zip(tuple(zip(variants_df.Chromosome,variants_df['Central variant position (hg19)']))
                  ,variants_df.AH_variant))
var_dict2=dict(zip(tuple(zip(variants_df.Chromosome,variants_df['Central variant position (hg19)']))
                  ,variants_df.MH_variant))
#%%
with open(r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\Tomas_Apes_Catalog_parsing_5"
                        r"\Human_ancestral_dict.pkl", 'wb') as f:
    pickle.dump(var_dict1, f)
with open(r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\Tomas_Apes_Catalog_parsing_5"
          r"\MH_derived_dict.pkl", 'wb') as f:
    pickle.dump(var_dict2, f)


#%%
with open(r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\Tomas_Apes_Catalog_parsing_5"
          r"\Human_ancestral_dict.pkl", 'rb') as f:
    variant_dictionary = pickle.load(f)


#%%
def H_ancestral_modify(row, d):
    seq = row['sequence_hg19']
    chrom = row['chr']
    start = row['start']
    end = row['end']
    new_d = {k: v for k, v in d.items() if k[0] == chrom and end >= k[1] >= start}
    print(len(new_d))
    if len(new_d) == 0:
        return seq
    H_anc_seq_list = list(seq)
    for k in new_d.keys():
        var_nuc = new_d[k]
        oligo_pos = k[1] - start
        H_anc_seq_list[oligo_pos] = var_nuc
    H_anc_seq = "".join(H_anc_seq_list)
    return H_anc_seq



#%%
chosen_oligos=variants_df.iloc[0:2].copy()
#%%
test_df=pd.DataFrame(data={"chr":chosen_oligos['Chromosome'],"start":chosen_oligos['Oligo start (hg19)'],
                           "end":chosen_oligos['Oligo end (hg19)'],"sequence_hg19":chosen_oligos['Modern sequence sequence']})
#%%
anc_seq=test_df.apply(H_ancestral_modify,axis=1,args=(var_dict1,))
