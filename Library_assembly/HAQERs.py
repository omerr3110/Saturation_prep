import pandas as pd
import datatable as dt
import os
# import xlrd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from functions import seq_retrieve
from functions import df_tiling
#%%
window_size=135
filter_threshold=750
# %%
pos_df = pd.read_excel(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly"
                       r"\HAQERs\HAQERs_supp_data.xlsx", sheet_name="HaqerBed")
starrseq_df = pd.read_excel(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly"
                            r"\HAQERs\HAQERs_supp_data.xlsx", sheet_name="scSTARR-seq Tested Sequences")
# %%
pos_df['START']=pos_df['START']-1 #switch to 0-format
pos_df.to_csv("to_lift_haqers.bed", header=False, sep="\t", index=False)
# %%
lifted_pos_df = pd.read_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly"
                            r"\HAQERs\lifted_haqers.bed", sep="\t", names=pos_df.columns.to_list() + ['score'])
lifted_pos_df['START']=lifted_pos_df['START']+1 #switch back to 1-based
pos_df['START']=pos_df['START']+1
# %%
bool_vec = lifted_pos_df.groupby("NAME")['score'].max() > 1
lifted_pos_df['final'] = lifted_pos_df.apply(func=lambda x: x['NAME'] in bool_vec[~bool_vec], axis=1)
# %%
final_pos_df = lifted_pos_df[lifted_pos_df['final']].copy()
# %%
comb_df = pd.merge(final_pos_df, pos_df, on="NAME")
comb_df['length_19'] = comb_df.apply(func=lambda x: (x[2] - x[1]),axis=1)
comb_df['length_38'] = comb_df.apply(func=lambda x: (x[8] - x[7]),axis=1)
comb_df['diff'] = comb_df.apply(func=lambda x: (x['length_19'] - x['length_38']),axis=1)
#%%
filtered_comb_df=comb_df[(comb_df['length_19']<=filter_threshold) & (comb_df['diff']==0)]
#%%
sns.histplot(data=comb_df[comb_df['length_19']<=1500],x='length_19').set(title="HAQERs Length")
plt.show()
# %%

haqer_tiles_df=df_tiling(filtered_comb_df,window_size,270)
print(haqer_tiles_df.shape)

#%%
final_df_haqers=pd.DataFrame(data={"chr":haqer_tiles_df['chromosome'], "start":haqer_tiles_df['oligo_start'],
                                   "end":haqer_tiles_df['oligo_end'],
                                   "group":"HAQERs","name":haqer_tiles_df['name'],"external_ID":haqer_tiles_df['ID']})
#%%
genome_file="hg19.fa"
records = list(SeqIO.parse(genome_file, "fasta"))
#%%
final_df_haqers['sequence_hg19']=final_df_haqers.apply(func=lambda x: seq_retrieve(x['chr'],x['start'],x['end'],records),axis=1)
#%%
final_df_haqers.to_csv(fr"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\output\pre_satmut_haqers_{filter_threshold}_{window_size}.csv")