import pandas as pd
import datatable as dt
import os
# import xlrd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from functions import seq_retrieve
# %%
#carly_df_esc = pd.read_excel(r"C:\Users\omerro\Dropbox (Weizmann Institute)\Gokhman lab general info"
                             #r"\USEFUL DATASETS\Expression\MPRAs\Modern human derived MPRA\elife-63713-supp1-v1.xlsx",
                             #header=2, sheet_name=0)
carly_df_ost = pd.read_excel(r"C:\Users\omerro\Dropbox (Weizmann Institute)\Gokhman lab general info"
                             r"\USEFUL DATASETS\Expression\MPRAs\Modern human derived MPRA\elife-63713-supp1-v1.xlsx",
                             header=2, sheet_name=1)
carly_df_ost['ost']=True
carly_df_npc = pd.read_excel(r"C:\Users\omerro\Dropbox (Weizmann Institute)\Gokhman lab general info"
                             r"\USEFUL DATASETS\Expression\MPRAs\Modern human derived MPRA\elife-63713-supp1-v1.xlsx",
                             header=2, sheet_name=2)
carly_df_npc['npc']=True

# %%
carly_df_ost_active=carly_df_ost.query('`Active - osteoblast` == "Yes"')
carly_df_npc_active=carly_df_npc.query('`Active - NPC` == "Yes"')
carly_df_all_active=pd.concat([ carly_df_ost_active, carly_df_npc_active], axis=1, join="outer")

#%%
keep = [ 'Differentially active - Osteoblast', 'Active - osteoblast',
        'Differential activity log2(fold-change) - modern vs archaic - osteoblast',
        'Modern activity alpha - osteoblast', 'Differentially active - NPC', 'Active - NPC',
        'Differential activity log2(fold-change) - modern vs archaic - NPC',
        'Modern activity alpha - NPC','npc','ost']
#%%
carly_df_all_active = carly_df_all_active[keep].copy()
carly_df_all_active['Sequence ID'] = carly_df_ost.loc[carly_df_all_active.index.to_list(), 'Sequence ID']
carly_df_all_active['Chromosome'] = carly_df_ost.loc[carly_df_all_active.index.to_list(), 'Chromosome']
carly_df_all_active['Oligo start (hg19)'] = carly_df_ost.loc[carly_df_all_active.index.to_list(), 'Oligo start (hg19)']
carly_df_all_active['Oligo end (hg19)'] = carly_df_ost.loc[carly_df_all_active.index.to_list(), 'Oligo end (hg19)']
carly_df_all_active['Number of variants in sequence'] = carly_df_ost.loc[
    carly_df_all_active.index.to_list(), 'Number of variants in sequence']
carly_df_all_active.npc=carly_df_all_active.npc.fillna(False)
carly_df_all_active.ost=carly_df_all_active.ost.fillna(False)
carly_df_all_active['cell_type_id']=carly_df_all_active.apply(func= lambda x: 2*x[8]+x[9],axis=1)  #npc+ost=3, npc=2, ost=1
carly_df_all_active.drop(['npc','ost'],axis=1,inplace=True)


#%%

#%%
carly_df_all_active['Oligo start (hg19)'] = carly_df_all_active['Oligo start (hg19)'] - 35
carly_df_all_active['Oligo end (hg19)'] = carly_df_all_active['Oligo end (hg19)'] + 35
#%%

# %%
oligo_type = ["MH divergent SNV" for _ in range(carly_df_all_active.shape[0])]
ori = [None for _ in range(carly_df_all_active.shape[0])]
lib = ["1" for _ in range(carly_df_all_active.shape[0])]
seq_name = carly_df_all_active['Sequence ID'].apply(func=lambda x: x.split("_")[1]+"_tile1")
seq_id = carly_df_all_active['Sequence ID'].apply(func=lambda x: x.split("_")[1])
# %%
df_carly_pre_MPRA = pd.DataFrame(
    data={"chr": carly_df_all_active['Chromosome'], "start": carly_df_all_active['Oligo start (hg19)'],
          "end": carly_df_all_active['Oligo end (hg19)'] ,"group": oligo_type,"name": seq_name,'external_ID':seq_id
          })

# %%


#%%
genome_file="hg19.fa"
records = list(SeqIO.parse(genome_file, "fasta"))
#%%
df_carly_pre_MPRA['sequence_hg19']=df_carly_pre_MPRA.apply(func=lambda x: seq_retrieve(x['chr'],x['start'],x['end'],records),axis=1)
#%%
df_carly_pre_MPRA.to_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\output\pre_satmut_carly.csv")
# %%
#ryder_df = pd.read_excel(r"C:\Users\omerro\Dropbox (Weizmann Institute)\Gokhman lab general info\USEFUL DATASETS"
#                         r"\Expression\MPRAs\Archaic human derived MPRA\archaic_derived_MPRA_oligos.xlsx")
#filtered_ryder = ryder_df.query('number_of_variants>1')
# %%
single_df_ost = carly_df_ost.query('`Active - osteoblast` == "Yes" & `Number of variants in sequence`==1 & `Differentially active - Osteoblast` =="Yes"')
single_df_npc = carly_df_npc.query(' `Active - NPC` == "Yes" & `Number of variants in sequence`==1 & `Differentially active - NPC` =="Yes"')
single_df_all = pd.concat([ single_df_ost, single_df_npc], axis=1, join="outer")
single_df_all = single_df_all[keep].copy()
single_df_all['Sequence ID'] = carly_df_ost.loc[single_df_all.index.to_list(), 'Sequence ID']
single_df_all['Chromosome'] = carly_df_ost.loc[single_df_all.index.to_list(), 'Chromosome']
single_df_all['Oligo start (hg19)'] = carly_df_ost.loc[single_df_all.index.to_list(), 'Oligo start (hg19)']
single_df_all['Oligo end (hg19)'] = carly_df_ost.loc[single_df_all.index.to_list(), 'Oligo end (hg19)']
single_df_all['Number of variants in sequence'] = carly_df_ost.loc[
    single_df_all.index.to_list(), 'Number of variants in sequence']
single_df_all.npc=single_df_all.npc.fillna(False)
single_df_all.ost=single_df_all.ost.fillna(False)
single_df_all['cell_type_id']=single_df_all.apply(func= lambda x: 2*x[8]+x[9],axis=1)  #npc+ost=3, npc=2, ost=1
single_df_all.drop(['npc','ost'],axis=1,inplace=True)
#%%
filtered_df_ost = carly_df_ost.query('`Active - osteoblast` == "Yes" & `Number of variants in sequence`>1')
filtered_df_npc = carly_df_npc.query(' `Active - NPC` == "Yes" & `Number of variants in sequence`>1')
filtered_df_all = pd.concat([ filtered_df_ost, filtered_df_npc], axis=1, join="outer")
filtered_df_all = filtered_df_all[keep].copy()
filtered_df_all['Sequence ID'] = carly_df_ost.loc[filtered_df_all.index.to_list(), 'Sequence ID']
filtered_df_all['Chromosome'] = carly_df_ost.loc[filtered_df_all.index.to_list(), 'Chromosome']
filtered_df_all['Oligo start (hg19)'] = carly_df_ost.loc[filtered_df_all.index.to_list(), 'Oligo start (hg19)']
filtered_df_all['Oligo end (hg19)'] = carly_df_ost.loc[filtered_df_all.index.to_list(), 'Oligo end (hg19)']
filtered_df_all['Number of variants in sequence'] = carly_df_ost.loc[
    filtered_df_all.index.to_list(), 'Number of variants in sequence']
filtered_df_all.npc=filtered_df_all.npc.fillna(False)
filtered_df_all.ost=filtered_df_all.ost.fillna(False)
filtered_df_all['cell_type_id']=filtered_df_all.apply(func= lambda x: 2*x[8]+x[9],axis=1)  #npc+ost=3, npc=2, ost=1
filtered_df_all.drop(['npc','ost'],axis=1,inplace=True)
# %%