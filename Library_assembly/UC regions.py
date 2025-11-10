# uc.482 chrX:137876095-137876389
# uc.467 chrX:24369780-24370510
# uc.329 chr11:32162301-32162607
# uc.248	chr9:964189-964410
# %%
import pandas as pd
from bs4 import BeautifulSoup
from functions import df_tiling
import regex as re
import numpy as np
from functions import seq_retrieve
from Bio import SeqIO
from functions import anc_retreive


# %%

with open(r'/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/UC/UC_supp_2004.html') as file:
    soup = BeautifulSoup(file, 'html.parser')
tables = pd.read_html(str(soup))
uc_table = tables[1]
uc_table_final = uc_table[uc_table[('ultra conserved element', 'name')].str.contains('uc.')]
uc_table_final.to_csv(r'/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/UC/UC_table.csv')
# %%
uc_df = pd.read_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\UC"
                    r"\UC_table.csv", header=[0, 1], index_col=0)
pos = uc_df[('ultra conserved element', 'Build 34 (hg16) coords')].apply(lambda x: re.split(r':|-', x))
uc_df_to_lift = pd.DataFrame(
    data={"chr": pos.apply(func=lambda x: x[0]), "start": pos.apply(func=lambda x: int(x[1]) - 1),
          "end": pos.apply(func=lambda x: int(x[2]))})
uc_df_to_lift.to_csv("to_list_uc.bed", header=None, sep='\t', index=None)
# %%
uc_df_lifted = pd.read_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\UC"
                           r"\lifted_uc.bed", sep="\t", names=['chr', 'start', 'end', 'prev', 'num'],
                           dtype={'start': np.int32,
                                  'end': np.int32})
# %%
uc_df_filter=uc_df[[('ultra conserved element', 'name'),('ultra conserved element', 'type')]].reset_index()
#%%
uc_df_final = pd.concat(
    [uc_df_lifted[['chr', 'start', 'end']], uc_df_filter.iloc[:,1:3]].copy(), axis=1)
#%%
uc_df_final.columns = ['chr', 'start', 'end', 'name','type']
# %%
uc_df_final=uc_df_final[uc_df_final['type']=='n']


# %%
uc_df_tiles = df_tiling(uc_df_final, 135, 270)
uc_df_tiles['group']='Ultra Conserved'
#%%
genome_file="hg19.fa"
records = list(SeqIO.parse(genome_file, "fasta"))

#%%
uc_df_tiles['sequence_hg19']=uc_df_tiles.apply(func=lambda x: seq_retrieve(x['chromosome'],x['oligo_start'],
                                                                           x['oligo_end'],records),axis=1)
#%%
final_df_uc=pd.DataFrame(data={"chr":uc_df_tiles['chromosome'],"start":uc_df_tiles['oligo_start'],"end":uc_df_tiles['oligo_end']
                                ,"group":uc_df_tiles['group'],"name":uc_df_tiles['name'],
                               "external_ID":uc_df_tiles['ID'],'sequence_hg19':uc_df_tiles['sequence_hg19']})
#%%
final_df_uc.to_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis"
                      r"\main\Library_assembly\output\pre_satmut_uc.csv")

#%%
uc_df_tiles['sequence_anc']=uc_df_tiles.apply(anc_retreive,axis=1)
