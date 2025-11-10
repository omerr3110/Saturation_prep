# %%
import pandas as pd
import datatable as dt
import os
#import xlrd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from functions import seq_retrieve
from functions import df_tiling
import Bio.Seq
#%%
#def seq_retrieve(chr,start,end,records):
    #curr_chr=[record for record in records if record.id==chr]
    #return curr_chr[0].seq[start-1:end]


# %%
all_zh_df = pd.read_excel(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly"
                          r"\HARs\supplementary data zoonomia\science.abm1696_table_s1.xlsx", sheet_name="zooHARs")
#zh_mpra_df = pd.read_excel(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\Saturation Mutagenesis\main\Library assembly"
#                           r"\HARs\supplementary data zoonomia\science.abm1696_table_s8.xlsx",
#                           sheet_name="zooHARMPRAreplicateMean")


#%%
positions=all_zh_df.loc[:,['chrom','start','end']].copy()
positions['start']=positions['start']-1
#%%
positions.to_csv("to_lift_hars.bed",header=False,sep="\t",index=False)
#%%
#lifted file is 0-based
lifted_df=pd.read_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\HARs"
                      r"\lifted_hars.bed", sep='\t',header=None)
# add 1 to start to create one based
all_zh_df['chrom']=lifted_df[0]
all_zh_df['start']=lifted_df[1]+1
all_zh_df['end']=lifted_df[2]

#%%
#active_zh=zh_mpra_df.sort_values(by="RNA:DNA_ratio",ascending=False)['HAR'][:139].to_list()
#active_zh_df=pd.merge(all_zh_df,zh_mpra_df,left_on='simple_name', right_on='HAR').sort_values(by="RNA:DNA_ratio",ascending=False)[:139]


#%%
all_zh_df=all_zh_df.drop(columns=['name'])
#%%
all_zh_df.rename(columns={"simple_name":"name"},inplace=True)
#%%
tiles_df=df_tiling(all_zh_df,135,270)
#%%
tiles_df['type']="HAR"

#%%

#%%
final_df_har=pd.DataFrame(data={"chr":tiles_df['chromosome'],"start":tiles_df['oligo_start'],"end":tiles_df['oligo_end']
                                ,"group":tiles_df['type'],"name":tiles_df['name'],'external_ID':tiles_df['ID']})

#%%
genome_file="hg19.fa"
records = list(SeqIO.parse(genome_file, "fasta"))
#%%
final_df_har['sequence_hg19']=final_df_har.apply(func=lambda x: seq_retrieve(x['chr'],x['start'],x['end'],records),axis=1)

#%%
final_df_har.to_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\output\pre_satmut_HARs.csv")
#%%
all_zh_df['Length']=all_zh_df.apply(lambda x: x[2]-x[1],axis=1)
plot=sns.histplot(data=all_zh_df,x='Length')
plot.set(title='Histogram of HARs length')
plt.show()
print(all_zh_df['Length'].max(),all_zh_df['Length'].min(),all_zh_df['Length'].median())

#%%
# test for ancestral translation
zoohar311=all_zh_df.iloc[311,:]
zoohar311_seq=seq_retrieve(zoohar311['chrom'],zoohar311['start'],zoohar311['end'],records)
#%%

final_df_har