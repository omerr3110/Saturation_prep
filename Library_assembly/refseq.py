import pandas as pd
#import datatable as dt
import os
#import xlrd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from functions import seq_retrieve
from functions import df_tiling
#%%
refseq_df=pd.read_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\Refseq"
                        r"\refseq_gene_list_fixed.txt",sep="\t")
#%%
filter_threshold=2000
window_size=135
#%%
filtered_refseq_df=refseq_df[refseq_df['length']<=filter_threshold]
#%%

#%%
tiles_df=df_tiling(filtered_refseq_df,window_size,270)
tiles_df['type']="refseq"
print(tiles_df.shape)
#%%
sns.histplot(data=refseq_df.query('length<=5000'),x='length')
plt.show()
#%%
final_df_refseq=pd.DataFrame(data={"chr":tiles_df['chromosome'],"start":tiles_df['oligo_start'],"end":tiles_df['oligo_end']
                                ,"group":"refseq","name":tiles_df['name'],'external_ID':tiles_df['ID']})
#%%
genome_file="hg19.fa"
records = list(SeqIO.parse(genome_file, "fasta"))
#%%
final_df_refseq['sequence_hg19']=final_df_refseq.apply(func=lambda x: seq_retrieve(x['chr'],x['start'],x['end'],records),axis=1)

#%%
final_df_refseq.to_csv(fr"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\output\pre_satmut_refseq_{filter_threshold}_{window_size}.csv")

