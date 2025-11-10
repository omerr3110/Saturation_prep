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
import regex as re
import Bio
from functions import anc_retreive
#%%
oligo_path=r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\disease" \
           r"\oligos_coordinates_updated.xlsx"
#%%
oligo_df=pd.read_excel(oligo_path)
oligo_df['reference_ID']=oligo_df['reference_ID'].apply(lambda x: x.replace("_","."))

#%%
filtered_oligos_df=oligo_df.query('`nt change`!="manual"').copy()
manual_oligos_df=oligo_df.query('`nt change`=="manual"').copy()
#%%
filtered_oligos_tiles=df_tiling(filtered_oligos_df,50,270)
#%%
#%%
final_df_disease=pd.DataFrame(data={"chr":filtered_oligos_tiles['chromosome'],
                                  "start":filtered_oligos_tiles['oligo_start'],"end":filtered_oligos_tiles['oligo_end'],
                                  "group":"disease linked","name":filtered_oligos_tiles['name'],
                                    "external_ID":filtered_oligos_tiles['ID']})
#%%
genome_file="hg19.fa"
records = list(SeqIO.parse(genome_file, "fasta"))
#%%
final_df_disease['sequence_hg19']=final_df_disease.apply(func=lambda x: seq_retrieve(x['chr'],x['start'],x['end'],records),axis=1)

#%%
variants=filtered_oligos_df['nt change'].apply(lambda x: re.split(r">|;",x))
variants=variants.reset_index()
#%%
final_df_disease['variants']=variants['nt change']
#%%
def variant_retreive(row):
    seq=row['sequence_hg19']
    possibilities=row['variants']
    if len(possibilities)==2:
        hg19=seq[135]
        hg19=hg19.upper()
        if hg19==possibilities[0]:
            variant=possibilities[1]
        elif hg19==possibilities[1]:
            variant = possibilities[0]
        else:
            variant="TBD"
            print(row['name'])
        seq_list=list(seq)
        seq_list[135]=variant
        new_seq = Bio.Seq.Seq("".join(seq_list))
        return new_seq
    else:
        # used to take care of oligo with two possible variants in the same position
        #print(possibilities)
        variant=possibilities[1]
        seq_list = list(seq)
        seq_list[135] = variant
        new_seq1 = Bio.Seq.Seq("".join(seq_list))
        variant=possibilities[3]
        seq_list = list(seq)
        seq_list[135] = variant
        new_seq2 = Bio.Seq.Seq("".join(seq_list))
        return (new_seq1,new_seq2)

#%%

final_df_disease['variant_sequence_hg19']=final_df_disease.apply(variant_retreive,axis=1)
#%%
double_df=final_df_disease[final_df_disease['variant_sequence_hg19'].apply(lambda x: len(x)==2)].copy()
single_df=final_df_disease[final_df_disease['variant_sequence_hg19'].apply(lambda x: len(x)==270)].copy()

#%%
seq1=double_df['variant_sequence_hg19'][87][0]
seq2=double_df['variant_sequence_hg19'][87][1]
#%%
double_df['variant_sequence_hg19']=[seq1]
double_df['variant2_sequence_hg19']=[seq2]

#%%
auto_df=pd.concat([single_df,double_df],ignore_index=True)
auto_df.drop(columns=['variants'],inplace=True)

#%%
manual_oligos_tiles=df_tiling(manual_oligos_df,50,270)
manual_df=pd.DataFrame(data={"chr":manual_oligos_tiles['chromosome'],
                                  "start":manual_oligos_tiles['oligo_start'],"end":manual_oligos_tiles['oligo_end'],
                                  "group":"disease linked","name":manual_oligos_tiles['name'],
                             'external_ID':manual_oligos_tiles['ID']})
#%%
manual_df['sequence_hg19']=manual_df.apply(func=lambda x: seq_retrieve(x['chr'],x['start'],x['end'],records),axis=1)
#%%
comments=manual_oligos_df.comments.reset_index()
#manual_df['variants']=comments['comments']
manual_df['variant_sequence_hg19']=manual_oligos_df.disease_seq.apply(lambda x: Bio.Seq.Seq(x)).reset_index()['disease_seq']

#%%
complete_df_disease=pd.concat([auto_df,manual_df],axis=0,ignore_index=True)
#%%
complete_df_disease.to_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis"
                      r"\main\Library_assembly\output\pre_satmut_disease.csv")
