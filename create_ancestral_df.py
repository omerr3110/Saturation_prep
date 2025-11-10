import io
import os
import Bio
import pandas as pd
import numpy as np
import time
import sys
from Bio import SeqIO
from functions import seq_retrieve
import regex as re
import math
import pickle
def ancestral_modify(row):
    seq=row['sequence_hg19']
    chrom=re.findall("\d+|X|Y",row['chr'])[0]
    start=row['start']
    end=row['end']
    mult_flag=0
    portion_st=start//(5*10**6)+1
    portion_end=end//(5*10**6)+1
    if portion_st!=portion_end:
        print("Multiple Portions")
        print(row['name'],chrom,portion_st,portion_end)
        mult_flag=1
    anc_path=rf"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/{chrom}/{portion_st}/major_{chrom}_{portion_st}_full_with_anc_HCG.txt"
    anc_path2=rf"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/{chrom}/{portion_end}/major_{chrom}_{portion_end}_full_with_anc_HCG.txt"
    #print(anc_path)
    try:
        anc_df = pd.read_csv(anc_path,sep="\t",index_col=0)
        if mult_flag==1:
            anc_df2 = pd.read_csv(anc_path2, sep="\t", index_col=0)
            anc_df=pd.concat([anc_df,anc_df2])
    except FileNotFoundError:
        print(anc_path, row['name'])
        print("file not found")
        return None, None, None, None, None,0
    anc_df['pos']=anc_df['pos'].apply(lambda x: int(x))
    anc_df=anc_df.query('pos<=@end and pos>=@start')
    #print(anc_df.shape)
    if len(anc_df)<=0:
        print(anc_path, row['name'])
        print("No changes")
        return None, None, None, None, None,0
    seq_list_pongo=list(seq)
    seq_list_gorilla = list(seq)
    seq_list_chimp = list(seq)
    seq_list_ancestral = list(seq)
    seq_list_hcg=list(seq)
    gen_per=len(anc_df)/270
    for i,nuc in enumerate(seq_list_ancestral):
        curr_nuc=nuc.upper() ###maybe can be deleted because all sequences are in upper case
        pos_df=anc_df.query('pos==@i+@start')
        if len(pos_df)>0:
            hc_anc=pos_df["HC_anc"].values[0]
            hcg_anc=pos_df["HCG_anc"].values[0]
            pon=pos_df["Pongo"].values[0]
            gor=pos_df["Gorilla"].values[0]
            chi = pos_df["Chimp"].values[0]
            #if chi==True:
                #print("here")
            if hc_anc!=curr_nuc:
                seq_list_ancestral[i]=hc_anc
            if hcg_anc!=curr_nuc:
                seq_list_hcg[i]=hcg_anc
            if pon != curr_nuc and type(pon) == str:  # need to validate with Simon
                seq_list_pongo[i] = pon
            if gor != curr_nuc and type(gor) == str:
                seq_list_gorilla[i] = gor
            if chi != curr_nuc and type(chi) == str:
                seq_list_chimp[i] = chi
    ancestral_seq="".join(seq_list_ancestral)
    pongo_seq="".join(seq_list_pongo)
    chimp_seq="".join(seq_list_chimp)
    gorilla_seq="".join(seq_list_gorilla)
    hcg_seq="".join(seq_list_hcg)
    return ancestral_seq,pongo_seq,chimp_seq,gorilla_seq,hcg_seq,gen_per

def ancestral_modify_local(row):
    seq=row['sequence_hg19']
    chrom=re.findall("\d+|X|Y",row['chr'])[0]
    start=row['start']
    end=row['end']
    portion_st=start//(5*10**6)+1
    portion_end=end//(5*10**6)+1
    if portion_st!=portion_end:
        print("Multiple Portions")
        print(row['Name'],chrom,portion_st,portion_end)
    anc_path=rf"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\output\{chrom}\{portion_st}\major_{chrom}_{portion_st}_full_with_anc_HCG.txt"
    print(anc_path)
    try:
        anc_df = pd.read_csv(anc_path,sep="\t",index_col=0)
    except FileNotFoundError:
        print(anc_path,row['name'])
        print("file not found")
        return seq
    anc_df['pos']=anc_df['pos'].apply(lambda x: int(x))
    anc_df=anc_df.query('pos<=@end and pos>=@start')
    #print(anc_df.shape)
    if len(anc_df)<=0:
        #print("no changes")
        return seq
    seq_list_pongo=list(seq)
    seq_list_gorilla = list(seq)
    seq_list_chimp = list(seq)
    seq_list_ancestral = list(seq)
    seq_list_hcg = list(seq)
    #print(anc_df)
    for i,nuc in enumerate(seq_list_ancestral):
        curr_nuc=nuc.upper() ###maybe can be deleted because all sequences are in upper case
        pos_df=anc_df.query('pos==@i+@start')
        if len(pos_df)>0:
            hc_anc = pos_df["HC_anc"].values[0]
            hcg_anc = pos_df["HCG_anc"].values[0]
            pon=pos_df["Pongo"].values[0]
            gor=pos_df["Gorilla"].values[0]
            chi = pos_df["Chimp"].values[0]
            if hc_anc != curr_nuc:
                seq_list_ancestral[i] = hc_anc
            if hcg_anc != curr_nuc:
                seq_list_hcg[i] = hcg_anc
            if pon!=curr_nuc and type(pon)==str: #need to validate with Simon
                seq_list_pongo[i]=pon
            if gor!=curr_nuc and type(gor)==str:
                seq_list_gorilla[i]=gor
            if chi!=curr_nuc and type(chi)==str:
                seq_list_chimp[i]=chi
    ancestral_seq="".join(seq_list_ancestral)
    pongo_seq="".join(seq_list_pongo)
    chimp_seq="".join(seq_list_chimp)
    gorilla_seq="".join(seq_list_gorilla)
    hcg_seq = "".join(seq_list_hcg)
    return ancestral_seq,pongo_seq,chimp_seq,gorilla_seq,hcg_seq

def archaic_modify(row, d):
    seq = row['sequence_hg19']
    chrom = row['chr']
    start = row['start']
    end = row['end']
    new_d = {k: v for k, v in d.items() if k[0] == chrom and end >= k[1] >= start}
    if len(new_d) == 0:
        return seq
    archaic_seq_list = list(seq)
    for k in new_d.keys():
        var_nuc = new_d[k]
        oligo_pos = k[1] - start
        archaic_seq_list[oligo_pos] = var_nuc
    archaic_seq = "".join(archaic_seq_list)
    return archaic_seq

def H_ancestral_modify(row, d):
    seq = row['sequence_hg19']
    chrom = row['chr']
    start = row['start']
    end = row['end']
    new_d = {k: v for k, v in d.items() if k[0] == chrom and end >= k[1] >= start}
    if len(new_d) == 0:
        return seq
    H_anc_seq_list = list(seq)
    for k in new_d.keys():
        var_nuc = new_d[k]
        oligo_pos = k[1] - start
        H_anc_seq_list[oligo_pos] = var_nuc
    H_anc_seq = "".join(H_anc_seq_list)
    return H_anc_seq


print(sys.argv)
path = sys.argv[1]
with open(r"/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5"
          r"/archaic_derived_variants_dict_str.pkl", 'rb') as f:
    dict1=pickle.load(f)
with open(r"/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5"
          r"/Human_ancestral_dict.pkl", 'rb') as f:
    dict2 = pickle.load(f)
#path = r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\output\pre_satmut_test.csv"
#path = r"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/satmut_pilot.csv"
df=pd.read_csv(path,index_col=0)
sequences=df.apply(ancestral_modify,axis=1)
df['sequence_ancestral']=sequences.apply(lambda x: x[0])
df['sequence_pongo']=sequences.apply(lambda x: x[1])
df['sequence_chimp']=sequences.apply(lambda x: x[2])
df['sequence_gorilla']=sequences.apply(lambda x: x[3])
df['sequence_HCG_ancestral']=sequences.apply(lambda x: x[4])
df['genotyped_percentage']=sequences.apply(lambda x: x[5])
df['sequence_AH'] = df.apply(archaic_modify, axis=1, args=(dict1|dict2,))
df['sequence_H_ancestral'] = df.apply(H_ancestral_modify, axis=1, args=(dict2,))

new_path=re.split(pattern=".csv",string=path)[0]+"_with_all_perc.csv"
df.to_csv(new_path)


