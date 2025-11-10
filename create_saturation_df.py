mport io
import os
import Bio
import pandas as pd
import numpy as np
import time
import sys
from Bio import SeqIO
from functions import seq_retrieve
import regex as re
import pickle
#path=r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\output\pre_satmut_test.csv"
#/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/pre_satmut_test.csv
path=sys.argv[1]
df=pd.read_csv(path,index_col=0)
with open(r"/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5"
          r"/archaic_derived_variants_dict_str.pkl", 'rb') as f:
    dict1=pickle.load(f)
with open(r"/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5"
          r"/MH_derived_dict.pkl", 'rb') as f:
    dict2 = pickle.load(f)
with open(r"/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_6_CGI_catalog"
          r"/FILTER_INCLUSIVE_STRICT_0_3_good_samples3_delta/INDEL_filtered"
          r"/Human_derived_dict.pkl", 'rb') as f:
    dict3 = pickle.load(f)
merged_dict=dict(dict1.items()|dict2.items()|dict3.items())


def position_calc(var_dict,chrom,start,end):
    new_var_dict={k: v for k, v in var_dict.items() if k[0] == chrom and end >= k[1] >= start}
    norm_positions= {k[1]-start: v for k,v in new_var_dict.items()}
    filtered_norm_positions={k:v for k,v in norm_positions.items() if k<35 or k>=235}
    return filtered_norm_positions


def satmut_df(row):
    nucs={"A","T","C","G"}
    oligo_change=[]
    sequences=[]
    oligo_type=["node" for _ in range(10)]+["saturation" for _ in range(800)]
    #print(curr_name)
    seq=list(row['sequence_ancestral'])
    oligo_change.append("HC-ancestral")
    oligo_change.append("HC-ancestral-rep")
    final_seq="".join(seq)
    sequences.extend([final_seq,final_seq])
    seq_gorilla=row['sequence_gorilla']
    oligo_change.append("Gorilla")
    sequences.append(seq_gorilla)
    seq_pongo = row['sequence_pongo']
    oligo_change.append("Pongo")
    sequences.append(seq_pongo)
    seq_chimp = row['sequence_chimp']
    oligo_change.append("Chimp")
    sequences.append(seq_chimp)
    seq_AH = row['sequence_AH']
    oligo_change.append("AH")
    sequences.append(seq_AH)
    seq_HCG_anc = row['sequence_HCG_ancestral']
    oligo_change.append("HCG-ancestral")
    sequences.append(seq_HCG_anc)
    seq_hg19 = row['sequence_hg19']
    oligo_change.append("Hg19")
    sequences.append(seq_hg19)
    seq_HCGP_anc = row['sequence_HCGP_ancestral']
    oligo_change.append("HCGP-ancestral")
    sequences.append(seq_HCGP_anc)
    seq_H_anc = row['sequence_H_ancestral']
    oligo_change.append("H-ancestral")
    sequences.append(seq_H_anc)
    ex_pos=position_calc(merged_dict,row['chr'],row['start'],row['end'])
    for i in range(35,235):
        curr_nuc=seq[i].upper()
        for nuc in nucs:
            if nuc!=curr_nuc:
                change=f"{i+1}({curr_nuc}>{nuc})" #previously was i and not i+1, means that satmut135 is actually the 136th position
                oligo_change.append(change)
                seq[i]=nuc
                final_seq = "".join(seq)
                sequences.append(final_seq)
                seq[i]=curr_nuc
        seq.pop(i)
        final_seq = "".join(seq)
        sequences.append(final_seq)
        seq.insert(i,curr_nuc)
        change =f"{i+1}({curr_nuc}.del)"
        oligo_change.append(change)
    for k,v in ex_pos.items():
        curr_nuc = seq[k]
        change = f"{k}({curr_nuc}>{v})"
        oligo_change.append(change)
        oligo_type.append("saturation")
        seq[k] = v
        final_seq = "".join(seq)
        sequences.append(final_seq)
        seq[k] = curr_nuc
    print(len(sequences),len(sequences),len(oligo_change),len(oligo_type))
    final_df=pd.DataFrame(data={"chr":row['chr'],"start":row['start'],"end":row['end'],"group":row['group'],
                                'ID_satmut':row['pilot_id'],"original_name":row['name'],'oligo_change':oligo_change,
                                'oligo_type':oligo_type,"sequence":sequences})
    return final_df



res=df.apply(satmut_df,axis=1)
full_df=pd.concat(res.to_list(),axis=0)
new_path=re.split(pattern=".csv",string=path)[0]+"_saturation_oligos.csv"
full_df.to_csv(new_path)