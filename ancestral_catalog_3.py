import io
import os
import pandas as pd
import numpy as np
import time
import sys
sys.path.append("/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5")
#import quality_functions

# reconstruct HCG allele

def ancestral_decider_hcg(row):
    HCG_anc=None
    HC_anc=row.loc['HC_anc']
    major_gorilla=row.loc['Gorilla']
    major_pongo=row.loc['Pongo']
    #print(major_human,major_gorilla,major_chimp,major_pongo,row.loc['Hg19'])
    if HC_anc==major_gorilla:
        #print("1")
        HCG_anc=HC_anc
    elif(major_pongo==major_gorilla or major_pongo==HCG_anc) and not(major_pongo is None):
        HCG_anc=major_pongo
    else:
        HCG_anc=HC_anc
        #print("4")
    return HCG_anc


print(sys.argv)
chromosome = sys.argv[1]
chromosome_portion = sys.argv[2]

catalog_path=f"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/{chromosome}/{chromosome_portion}/major_{chromosome}_{chromosome_portion}_full_with_anc.txt"
#catalog_path="/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/major_1_1_1.txt"
try:
    catalog_df=pd.read_csv(catalog_path,sep="\t",index_col=0,na_values=[None],header=[0])
    #if len(catalog_df)>0:
        #catalog_df.columns=['chr','pos','Human','Chimp','Gorilla','Pongo','Hg19',"HC_anc"]
    catalog_df=catalog_df.replace({np.nan:None})
    new_catalog_path=f"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/{chromosome}/{chromosome_portion}/major_{chromosome}_{chromosome_portion}_full_with_anc_HCG.txt"
    catalog_df['HCG_anc']=catalog_df.apply(ancestral_decider_hcg,axis=1)
    catalog_df.to_csv(new_catalog_path,sep="\t")
except FileNotFoundError:
    print("No file")
except pd.errors.EmptyDataError:
    print("empty file")
