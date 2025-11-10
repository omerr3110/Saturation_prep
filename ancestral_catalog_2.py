import io
import os
import pandas as pd
import numpy as np
import time
import sys
sys.path.append("/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5")
#import quality_functions

#determine HC ancestral allele

def ancestral_decider(row):
    anc=None
    major_human=row.loc['Human']
    major_gorilla=row.loc['Gorilla']
    major_chimp=row.loc['Chimp']
    major_pongo=row.loc['Pongo']
    #print(major_human,major_gorilla,major_chimp,major_pongo,row.loc['Hg19'])
    if major_human is None: ### for python 2.7, pd.insa() on current version
        #print("none")
        major_human=row.loc['Hg19']
    #print(major_thuman)
    if major_human==major_chimp:
        #print("1")
        anc=major_human
    elif(major_gorilla==major_human or major_gorilla==major_chimp) and not(major_gorilla is None):
        #print(~(major_gorilla is None))
        anc=major_gorilla
    elif(major_pongo==major_human or major_pongo==major_chimp) and not(major_pongo is None):
        anc=major_pongo
        #print("3")
    else:
        anc=major_human
        #print("4")
    return anc

# get from cmd: 1. chromosome; 2. chr. portion 3. settings file name
print(sys.argv)
chromosome = sys.argv[1]
chromosome_portion = sys.argv[2]

catalog_path=f"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/{chromosome}/{chromosome_portion}/major_{chromosome}_{chromosome_portion}_full.txt"
#catalog_path="/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/1/1/major_1_1_1.txt"
try:
    catalog_df=pd.read_csv(catalog_path,sep="\t",index_col=0,na_values=[None],header=None)
    if len(catalog_df)>0:
        catalog_df.columns=['chr','pos','Human','Chimp','Gorilla','Pongo','Hg19']
    catalog_df=catalog_df.replace({np.nan:None})
    new_catalog_path=f"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/{chromosome}/{chromosome_portion}/major_{chromosome}_{chromosome_portion}_full_with_anc.txt"
    catalog_df['HC_anc']=catalog_df.apply(ancestral_decider,axis=1)
    catalog_df.to_csv(new_catalog_path,sep="\t")
except pd.errors.EmptyDataError:
    print("empty file")




