import io
import os
import pandas as pd
import numpy as np
import time
import sys
sys.path.append("/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5")
#import quality_functions

# reconstruct HCGP allele

def ancestral_decider_hcgp(row):
    hcg_anc=row.loc['HCG_anc']
    major_pongo=row.loc['Pongo']
    gibbon_allele=row.loc['gibbon']
    if hcg_anc==major_pongo:
        hcgp_anc=hcg_anc
    elif gibbon_allele==major_pongo:
        hcgp_anc = major_pongo
    elif gibbon_allele==hcg_anc:
        hcgp_anc = hcg_anc
    else:
        hcgp_anc="TBD"
        print("TBD",row)
    return hcgp_anc


print(sys.argv)
chromosome = sys.argv[1]
chromosome_portion = sys.argv[2]

catalog_path=f"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/{chromosome}/{chromosome_portion}/major_{chromosome}_{chromosome_portion}_full_with_gibbon.txt"
#catalog_path="/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/1/1/major_1_1_full_with_gibbon.txt.txt"
try:
    catalog_df=pd.read_csv(catalog_path,sep="\t",index_col=0,na_values=[None])
    if len(catalog_df)>0:
        catalog_df.columns=['chr','pos','Human','Chimp','Gorilla','Pongo','Hg19',"HC_anc",'HCG_anc','gibbon']
    catalog_df=catalog_df.replace({np.nan:None})
    new_catalog_path=f"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/{chromosome}/{chromosome_portion}/major_{chromosome}_{chromosome_portion}_full_with_anc_HCGP.txt"
    catalog_df['HCGP_anc']=catalog_df.apply(ancestral_decider_hcgp,axis=1)
    catalog_df.to_csv(new_catalog_path,sep="\t")
except FileNotFoundError:
    print("No file")
except pd.errors.EmptyDataError:
    print("empty file")