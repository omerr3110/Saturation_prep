#%%
import pandas as pd
import regex as re
import numpy as np
import pickle
#%%
seq_df=pd.read_excel(r"C:\Users\omerro\Dropbox (Weizmann Institute)\Omer-David\L4\L4_v5.xlsx")
#%%
oligo_class=seq_df['class']
#%%
pre_satmut_df=seq_df[seq_df['class'].str.contains("satmut\\.")]
pilot_df=seq_df[seq_df['class'].str.contains("Pilot")]
#%%
pre_satmut_df['class'].value_counts()

#%%
#QA on haqers:
haqer_df=pre_satmut_df[pre_satmut_df['class']=='satmut.preMPRA.HAQERs']
print(sum(haqer_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1))) #1012 as Simon's
print(sum(haqer_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1))) #104

#%%
har_df=pre_satmut_df[pre_satmut_df['class']=='satmut.preMPRA.HAR']
print(sum(har_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1))) #33 as Simon's
print(sum(har_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1))) #104 as Simon's

#%%
ref_df=pre_satmut_df[pre_satmut_df['class']=='satmut.preMPRA.refseq']
print(sum(ref_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1))) #0
print(sum(ref_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1))) #1553
#%%
SNV_df=pre_satmut_df[pre_satmut_df['class']=='satmut.preMPRA.MH.divergent.SNV']
print(sum(SNV_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1))) #0
print(sum(SNV_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1))) #990

#%%
pheno_df=pre_satmut_df[pre_satmut_df['class']=='satmut.preMPRA.Phenotype.related']
print(sum(pheno_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1))) #0
print(sum(pheno_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1))) #549

#%%
psych_df=pre_satmut_df[pre_satmut_df['class']=='satmut.preMPRA.random.neurons']
print(sum(psych_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1))) #1 - checks out
print(sum(psych_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1))) #763

#%%
chond_df=pre_satmut_df[pre_satmut_df['class']=='satmut.preMPRA.Random.chond']
print(sum(chond_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1))) #0
print(sum(chond_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1))) #362

#%%
disease_df=pre_satmut_df[pre_satmut_df['class']=='satmut.preMPRA.disease.linked']
print(sum(disease_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1)))  # 9
print(sum(disease_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1)))  # 35
print(sum(disease_df.apply(lambda x: pd.isna(x['sequence1_with_adapters']),axis=1)))  # 0
print(sum(disease_df.apply(lambda x: pd.isna(x['sequence2_with_adapters']),axis=1)))  # 109
print(sum(disease_df.apply(lambda x: pd.isna(x['sequence3_with_adapters']),axis=1)))  # 110
print(sum(disease_df.apply(lambda x: pd.isna(x['sequence4_with_adapters']),axis=1)))  # 110
#%%
uc_df=pre_satmut_df[pre_satmut_df['class']=='satmut.preMPRA.Ultra.Conserved']
print(sum(uc_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1))) #378
print(sum(uc_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1))) #0

changes_df=uc_df.loc[uc_df.apply(lambda x: not(pd.isna(x['derived_sequence_with_adapters'])),axis=1)]
uc_ids=changes_df['name'].apply(lambda x: re.findall("uc.\d+",x)[0])

#%%
print(sum(pilot_df.name.apply(lambda x: x.count("del"))))
print(sum(pilot_df.name.apply(lambda x: bool(re.findall("del",x)))))
print(sum(pilot_df.name.apply(lambda x: bool(re.findall(">",x)))))
print(sum(pilot_df.name.apply(lambda x: bool(re.findall("node",x)))))
print(sum(pilot_df.name.apply(lambda x: bool(re.findall(">",x)) and bool(re.findall("node",x)))))
node_perm_df=pilot_df[pilot_df.name.apply(lambda x: bool(re.findall(">",x)) and bool(re.findall("node",x)))]
#%%
missing_df=seq_df[seq_df['class']=="hh.missing.oligos"]
print(sum(missing_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1))) #0
print(sum(missing_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1))) #0
print(sum(missing_df.apply(lambda x: pd.isna(x['sequence1_with_adapters']),axis=1))) #0

#%%
hh_screen_df=seq_df[seq_df['class']=="hh.SCREEN"]
print(sum(hh_screen_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1))) #0
print(sum(hh_screen_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1))) #0
print(sum(hh_screen_df.apply(lambda x: pd.isna(x['sequence1_with_adapters']),axis=1))) #0
#%%

FABP7_df=seq_df[seq_df['class']=="FABP7"]
print(sum(FABP7_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1))) #16
print(sum(FABP7_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1))) #16
print(sum(FABP7_df.apply(lambda x: pd.isna(x['sequence1_with_adapters']),axis=1))) #0
print(sum(FABP7_df.apply(lambda x: pd.isna(x['sequence2_with_adapters']),axis=1))) #0
print(sum(FABP7_df.apply(lambda x: pd.isna(x['sequence3_with_adapters']),axis=1))) #2
print(sum(FABP7_df.apply(lambda x: pd.isna(x['sequence4_with_adapters']),axis=1))) #2

#%%
CGI_df=seq_df[seq_df['class']=="CGI.fix"]
print(sum(CGI_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1)))
print(sum(CGI_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1)))
print(sum(CGI_df.apply(lambda x: pd.isna(x['sequence1_with_adapters']),axis=1)))

#%%
ori_df=seq_df[seq_df['class']=="promoter.orientation.fix"]
print(sum(ori_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1)))
print(sum(ori_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1)))
print(sum(ori_df.apply(lambda x: pd.isna(x['sequence1_with_adapters']),axis=1)))
#%%
indel_df=seq_df[seq_df['class']=="MH.AH.INDELs"]
print(sum(indel_df.apply(lambda x: pd.isna(x['derived_sequence_with_adapters']),axis=1)))
print(sum(indel_df.apply(lambda x: pd.isna(x['ancestral_sequence_with_adapters']),axis=1)))
print(sum(indel_df.apply(lambda x: pd.isna(x['sequence1_with_adapters']),axis=1)))
#%%
seq1_df=pilot_df[pilot_df['external_ID']==1]
#%%
seq1_nodes_df=seq1_df.loc[pilot_df.name.apply(lambda x: bool(re.findall("node",x)))]
#%%
seq1_cat=pd.read_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main"
                     r"\Library_assembly\output\12\22\major_12_22_full_with_anc_HCG.txt",sep="\t")
#%%
seq1_cat=seq1_cat.query('108511253<pos<108511522')
#%%
seq1_cat=seq1_cat.loc[seq1_cat.apply(lambda x: not(x['Chimp']== x['Gorilla']== x['Pongo']== x['Hg19']== x['HC_anc']==
                                               x['HCG_anc']),axis=1)]
#%%
seq1_cat['relative_pos']=seq1_cat['pos']-108511253
#%%
with open(r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\Tomas_Apes_Catalog_parsing_6_CGI_catalog"
          r"\FILTER_INCLUSIVE_STRICT_0_3_good_samples3_delta\INDEL_filtered"
          r"\Human_derived_dict.pkl", 'rb') as f:
    dict3 = pickle.load(f)

with open(r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\Tomas_Apes_Catalog_parsing_5"
          r"\archaic_derived_variants_dict_str.pkl", 'rb') as f:
    dict1=pickle.load(f)
with open(r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\Tomas_Apes_Catalog_parsing_5"
          r"\Human_ancestral_dict.pkl", 'rb') as f:
    dict2 = pickle.load(f)
#%%
merged_dict=dict1|dict2|dict3
#%%
new_d = {k: v for k, v in merged_dict.items() if k[0] == "chr12" and 108511522 >= k[1] >= 108511253}
#%%
scrambled_df=pd.read_excel(r"C:\Users\omerro\Dropbox (Weizmann Institute)\Omer-David\L4\L4_v5.xlsx"
                           ,sheet_name="ctrl_scrambled")
#%%
def GC_calc(seq):
    #print(seq)
    return (seq.count("g")+seq.count("c")+seq.count("G")+seq.count("C"))/270
def origin_calc(row):
    ori_name=row['source_oligo_name']
    if ori_name==row['name_ancestral_sequence']:
        gc_origin=GC_calc(row['ancestral_sequence'])
    elif ori_name == row['name_derived_sequence']:
        gc_origin = GC_calc(row['derived_sequence'])
    else:
        gc_origin = GC_calc(row['sequence1'])
    gc_scrambled=GC_calc(row['shuffled_seq'])
    return gc_scrambled==gc_origin
#%%
print(sum(scrambled_df.apply(origin_calc,axis=1)))