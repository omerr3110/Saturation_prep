import pandas as pd
#import datatable as dt
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
genome_file="hg19.fa"
records = list(SeqIO.parse(genome_file, "fasta"))
#%%

#%%
evc2_df=pd.read_excel(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis"
                      r"\main\Library_assembly\Divergent SNV Phenotypes\evc2_supp_data.xlsx",sheet_name="S28",
                      header=2)
#%%
evc2_df=evc2_df[evc2_df['Linked Gene']=='EVC2']
#%%
# positions are the 1-based lifted (hg19) original positions (hg38) from evc2_df
evc2_df_final=pd.DataFrame(data={"chr":["chr4","chr4","chr4"],"start":[4991950,5366728,5671799],
                                 "end":[4994127,5369127,5672167],
                                 'name':["EVC2enhancer.1","EVC2enhancer.2","EVC2enhancer.3"]},index=[0,1,2])
evc2_df_tiles=df_tiling(evc2_df_final,step_size=135,oligo_length=270)

#%%
# positions are taken from supplementary data of fzd8 paper, already in hg19, one-based
hare5_df=pd.DataFrame(data={"chr":"chr10","start":36238121,"end":36239339,'name':"HARE5"},index=[0])
#%%
hare5_df_tiles=df_tiling(hare5_df,135,270)
#%%

#%%
# positions are the 1-based lifted (hg19) original positions (hg38) from kiaa1217 paper supplementary data
kiaa1217_df=pd.DataFrame(data={"chr":"chr10","start":24447422,"end":24447722,'name':"KIAA1217enhancer"},index=[0])
kiaa1217_df_tiles=df_tiling(kiaa1217_df,135,270)

#%%
# positions are taken from supplementary data of sox9 paper - the DMR that overlaps ec1.45 enhancer, already in hg19, one-based
sox9_df=pd.DataFrame(data={"chr":"chr17","start":68668482,"end":68674772},index=[0])
#%%
# run on server - currently not in use
sox9_bed_df=sox9_df
sox9_bed_df['start']=sox9_bed_df['start']-1
sox9_bed=pybedtools.BedTool.from_dataframe(sox9_bed_df)
var_cat_human_bed=pybedtools.BedTool(r"/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5/"
                          r"Human_derived_0.3_genotyped_good_samples3/INDEL_filtered/catalog_gnomad_0_95_filtered.bed")
intersect_human_bed= sox9_bed.intersect(var_cat_human_bed,wb=True)
intersect_human_df = intersect_human_bed.to_dataframe()
intersect_human_df=intersect_human_df.iloc[:,[0,1,2,6,7]]
intersect_human_df['derive']='Human'
var_cat_chimp_bed=pybedtools.BedTool(r"/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5"
                          r"/Chimp_derived_0.3_genotyped_good_samples3/INDEL_filtered/catalog.bed")
intersect_chimp_bed= sox9_bed.intersect(var_cat_chimp_bed,wb=True)
intersect_chimp_df = intersect_chimp_bed.to_dataframe()
intersect_chimp_df=intersect_chimp_df.iloc[:,[0,1,2,6,7]]
intersect_chimp_df['derive']='Chimp'



var_cat_gorilla_bed=pybedtools.BedTool(r"/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5/"
                                       r"Gorilla_derived_0.3_genotyped_good_samples3/INDEL_filtered/catalog.bed")
intersect_gorilla_bed= sox9_bed.intersect(var_cat_gorilla_bed,wb=True)
intersect_gorilla_df = intersect_gorilla_bed.to_dataframe()
intersect_gorilla_df=intersect_gorilla_df.iloc[:,[0,1,2,6,7]]
intersect_gorilla_df['derive']='Gorilla'

var_cat_pongo_bed=pybedtools.BedTool(r"/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5/"
                                     r"Orangutan_derived_0.3_genotyped_good_samples3/INDEL_filtered/catalog.bed")
intersect_pongo_bed= sox9_bed.intersect(var_cat_pongo_bed,wb=True)
intersect_pongo_df = intersect_pongo_bed.to_dataframe()
intersect_pongo_df=intersect_pongo_df.iloc[:,[0,1,2,6,7]]
intersect_pongo_df['derive']='Pongo'
var_df=pd.concat([intersect_human_df,intersect_chimp_df,intersect_gorilla_df,intersect_pongo_df],ignore_index=True)
var_df.columns=["chr","start","end","mut","prob","origin"]
var_df.to_csv(r"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/"
              r"Divergent SNV Phenotypes/sox9_DMR_variants.csv")

#%%
var_df=pd.read_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main"
                   r"\Library_assembly\Divergent SNV Phenotypes\sox9_DMR_variants.csv",index_col=0)
#%%
var_df=var_df.sort_values(by="start")
#%%
# use sox9_df and not test_df because it overlaps with all the area with less oligos
#test_df=df_tiling(var_df,50,270)
sox9_df['name']='EC1.45.DMR'
sox9_df_tiles=df_tiling(sox9_df,135,270)

#%%
# one-based, hg19 can be filtered if needed
foxp2_df1=pd.read_excel(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main"
                   r"\Library_assembly\Divergent SNV Phenotypes\foxp2_enhancers.xlsx",header=1)
#%%
pos=foxp2_df1['Coordinates'].apply(lambda x: re.split(r':|-',x))
foxp2_df_final=pd.DataFrame(data={"chr":pos.apply(func=lambda x: x[0]),"start":pos.apply(func=lambda x: int(x[1])),
                                  "end":pos.apply(func=lambda x: int(x[2])),"name":foxp2_df1['HARs']})
foxp2_df_tiles1=df_tiling(foxp2_df_final,135,270)
#%%
# taken from tables 1 and 4 in the paper, one-based, hg19
foxp2_df2=pd.DataFrame(data={"chr":["chr7","chr7","chr7","chr7","chr7","chr7"],
                             "start":[113724817,114051220,114055454,113688009,114056845,114568454],
                             "end":[113726609,114055324,114056459,113688782,114058646,114572411],
                             "name":["FOXP2promoter.TSS1","FOXP2promoter.TSS2","FOXP2promoter.TSS3",
                                     "FOXP2enhancer.37","FOXP2enhancer.330","FOXP2enhancer.843"]})
foxp2_df_tiles2=df_tiling(foxp2_df2,135,270)
#%%
# Taken from the introduction section of the paper, one-based, hg19
foxp2_df3=pd.DataFrame(data={"chr":["chr7","chr7"],
                             "start":[114456873,114541370],
                             "end":[114463136,114543683],
                             "name":["FOXP2enhancerP","FOXP2enhancerD"]})
foxp2_df_tiles3=df_tiling(foxp2_df3,135,270)
#%%
# positions are the 3k boundries of a 1.8k conserved region that was tested, hg18->hg19, 1-based

foxp2_df4=pd.DataFrame(data={"chr":["chr7"],
                             "start":[114288164],
                             "end":[114291164],
                             "name":["POU3F2BSFOXP2"]})

foxp2_df_tiles4=df_tiling(foxp2_df4,135,270)

#%%
foxp2_df_tiles_all=pd.concat([foxp2_df_tiles1,foxp2_df_tiles2,foxp2_df_tiles3,foxp2_df_tiles4])


#%%
auts_df=pd.read_excel(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis"
                      r"\main\Library_assembly\Divergent SNV Phenotypes\AUTS2.xlsx",header=2)
auts_df=auts_df.iloc[:40,:]
#%%
filtered_auts2_df=auts_df.query("AEC in [24,40,10,13,21,29,32,3,4,5,6]").reset_index()
#evolutionary conserved noncoding intronic regions in AUTS2:HAR31(AEC24), HACNS174(AEC40) and HACNS369(AEC10)
#Positive mouse enhancers: 21,13
#Human: 32, 29
#Dyslexia: 3-6
#%%
position=filtered_auts2_df['Chromosomal position (hg18)']
position=position.str.replace(",","")
position=position.str.split(":|-")
#%%
#1-based after hg18->hg19 lift
lifted_auts=pd.read_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis"
                      r"\main\Library_assembly\Divergent SNV Phenotypes\auts2_lifted_1based.bed",sep="\t",
                        names=["position"])
#%%
position=lifted_auts.loc[:,'position']
position=position.str.replace(",","")
position=position.str.split(":|-")
final_lifted_auts=pd.DataFrame(data={"chr":position.apply(func= lambda x: x[0]),"start":position.apply(func= lambda x: int(x[1])),
                               "end":position.apply(func= lambda x: int(x[2]))})
#%%
final_auts_df=pd.concat([final_lifted_auts,filtered_auts2_df.loc[:,'AEC']],axis=1)
final_auts_df['AEC']=final_auts_df['AEC'].apply(lambda x: "AEC"+str(x))
#%%
auts2_df_tiles=df_tiling(final_auts_df,135,270)
#%%
# one-based, already hg19
cux1_df=pd.DataFrame(data={"chr":["chr7"],
                             "start":[101249296],
                             "end":[101249680],
                             "name":["CUX1HAR426"]})

#%%
cux1_df_tiles=df_tiling(cux1_df,135,270)

#%%
full_df=pd.concat([sox9_df_tiles,kiaa1217_df_tiles,evc2_df_tiles,hare5_df_tiles,foxp2_df_tiles_all,auts2_df_tiles,cux1_df_tiles])
final_df_pheno=pd.DataFrame(data={"chr":full_df['chromosome'],"start":full_df['oligo_start'],"end":full_df['oligo_end']
                                ,"group":"Phenotype_related","name":full_df['name'],
                                  'external_ID':full_df['ID']})
#%%
final_df_pheno['sequence_hg19']=final_df_pheno.apply(func=lambda x: seq_retrieve(x['chr'],x['start'],x['end'],records),axis=1)
#%%
final_df_pheno.to_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis"
                      r"\main\Library_assembly\output\pre_satmut_phenotype.csv")
#%%
pheno_df=pd.read_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis"
                      r"\main\Library_assembly\output\pre_satmut_phenotype.csv",index_col=0)
#%%

#%%
test=pheno_df.iloc[:3].copy()
#%%


#%%
