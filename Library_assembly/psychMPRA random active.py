#%%
import pandas as pd
#import pybedtools
from functions import df_tiling
import regex as re
import numpy as np
from functions import seq_retrieve
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
#%%
genome_file="hg19.fa"
records = list(SeqIO.parse(genome_file, "fasta"))
#%%
#psych_df=pd.read_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\psychMPRA"
#                     r"\annotations-da-primary.tsv",sep="\t")
psych_df=pd.read_csv(r"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/psychMPRA/annotations-da-primary.tsv",sep="\t")

#%%
active_df=psych_df.query("is_active==1").copy()

#%%
'''Chondrocytes analysis, run on server'''
# read atac seq data for chondrocytes
#atac_seq_path=r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\psychMPRA\ATC.Bon.50.AllAg.Chondrocytes.bed"
atac_seq_path=r"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/psychMPRA/ATC.Bon.50.AllAg.Chondrocytes.bed"

ch_atac_data_df=pd.read_csv(atac_seq_path,sep="\t",header=0)
ch_atac_data = pybedtools.BedTool(atac_seq_path)
#merged_atac=ch_atac_data.merge()
#%%
# read screen annotations for chondrocytes
#screen_file_hg19 = r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\GenomeAnnotation\Human\SCREEN\V3\GRCh37-lifted-cCREs.bed"
screen_file_hg19 = '/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/SCREEN/V3/GRCh37-lifted-cCREs.bed'
#screen_annotation_ch_file = r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\psychMPRA\osteoblast_hg38_v3.bed"
screen_annotation_ch_file = '/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/psychMPRA/osteoblast_hg38_v3.bed'
# the annotations were truncated in the liftover, need to join with the hg38 file
#%%
screen_hg19_df = pd.read_csv(screen_file_hg19, sep='\t', header=None)
screen_annotation_ch_df = pd.read_csv(screen_annotation_ch_file, sep='\t', header=None)
#%%
screen_annotation_ch_df.drop(columns=[0, 1, 2, 4], inplace=True)
#%%
screen_df = pd.merge(screen_hg19_df, screen_annotation_ch_df, left_on=4, right_on=3)
#%%
screen_df.drop(columns=['5_x','3_y','5_y',6,7,8], inplace=True)
screen_df.rename(columns={9: "regulatoryClass"}, inplace=True)
#%%
# keep only regulatory regions. This is the list of classes:
# Classes are concatenated, only records where the only class is CTCF - should be removed
# PLS
# PLS,CTCF-bound
# pELS
# pELS,CTCF-bound
# dELS
# dELS,CTCF-bound
# DNase-H3K4me3
# DNase-H3K4me3,CTCF-bound
# CTCF-only,CTCF-bound
vals=['CTCF-only,CTCF-bound','Low-DNase','DNase-only']
screen_df_filtered=screen_df[screen_df['regulatoryClass'].isin(values=vals)==False]
#%%
screen_df_filtered = screen_df_filtered.sort_values([0, 1], ascending=[True, True])
screen_df_filtered['regulatoryClass'] = screen_df_filtered['regulatoryClass'].str.replace(',CTCF-bound', '')
#%%
screen_df_filtered.rename(columns={0: 'chr', 1: 'start', 2: 'end', '3_x': 'ID2', 4: 'ID1',10:"dtype"}, inplace=True)
screen_bed=pybedtools.BedTool.from_dataframe(screen_df_filtered)
#%%
# intersect atac seq and screen to get overlapping oligos
intersect_bed = ch_atac_data.intersect(screen_bed,wo=True)
#intersect_bed = merged_atac.intersect(screen_bed,wo=True)
intersect_df = intersect_bed.to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart',
                                                 'thickEnd', 'itemRgb', 'screen_chr','screen_start','screen_end', 'ID2',
                                                 'ID1',"regulatoryClass","dtype","overlap"])
#intersect_df = intersect_bed.to_dataframe(names=['chrom', 'start', 'end', 'screen_chr','screen_start','screen_end', 'ID2',
#                                                 'ID1',"regulatoryClass","dtype","overlap"])
intersect_df.sort_values(by='overlap',axis=0,ascending=False,inplace=True)
#intersect_df=intersect_df.drop_duplicates(subset=['ID1'])
intersect_df.drop_duplicates(subset=['chrom','start','end'],inplace=True)
intersect_df['peak length']=intersect_df['end']-intersect_df['start']
intersect_df['screen length']=intersect_df['screen_end']-intersect_df['screen_start']
intersect_df['peak coverage']=intersect_df['overlap']/intersect_df['peak length']
intersect_df['screen coverage']=intersect_df['overlap']/intersect_df['screen length']

#%%

#%%
intersect_df.to_csv("/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/psychMPRA/"
                    "chondrocytes_atac_seq_filtered_screen.csv")
#%%
intersect_df=pd.read_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly"
                         r"\psychMPRA\chondrocytes_atac_seq_filtered_screen.csv",index_col=0)
#switch to 1-base
intersect_df['start']=intersect_df['start']+1
#%%
filtered_ch_peaks_df=intersect_df[(intersect_df['peak coverage']==1) | ((intersect_df['screen coverage']==1) & (intersect_df['peak coverage']>=0.85))].copy()
#filtered_ch_peaks_df=intersect_df[intersect_df['peak coverage']==1]
filtered_ch_peaks_df['ATAC_seq_ID']=filtered_ch_peaks_df.name.apply(lambda x: re.match(r"^ID=(SRX\d+)",x)[1])
#%%
filtered_ch_peaks_df['old_name']=filtered_ch_peaks_df['name']
filtered_ch_peaks_df['name']=["ATAC-seq-chond."+str(i) for i in range(len(filtered_ch_peaks_df))]

#%%

#atac_seq_df_tiles=df_tiling(intersect_df,50,270)
atac_seq_df_tiles=df_tiling(filtered_ch_peaks_df,135,270)
#%%
final_df_atac_ch=pd.DataFrame(data={"chr":atac_seq_df_tiles['chromosome'], "start":atac_seq_df_tiles['oligo_start'],
                                   "end":atac_seq_df_tiles['oligo_end'],
                                   "group":"Random_chond","name":atac_seq_df_tiles['name'],
                                    'external_ID':atac_seq_df_tiles['ID']})
#%%
'''Neurons Analysis, run on server'''
# read screen annotations - across all cell types
#screen_annotation_file = r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\psychMPRA\osteoblast_hg38_v3.bed"
screen_annotation_file = '/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/SCREEN/V3/GRCh38-cCREs.bed'
#%%
screen_annotation_df = pd.read_csv(screen_annotation_file, sep='\t', header=None)
#%%
screen_annotation_df.drop(columns=[0, 1, 2, 4], inplace=True)
#%%
screen_df = pd.merge(screen_hg19_df, screen_annotation_df, left_on=3, right_on=3)
#%%
screen_df.drop(columns=['5_x'], inplace=True)
screen_df.rename(columns={0: 'chr', 1: 'start', 2: 'end', 3: 'ID2', 4: 'ID1','5_y':'regulatoryClass'}
                 , inplace=True)
vals=['CTCF-only,CTCF-bound','Low-DNase','DNase-only']
screen_df_filtered=screen_df[screen_df['regulatoryClass'].isin(values=vals)==False]
screen_df_filtered = screen_df_filtered.sort_values(['chr', 'start'], ascending=[True, True])
screen_df_filtered['regulatoryClass'] = screen_df_filtered['regulatoryClass'].str.replace(',CTCF-bound', '')
screen_bed=pybedtools.BedTool.from_dataframe(screen_df_filtered)
active_df=active_df.sort_values(['insert_chrom', 'insert_start'], ascending=[True, True])
active_bed=pybedtools.BedTool.from_dataframe(active_df)
#active_bed_merged=active_bed.merge()
# intersect active ataq seq data in neurons with screen data to get overlapping oligos
intersect_psych_bed = active_bed.intersect(screen_bed,wo=True)
intersect_psych_df = intersect_psych_bed.to_dataframe(names=range(59))
intersect_psych_df.sort_values(by=58,axis=0,ascending=False,inplace=True)
intersect_psych_df.drop_duplicates(subset=[0,1,2],inplace=True)
intersect_psych_df_short=intersect_psych_df[[0,1,2,3,4,52,53,54,55,56,57,58]].copy()
intersect_psych_df_short.columns=['chrom','start','end',"name",'RNA/DNA_ratio','screen_chr','screen_start','screen_end'
    ,"ID2","ID1","regulatoryClass","overlap"]
intersect_psych_df.to_csv("/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/psychMPRA/"
                    "psychMPRA_filtered_oligos.csv")
intersect_psych_df_short['peak length']=intersect_psych_df_short['end']-intersect_psych_df_short['start']
intersect_psych_df_short['screen length']=intersect_psych_df_short['screen_end']-intersect_psych_df_short['screen_start']
intersect_psych_df_short['peak coverage']=intersect_psych_df_short['overlap']/intersect_psych_df_short['peak length']
intersect_psych_df_short['screen coverage']=intersect_psych_df_short['overlap']/intersect_psych_df_short['screen length']
intersect_psych_df_short.to_csv("/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/psychMPRA/"
                    "psychMPRA_filtered_oligos_short.csv")
#%%
intersect_psych_df=pd.read_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly"
                               r"\psychMPRA\psychMPRA_filtered_oligos_short.csv",index_col=0)
#%%
intersect_psych_df.drop_duplicates(subset=['ID1'],inplace=True)
#%%
intersect_psych_df.rename(columns={"0":'chr',"1":'start',"2":'end',"3":"name","4":"ratio"},inplace=True)
#%%
filtered_psych_df=intersect_psych_df[(intersect_psych_df['peak coverage']==1) | (intersect_psych_df['screen coverage']==1)].copy()
#filtered_psych_df=intersect_psych_df[(intersect_psych_df['peak coverage']==1) | (intersect_psych_df['screen coverage']==1)].copy()

#%%
filtered_psych_df['start']=filtered_psych_df['start']+1
filtered_psych_df['old_name']=filtered_psych_df['name']
filtered_psych_df['name']=["PsychMPRA."+str(i) for i in range(len(filtered_psych_df))]
#%%
psych_tiles=df_tiling(filtered_psych_df,135,270)

#%%
final_df_psych=pd.DataFrame(data={"chr":psych_tiles['chromosome'], "start":psych_tiles['oligo_start'],
                                   "end":psych_tiles['oligo_end'],
                                   "group":"random_neurons","name":psych_tiles['name'],'external_ID':psych_tiles['ID']})

#%%
final_random_df=pd.concat([final_df_atac_ch,final_df_psych],axis=0)
#%%
final_random_df['sequence_hg19']=final_random_df.apply(func=lambda x: seq_retrieve(x['chr'],x['start'],x['end'],records),axis=1)

#%%
final_random_df.to_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis"
                      r"\main\Library_assembly\output\pre_satmut_random_active.csv")
#%%
'''
screen_file_hg19 = r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\GenomeAnnotation\Human\SCREEN\V3\GRCh37-lifted-cCREs.bed"
# screen_file_hg19 = '/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/SCREEN/V3/GRCh37-lifted-cCREs.bed'
screen_annotation_file = r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\GenomeAnnotation\Human\SCREEN\V3\GRCh38-cCREs.bed"
# screen_annotation_file = '/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/SCREEN/V3/GRCh38-cCREs.bed'
screen_hg19_df = pd.read_csv(screen_file_hg19, sep='\t', header=None)
screen_annotation_df = pd.read_csv(screen_annotation_file, sep='\t', header=None)
#%%
screen_annotation_df.drop(columns=[0, 1, 2, 4], inplace=True)
#%%
screen_df = pd.merge(screen_hg19_df, screen_annotation_df, left_on=3, right_on=3)
#%%
screen_df.drop(columns=['5_x'], inplace=True)
screen_df.columns=['chr','start','end',"ID2","ID1","regulatoryClass"]
#%%
vals=['CTCF-only,CTCF-bound','Low-DNase','DNase-only']
screen_df_filtered=screen_df[screen_df['regulatoryClass'].isin(values=vals)==False]
screen_df_filtered = screen_df_filtered.sort_values(['chr', 'start'], ascending=[True, True])
screen_df_filtered['regulatoryClass'] = screen_df_filtered['regulatoryClass'].str.replace(',CTCF-bound', '')
screen_df_filtered['Length']=screen_df_filtered['end']-screen_df_filtered['start']
'''