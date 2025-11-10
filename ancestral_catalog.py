import io
import os
import pandas as pd
import vcf
import numpy as np
import time
import sys
sys.path.append("/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5")
import quality_functions

# this script is based on the major ape catalog script FILTER_AF_1_advanced.py
# it is modified as follows:
# It runs on the entire catalog and counts for each position the number of non-reference alleles of the major non reference allele, and the number of non reference alleles
# create a df for each chromosome portion - major allele in Human, Chimp, Gorilla, Pongo and hg19
# python ./../../simonf/backup/MonogenicDiseases/scripts/ApeCatalogAF.py "1" "1"


def row_df_creator(rec):
    curr_group_df = group_df.copy()
    # switching to using a dictionary - better efficiency
    group_dict = dict()
    for group, group_data in group_df.iterrows():
        group_dict[group] = dict()
        group_dict[group]['ALT_count'] = 0
        group_dict[group]['genotyped'] = 0
        group_dict[group]['alt_minus_1'] = 0
        group_dict[group]['alt_0'] = 0
        group_dict[group]['alt_1'] = 0
        group_dict[group]['alt_2'] = 0
        group_dict[group]['alt_3'] = 0
    for sample in rec.samples:
        if sample.sample not in bad_samples_df.index:
            # identify the group of the sample
            curr_group = samples_dict2[sample.sample]
            # populate the curr_group_df
            alt_allele_count, genotyped_allele_count, valid_GT_arr = quality_functions.quality_function(sample, params_dict,DP_stats_df)
            group_dict[curr_group]['ALT_count'] += alt_allele_count
            group_dict[curr_group]['genotyped'] += genotyped_allele_count
            group_dict[curr_group]["alt_" + str(valid_GT_arr[0]).replace("-", "minus_")] += 1
            group_dict[curr_group]["alt_" + str(valid_GT_arr[1]).replace("-", "minus_")] += 1
    # placing the dictionary data in the df
    for temp_group, group_data in curr_group_df.iterrows():
        curr_group_df.at[temp_group, 'ALT_count'] = group_dict[temp_group]['ALT_count']
        curr_group_df.at[temp_group, 'genotyped'] = group_dict[temp_group]['genotyped']
        curr_group_df.at[temp_group, 'alt_minus_1'] = group_dict[temp_group]['alt_minus_1']
        curr_group_df.at[temp_group, 'alt_0'] = group_dict[temp_group]['alt_0']
        curr_group_df.at[temp_group, 'alt_1'] = group_dict[temp_group]['alt_1']
        curr_group_df.at[temp_group, 'alt_2'] = group_dict[temp_group]['alt_2']
        curr_group_df.at[temp_group, 'alt_3'] = group_dict[temp_group]['alt_3']
    return curr_group_df


def major_calc(row,rec):
    """

    :param row: an entry of the df - allele counts
    :param rec: the catalog record
    :return: major allele
    """
    #print(row)
    col=row.loc[['alt_0','alt_1','alt_2','alt_3']].argmax()
    ind=int(col.split("_")[1])
    if row.loc[col]==0:
        allele=None
    elif ind==0:
        allele=rec.REF
    else:
        allele=rec.ALT[ind-1]
    return allele


def table_builder(row_df,rec,gn_df):
    """

    :param row_df: counts df
    :param rec: the catalog record
    :param gn_df: matching gnomad df
    :return: dictionary with major alleles
    """
    chr_id=rec.CHROM
    position=rec.POS
    ref=rec.REF
    human_df=gn_df.query('pos==@position') ####check if can be improved
    if len(human_df)>0:
        ref_freq=1-human_df['MAF'].sum()
        alt_freq=human_df['MAF'].max()
        if ref_freq>=alt_freq:
            human_major=human_df['ref'].iloc[0]
        else:
            human_major=human_df[human_df['MAF']==human_df['MAF'].max()].alt.iloc[0]
    else:
        human_major=None
    majors=row_df.apply(func=major_calc,axis=1,rec=rec)
    chimp_major=majors[2]
    gorilla_major=majors[0]
    pongo_major=majors[3]
    row_dict={"chr":chr_id,"pos":position,"Human":human_major,"Chimp":chimp_major,"Gorilla":gorilla_major,
    "Pongo":pongo_major,"Hg19":ref}
    return row_dict



# get from cmd: 1. chromosome; 2. chr. portion 3. settings file name
print(sys.argv)
chromosome = sys.argv[1]
chromosome_portion = sys.argv[2]
filtering_params_file= sys.argv[3]

#filtering_params_file='/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5/AF_STATS_config_Omer.txt'

# import settings file
filtering_params_df = pd.read_csv(filtering_params_file, sep='\t')
filtering_params_df = filtering_params_df.set_index('parameter')

# params dictionary
params_dict={}
params_dict['DP_cutoff_bottom'] = int(filtering_params_df['value']['DP hard filtering bottom'])
params_dict['DP_cutoff_top'] = int(filtering_params_df['value']['DP hard filtering top'])
params_dict['DP_cutoff_per_sample'] = int(filtering_params_df['value']['DP 2.5% 7.5%'])
params_dict['GQ_cutoff'] = int(filtering_params_df['value']['GQ hard filtering bottom'])
params_dict['samples_metadata_file'] = filtering_params_df['value']['group file path']
params_dict['bad_samples_file'] = filtering_params_df['value']['bad sample names file path']
params_dict['vcf_folder_input'] = filtering_params_df['value']['input vcf path']
params_dict['vcf_path'] = params_dict['vcf_folder_input'] + '/' + chromosome + '_' + chromosome_portion + '_only_ALT.vcf.gz'
params_dict['DP_sample_stats_file'] = filtering_params_df['value']['DP sample stats']
params_dict['ALLELIC_BALANCE_HOMOZ'] = int(filtering_params_df['value']['allelic balance homozygous'])
params_dict['ALLELIC_BALANCE_HETERO'] = int(filtering_params_df['value']['allelic balance heterozygous'])

# params table
params_table = pd.DataFrame.from_dict(params_dict, "index")

print("entered main")

# import sample metadata and groups
samples_metadata_df = pd.read_csv(params_dict['samples_metadata_file'], sep='\t')
samples_metadata_df2 = samples_metadata_df.set_index('sample')
# having a dictionary for efficiency when querying
samples_dict = samples_metadata_df2.to_dict()
samples_dict2 = samples_dict['group']

# import dp stats per sample
DP_stats_df = pd.read_csv(params_dict['DP_sample_stats_file'], sep='\t')
DP_stats_df = DP_stats_df.set_index('sample')
# remove not used columns
DP_stats_df = DP_stats_df[["group","DP_2_5","DP_97_5"]]

# import bad sample names
bad_samples_df = pd.read_csv(params_dict['bad_samples_file'], sep='\t')
bad_samples_df = bad_samples_df.set_index('sampleName')

# remove bad samples from the group df
for sample in bad_samples_df.index:
    samples_metadata_df = samples_metadata_df.drop(samples_metadata_df[samples_metadata_df['sample'] == sample].index)

# a group dataframe
# this will be the basis for each variant filtering
group_df = samples_metadata_df.groupby(['group']).count()
#group_df = group_df.drop(group_df[group_df.index=='Homo sapiens'].index)

group_df['ALT_count'] = 0
group_df['genotyped'] = 0
group_df['alt_minus_1'] = 0
group_df['alt_0'] = 0
group_df['alt_1'] = 0
group_df['alt_2'] = 0
group_df['alt_3'] = 0
group_df['non_ref_AF'] = 0
group_df['MAF'] = 0

alt_column_list=['alt_1','alt_2','alt_3']

MAF_values = np.arange (0,159*2+1,1)

# initiate a dictionary that has for each possible AF count an array
MAF_hist_dict = {}
for bin in MAF_values:
    MAF_hist_dict[bin] = 0

non_ref_AF_hist_dict = {}
for bin in MAF_values:
    non_ref_AF_hist_dict[bin] = 0


if chromosome!="Y":
    gnomad_path="/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/Gnomad/hg19_v2.1.1/extracts/gnomad_extract_"+chromosome+"_full.txt"
    gnomad_df=pd.read_csv(gnomad_path,sep="\t",header=None, names=['chr','pos','name','ref','alt','score','bool','MAF'])
else:
    gnomad_df=pd.DataFrame(columns=['pos'])
# import vcf
vcf_reader = vcf.Reader(open(params_dict['vcf_path'], 'r'))
# start_time = time.time()
# k=0
#for record in vcf_reader:
#    if record.POS == 69682:
#        curr_record = record
#        print(record.POS)
#record = curr_record

dict_list=[]
df_list=[]
jump=10000

for i,record in enumerate(vcf_reader):
    if i%jump==0 and i>0:
        # every jump - save the df and start over
        print(i)
        curr_table=pd.DataFrame(dict_list)
        curr_table = curr_table[['chr', 'pos', 'Human',"Chimp","Gorilla",'Pongo','Hg19']]
        path="/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/{}/{}/major_{}_{}_{}.txt".format(chromosome,chromosome_portion,chromosome,chromosome_portion,str(i//jump))
        curr_table.to_csv(path,sep="\t")
        dict_list=[]
    #print(record)
    # set a temp group dataframe for this variant
    curr_row_df=row_df_creator(record)
    #print(record.POS)
    curr_dict=table_builder(curr_row_df,record,gnomad_df)
    dict_list.append(curr_dict)

print("done")
if len(dict_list)>0:
    curr_table=pd.DataFrame(dict_list)
    curr_table = curr_table[['chr', 'pos', 'Human',"Chimp","Gorilla",'Pongo','Hg19']]
    path="/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/{}/{}/major_{}_{}_{}.txt".format(chromosome,chromosome_portion,chromosome,chromosome_portion,str(i//jump+1))
    curr_table.to_csv(path,sep="\t")

