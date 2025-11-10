#%%
from Bio import AlignIO
import pandas as pd
import re
import numpy as np
import sys
import time

#msa_gen=AlignIO.parse(r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\GenomeAnnotation\Human\UCSC_maf\chr22.maf", "maf")
#catalog_df=pd.read_csv(r"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\output\22\4\major_22_4_full_with_anc_HCG.txt",sep="\t",index_col=0)
#for rec in msa_gen:
#    if not re.search("hg19",rec[0].name):
#        print(rec)

#AlignIO.MafIO.MafIndex(r"\\data.wexac.weizmann.ac.il\davidgo\Collaboration\GenomeAnnotation\Human\UCSC_maf\chr22.maf")

#results = idx.search([3010000], [3013000])

#df=pd.DataFrame(data={"n":[5,5],"c":['a','b']})

def gen_test(n,path):
    ind=[]
    generator=AlignIO.parse(path, "maf")
    for i,inst in enumerate(generator):
        if i>=n:
            return ind
        else:
            ind.append(inst[0].annotations['start'])



def gibb_search(positions,record):
    """

    :param positions: series of the positions that are in the catalog
    :param record: current MSA from the MAF file
    :return: a list of the alleles corresponding to the catalog positions
    """
    #print(pos)
    #pos=int(pos)
    rec_start=record[0].annotations['start']+1  # file is 0 based - switch to 1 based to match the catalog
    rec_len=record.get_alignment_length()  # size of alignment (including None)
    flag=0
    gibb_record=None
    for spec_rec in record:  # search for the gibbon record
            if re.search("nomLeu3",spec_rec.name):
                gibb_record=spec_rec
                flag=1
                break
    if flag==0:  # if the gibbon record was missing in the MSA
        length=catalog_df.query('@rec_start<=pos<@end+1').shape[0]
        #print("missing")
        return ["M" for _ in range(length)]  # M for missing
    gaps = [1 if record[0].seq[i]=="-" else 0 for i in range(rec_len)]
    sum_gaps=np.cumsum(gaps)
    indices=[i if record[0].seq[i]!="-" else np.inf for i in range(rec_len)]
    real_indices=np.array(indices)-np.array(sum_gaps)+rec_start
    #print(real_indices)
    bool_vec=np.isin(real_indices,positions)
    gib_nucs=list(filter(None,[gibb_record.seq[i] if bool_vec[i] else 0 for i in range(len(bool_vec))]))
    return gib_nucs


#sample_df=catalog_df.iloc[0:2].copy()
chrom=sys.argv[1]
portion=sys.argv[2]
file_path=fr"/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/UCSC_maf/chr{chrom}.maf"
generator = AlignIO.parse(file_path, "maf")
catalog_path=fr"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/{chrom}/{portion}/major_{chrom}_{portion}_full_with_anc_HCG.txt"
try:
    catalog_df=pd.read_csv(catalog_path,sep="\t",index_col=0,na_values=[None])
    catalog_df.sort_values(by="pos",inplace=True)
    #print(catalog_df.head())
    catalog_df=catalog_df.replace({np.nan:None})
    nuc_array = ["" for _ in range(catalog_df.shape[0])]
    lower_bound = catalog_df.pos.min()
    upper_bound = catalog_df.pos.max()
    curr_index = 0
    new_catalog_path=fr"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/{chrom}/{portion}/major_{chrom}_{portion}_full_with_gibbon.txt"
    for i, msa in enumerate(generator):  # iterate on each generator entry
        start = msa[0].annotations['start'] + 1
        end = msa[0].annotations['start'] + msa[0].annotations['size']
        if end < lower_bound:  # we haven't reached our interval yet
            continue
        if start > upper_bound: # we already passed our interval - can stop iterating
            break
        # start_time = time.time()
        result = gibb_search(catalog_df.pos, msa)
        extend = len(result)
        if catalog_df.query('@start<=pos<@end+1').shape[0] != len(result):  # sanity check - test matching lengths
            print(extend, catalog_df.query('@start<=pos<@end+1').shape[0])
        nuc_array[curr_index:extend + curr_index] = result  # insert corresponding alleles to the final array
        curr_index += extend
    catalog_df['gibbon'] = nuc_array
    catalog_df['gibbon']=catalog_df['gibbon'].apply(str.upper)
    catalog_df.to_csv(new_catalog_path,sep="\t")
    catalog_df[['pos','gibbon']].to_csv(
        fr"/home/labs/davidgo/omerro/SaturationMutagenesis/main/Library_assembly/output/{chrom}/{portion}/major_{chrom}_{portion}_gibbon.txt",
        sep="\t")
except FileNotFoundError:
    print("No file")
except pd.errors.EmptyDataError:
    print("empty file")
