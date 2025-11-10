import pandas as pd
#import pybedtools
import pickle
#%%
'''run on server - creating a Human derived variants dictionary'''
variant_bed=pybedtools.BedTool(r"/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_6_CGI_catalog"
                              r"/FILTER_INCLUSIVE_STRICT_0_3_good_samples3_delta/INDEL_filtered"
                              r"/new_flanking_catalog_unique.bed")
variant_df=variant_bed.to_dataframe()
variant_df['variant']=variant_df['name'].apply(lambda x: x.split("->")[1])
variant_dict=dict(zip(tuple(zip(variant_df.chrom,variant_df['end']))
                  ,variant_df.variant))

with open(r"/home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_6_CGI_catalog"
                              r"/FILTER_INCLUSIVE_STRICT_0_3_good_samples3_delta/INDEL_filtered"
                              r"/Human_derived_dict.pkl", 'wb') as f:
    pickle.dump(variant_dict, f)
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
print(len(merged_dict))
print(len(dict1)+len(dict2)+len(dict3))