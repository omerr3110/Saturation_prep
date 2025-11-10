#%%
import numpy as np
import pandas as pd
import Bio.Seq
import re
import pickle
#%%
def seq_retrieve(chr,start,end,records_file):
    """
    always enter with 1-based
    :param chr:
    :param start:
    :param end:
    :param records_file:
    :return:
    """
    curr_chr=[record for record in records_file if record.id==chr]
    return curr_chr[0].seq[start-1:end].upper()

def df_tiling(df,step_size,oligo_length):
    tiles_df = pd.DataFrame()
    df_list = []
    for i, row in enumerate(df.itertuples()):
        #print(row)
        curr_length = abs(row[3] - row[2])
        #print(curr_length)
        if curr_length > 270:
            tiles_number = np.ceil(curr_length / step_size) + 1
            tiling_region_length = (tiles_number - 1) * step_size
        else:
            tiles_number = 1
            tiling_region_length = 270
        midpoint = np.ceil((row[3] + row[2]) / 2)
        tiling_start = midpoint - tiling_region_length / 2
        tiling_end = midpoint + tiling_region_length / 2 - 1
        # create tiles
        for tile in np.arange(1, tiles_number + 1):
            curr_tile_dict = dict()
            curr_tile_dict['chromosome'] = row[1]
            if tiles_number > 1:
                oligo_midpoint = tiling_start + (tile - 1) * step_size
                curr_tile_dict['oligo_start'] = oligo_midpoint - oligo_length / 2
                curr_tile_dict['oligo_end'] = oligo_midpoint + oligo_length / 2 - 1
                curr_tile_dict['tileID'] = tile
            else:
                oligo_midpoint = midpoint
                curr_tile_dict['oligo_start'] = oligo_midpoint - oligo_length / 2
                curr_tile_dict['oligo_end'] = oligo_midpoint + oligo_length / 2 - 1
                curr_tile_dict['tileID'] = tile
            curr_tile_dict['name'] = str(row[4]) + '_tile' + str(int(tile))
            curr_tile_dict['ID'] = str(row[4])
            curr_tile_df = pd.DataFrame.from_dict([curr_tile_dict])
            df_list.append(curr_tile_df)
            # tiles_df = pd.concat([tiles_df, curr_tile_df], axis=0)
    tiles_df = pd.concat(df_list, axis=0)
    tiles_df.reset_index(inplace=True)
    tiles_df.drop(columns=['index'], inplace=True)
    tiles_df = tiles_df.astype({'oligo_start': 'int64', 'oligo_end': 'int64'})
    return tiles_df


def anc_retreive(row):
    seq=row['sequence_hg19']
    chrom=re.findall("\d+|X|Y",row['chr'])[0]
    start=row['start']
    #print(type(start))
    end=row['end']
    portion_st=start//(5*10**6)+1
    portion_end=end//(5*10**6)+1
    #print(chrom,portion_end,portion_st,start,end)
    if portion_st!=portion_end:
        print("Multiple Portions")
        print(row['Name'],chrom,portion_st,portion_end)
    anc_path=rf"\\data.wexac.weizmann.ac.il\davidgo\omerro\SaturationMutagenesis\main\Library_assembly\output\{chrom}\{portion_st}\major_{chrom}_{portion_st}_full_with_anc.txt"
    try:
        anc_df = pd.read_csv(anc_path,sep="\t",index_col=0)
    except FileNotFoundError:
        return seq
    anc_df['pos']=anc_df['pos'].apply(lambda x: int(x))
    anc_df=anc_df.query('pos<=@end and pos>=@start')
    if len(anc_df)<=0:
        return seq
    seq_list=list(seq)
    #print("here")
    for i,nuc in enumerate(seq_list):
        #print(i+start)
        pos_df=anc_df.query('pos==@i+@start')
        #print(pos_df)
        if len(pos_df)>0:
            anc=pos_df['HC_anc'].values[0]
            if anc!=nuc:
                print("here")
                print(anc,nuc)
                seq_list[i]=anc
    new_seq=Bio.Seq.Seq("".join(seq_list))
    return new_seq


#%%
def position_calc(var_dict,chrom,start,end):
    new_var_dict={k: v for k, v in var_dict.items() if k[0] == chrom and end >= k[1] >= start}
    norm_positions= {k[1]-start: v for k,v in new_var_dict.items()}
    filtered_norm_positions={k:v for k,v in norm_positions.items() if k<35 or k>=235}
    return filtered_norm_positions
