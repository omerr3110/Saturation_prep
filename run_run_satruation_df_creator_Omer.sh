#!/bin/bash

param_file=$1
oligo_file=$2
file_name=$3
#param_file=AF_STATS_config_Omer.txt
folder=$(awk -F'\t' '$1 == "output vcf path"{print$2}' $param_file)
#mkdir $folder
#mkdir $folder"/log/"
#mkdir $folder"/log/run_saturation_df_creator_Omer"
echo $folder
echo $oligo_file
echo "$param_file"


bsub -q molgen-q -R "rusage[mem=32000]" -e $folder"/log/run_saturation_df_creator_Omer/run_saturation_df_creator_Omer."$file_name"_%J.e.txt" -o $folder"/log/run_saturation_df_creator_Omer/run_saturation_df_creator_Omer."$file_name"_%J.o.txt" ./run_saturation_df_creator_Omer.sh $oligo_file
#bsub -q molgen-q -e $folder"/log/run_Ancestral_allele_Omer/run_Ancestral_allele_Omer.2.19_%J.e.txt" -o $folder"/log/run_Ancestral_allele_Omer/run_Ancestral_allele_Omer.2.19_%J.o.txt" ./run_Ancestral_allele_Omer.sh 2 19
#replace molgen next time


