#!/bin/bash
# make sure you're in /home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5 and ml python/2.7

param_file=$1
oligo_file=$2
file_name=$3
#param_file=AF_STATS_config_Omer.txt
folder=$(awk -F'\t' '$1 == "output vcf path"{print$2}' $param_file)
#mkdir $folder
#mkdir $folder"/log/"
#mkdir $folder"/log/run_Ancestral_modifier_Omer"
echo $folder
echo $oligo_file
echo "$param_file"


bsub -q molgen-q -e $folder"/log/run_Ancestral_modifier_Omer/run_Ancestral_modifier_Omer."$file_name"_%J.e.txt" -o $folder"/log/run_Ancestral_modifier_Omer/run_Ancestral_modifier_Omer."$file_name"_%J.o.txt" ./run_Ancestral_modifier_Omer.sh $oligo_file
#bsub -q molgen-q -e $folder"/log/run_Ancestral_allele_Omer/run_Ancestral_allele_Omer.2.19_%J.e.txt" -o $folder"/log/run_Ancestral_allele_Omer/run_Ancestral_allele_Omer.2.19_%J.o.txt" ./run_Ancestral_allele_Omer.sh 2 19
#replace molgen next time


