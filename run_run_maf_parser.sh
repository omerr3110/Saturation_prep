#!/bin/bash
# make sure you're in /home/labs/davidgo/Collaboration/Tomas_Apes_Catalog_parsing_5 and ml python/2.7

param_file=$1
#param_file=AF_STATS_config_Omer.txt
folder=$(awk -F'\t' '$1 == "output vcf path"{print$2}' $param_file)
#mkdir $folder
#mkdir $folder"/log/"
#mkdir $folder"/log/run_maf_parser"


for chr in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}
do
echo "$chr"
for portion in {1..52}
do

echo "$chr $portion"
bsub -q molgen-q -e $folder"/log/run_maf_parser/run_maf_parser."$chr"."$portion"_%J.e.txt" -o $folder"/log/run_maf_parser/run_maf_parser."$chr"."$portion"_%J.o.txt" ./run_maf_parser.sh $chr $portion
#bsub -q molgen-q -e $folder"/log/run_maf_parser/run_maf_parser.22.4_%J.e.txt" -o $folder"/log/run_maf_parser/run_maf_parser.22.4_%J.o.txt" ./run_maf_parser.sh 22 4

done
done



