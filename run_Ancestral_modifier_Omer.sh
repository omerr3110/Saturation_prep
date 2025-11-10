#!/bin/bash
ml python/3.7.0
#module purge
#module load python/2.7
echo "entered shell"
python "./create_ancestral_df.py" "$1"
