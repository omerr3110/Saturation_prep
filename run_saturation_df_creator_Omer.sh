#!/bin/bash

#module purge
#module load python/2.7
echo "entered shell"
python "./create_saturation_df.py" "$1" 
