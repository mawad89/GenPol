#!/bin/sh -login

software=/home/awad/python/v_c/variant_calling_software/31_1_2020/run_5_4_2020_only_for_hirsuta .py
bed=
ref=
draft=
output=$(pwd)
output_file_name=variant_calling 
mapping_quality=30
job_name=variant_calling
submit=True
hits=
gather_name=
python $software -b $bed -r $ref -g $draft -o $output -f $output_file_name -q $mapping_quality -n $job_name -t $submit -a $hits -c $gather_name

