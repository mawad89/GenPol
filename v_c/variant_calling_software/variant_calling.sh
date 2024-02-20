#!/bin/sh -login

software=/home/awad/python/v_c/variant calling_software/run.py
bed=
ref=/biodata/dep_tsiantis/common/chigenome/chi_v1.fa
draft=
output=$(pwd)
output_file_name=variant_calling 
mapping_quality=30
job_name=variant_calling
submit=True
python $software -b $bed -r $ref -g $draft -o $output -f $output_file_name -q $mapping_quality -n $job_name -t $submit

