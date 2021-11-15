#!/bin/bash
#$ -cwd
#$ -l s_core=1
#$ -l h_rt=0:10:00
#$ -j y

target_pdb_dir=${1}
output_csv_path=${2}

source /home/4/16B09097/.bashrc
conda activate mypython

python DOPE.py ${target_pdb_dir} ${output_csv_path}
