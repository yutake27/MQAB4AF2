#!/bin/bash
#$ -cwd
#$ -l s_core=1
#$ -l h_rt=1:30:00
#$ -j y

target=${1}
length=${2}

source /home/4/16B09097/.bashrc
python get_global_score_resolved.py ${target} ${length}