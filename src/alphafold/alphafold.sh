#!/bin/bash

# For TSUBAME scheduling system
#$ -cwd
#$ -l q_node=1
#$ -j y
# The execution time is specified when submitting a job according to the sequence length.

# Init (For TSUBAME)
. /etc/profile.d/modules.sh
module purge
module load cuda/11.0.3 alphafold/2.0.0
module li

# Set Alphafold execution file directory
# alphafold_path="SET YOUR OWN PATH"
alphafold_path=/gs/hs0/tga-ishidalab/takei/bin/alphafold

# Args
fasta=${1}
output=${2}
date=${3}

# Convert relative path to absolute path
fasta=`readlink -f ${fasta}`
output=`readlink -f ${output}`

# Execution (non-docker)
date
cd $alphafold_path
# Use a single GPU
./run_alphafold.sh -g -d ${ALPHAFOLD_DATA_DIR} -o ${output} -a 0 -f ${fasta} -t ${date} \
    -m model_1,model_2,model_3,model_4,model_5,model_1_ptm,model_2_ptm,model_3_ptm,model_4_ptm,model_5_ptm \

date