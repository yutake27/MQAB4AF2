#!/bin/bash

# For TSUBAME scheduling system
#$ -cwd
#$ -l q_node=1
#$ -j y
# The execution time is specified when submitting a job according to the sequence length.


date
# Init (For TSUBAME)
. /etc/profile.d/modules.sh
module purge
module load cuda/11.2.146
module li

# Set Alphafold execution file directory
# alphafold_path="SET YOUR OWN PATH"
alphafold_path=alphafold

# Args
fasta=${1}  # absolute path
output=${2} # absolute path

# precomputed MSA path if exists
msa_pickle_path=${output}/msa.pickle

# parameters
max_recycles="10"
num_samples="2"
model="both"
num_CASP14_models="5"
num_pTM_models="5"
num_ensembles="1 8"

# Run colabfold
date
cd ${alphafold_path}
if [ ! -e ${msa_pickle_path} ]; then
    echo precomputed MSA does not exist
    colabfold-conda/bin/python3.7 runner_af2advanced.py \
    -i ${fasta} \
    -o ${output} \
    --max_recycles ${max_recycles} \
    --num_samples ${num_samples} \
    --model ${model} \
    --num_CASP14_models ${num_CASP14_models} \
    --num_pTM_models ${num_pTM_models} \
    --noranking \
    --output_all_cycle \
    --num_ensembles ${num_ensembles}
else
    echo precomputed MSA exists
    colabfold-conda/bin/python3.7 runner_af2advanced.py \
    -i ${fasta} \
    -o ${output} \
    --max_recycles ${max_recycles} \
    --num_samples ${num_samples} \
    --model ${model} \
    --num_CASP14_models ${num_CASP14_models} \
    --num_pTM_models ${num_pTM_models} \
    --noranking \
    --output_all_cycle \
    --num_ensembles ${num_ensembles} \
    --msa_method precomputed \
    --precomputed ${msa_pickle_path}
fi
date