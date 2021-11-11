#!/bin/bash

# For TSUBAME scheduling system
#$ -cwd
#$ -l q_node=1
#$ -j y
# The execution time is specified when submitting a job according to the sequence length.

# Parse the command line arguments
usage() {
        echo ""
        echo "Usage: $0 -f fasta_path -o output_dir <OPTIONS>"
        echo "Required Parametes:"
        echo "-f <fasta_path> Absolute fasta path"
        echo "-o <output_dir> Absolute output directory"
        echo "Optional Parameters:"
        echo "-m <model> Models to run. Choices=[CASP14, pTM, both]. Default is both"
        echo "-c <num_CASP14_models> Number of CASP14 models to run. Defaults is 5."
        echo "-p <num_pTM_models> Number of pTM models to run. Defaults is 5."
        echo "-e <num_ensemble> Number of ensembles. Choices=['1', '8', '1 8']. Default is '1 8'."
        echo "-s <num_samples> Number of random seeeds. Defaults is 2."
        echo "-r <max_recycles> Number of max recycles. Defaults is 10."
        echo "-a <num_relax> Number of relax using amber. Choices=['None', 'Top1', 'Top5', 'All'], Defaults is None."
        echo ""
        exit 1
}

while getopts ":f:o:m:c:p:e:s:r:a:h" i; do
        case "${i}" in
        f)
                fasta_path=$OPTARG
        ;;
        o)
                output_dir=$OPTARG
        ;;
        m)
                model=$OPTARG
        ;;
        c)
                num_CASP14_models=$OPTARG
        ;;
        p)
                num_pTM_models=$OPTARG
        ;;
        e)
                num_ensembles=$OPTARG
        ;;
        s)
                num_samples=$OPTARG
        ;;
        r)
                max_recycles=$OPTARG
        ;;
        a)
                num_relax=$OPTARG
        ;;
        h) usage
        ;;
        \?) usage
        ;;
        esac
done

# Parse input and set defaults
if [[ "${fasta_path}" == "" || "${output_dir}" == "" ]]; then
        echo "ERROR: Missing required parameters"
        usage
fi

if [[ "${model}" == "" ]]; then
        model="both"
elif [[ "${model}" != "CASP14" && "${model}" != "pTM" && "${model}" != "both" ]]; then
        usage
fi

if [[ "${num_CASP14_models}" == "" ]]; then
        num_CASP14_models=5
fi

if [[ "${num_pTM_models}" == "" ]]; then
        num_pTM_models=5
fi

if [[ "${num_ensembles}" == "" ]]; then
        num_ensembles="1 8"
fi

if [[ "${num_samples}" == "" ]]; then
        num_samples=2
fi

if [[ "${max_recycles}" == "" ]]; then
        max_recycles=10
fi

if [[ "${num_relax}" == "" ]]; then
        num_relax="None"
fi


date
# Init (For TSUBAME)
. /etc/profile.d/modules.sh
module purge
module load cuda/11.2.146
module li

# Set Alphafold execution file directory
# alphafold_path="SET YOUR OWN PATH"
alphafold_path=alphafold

# precomputed MSA path if exists
msa_pickle_path=${output_dir}/msa.pickle

# Run colabfold
date
cmd="\
colabfold-conda/bin/python3.7 runner_af2advanced.py \
-i ${fasta_path} \
-o ${output_dir} \
--model ${model} \
--num_CASP14_models ${num_CASP14_models} \
--num_pTM_models ${num_pTM_models} \
--noranking \
--max_recycles ${max_recycles} \
--output_all_cycle \
--num_ensembles ${num_ensembles} \
--num_samples ${num_samples} \
--num_relax ${num_relax} \
"
pushd ${alphafold_path}
if [ ! -e ${msa_pickle_path} ]; then
    echo precomputed MSA does not exist
else
    echo precomputed MSA exists
    cmd+="--msa_method precomputed --precomputed ${msa_pickle_path}"
fi
echo ${cmd}
mkdir -p ${output_dir}
echo "`date` ${cmd}" >> ${output_dir}/colabfold.log
${cmd} # Execute colabfold
date
popd

# post processing
if [ ! -e ${output_score_path} ]; then
    echo "AlphaFold prediction is not complete."
    exit 1
fi

module load python/3.9.2
source ../../.venv/bin/activate
python post_process.py ${output_dir}
date