#!/bin/bash
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=3:00:00
#$ -j y

. /etc/profile.d/modules.sh
module purge
module load cuda/11.2.146
module li

# args
alphafold_output_dir=$1 # alphafold_output/target
native_pdb=$2 # native_pdb/target.pdb
output_score_path=$3 # score/label_relaxed/target.csv
output_pdb_dir=${alphafold_output_dir}/relaxed

# decomprss model_pickle.tar.gz
pushd ${alphafold_output_dir}
tar -xzf model_pickle.tar.gz
popd

# run relax
echo "Running relax..."
python_path="../../alphafold/colabfold-conda/bin/python3.7" # Set python path for alphafold
${python_path} relax_prediction.py -p ${alphafold_output_dir} -o ${output_pdb_dir}

# Calculate accuracy
echo "Calculating accuracy..."
source /home/4/16B09097/.bashrc
module load python/3.9.2
source ../../../../.venv/bin/activate
python calc_relaxed_accuracy.py -i ${output_pdb_dir} -n ${native_pdb} -o ${output_score_path}


# remove model_*.pickle in ${alphafold_output_dir}
rm -f ${alphafold_output_dir}/model_*.pickle

echo "Finished"
