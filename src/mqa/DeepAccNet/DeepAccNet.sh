#!/bin/bash
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=3:00:00
#$ -j y

input_dir=${1}
output_csv_path=${2}
output_dir=${output_csv_path%.*}

source /home/4/16B09097/.bashrc
. /etc/profile.d/modules.sh
module load cuda/11.0.194 cudnn/7.6

DeepAccNet_dir=/gs/hs0/tga-ishidalab/takei/bin/DeepAccNet # set the directory of DeepAccNet
pushd ${DeepAccNet_dir}
source venv/bin/activate

mkdir -p ${output_dir}

# generate bert features
echo "generate bert features"
python extractBert.py ${input_dir} ${output_dir} --modelpath ProtBert-BFD

# DeepAccNet (standard)
echo "Running DeepAccNet (standard)"
python DeepAccNet.py -v -lt ${input_dir} ${output_dir}
popd
python mv_local_score.py ${output_dir}/DeepAccNet

# DeepAccNet-Bert
echo "Running DeepAccNet-Bert"
pushd ${DeepAccNet_dir}
python DeepAccNet.py -v --bert ${input_dir} ${output_dir}
popd ${exec_dir}
python mv_local_score.py ${output_dir}/DeepAccNet-Bert

# Get global scores
echo "Convert local scores to global scores"
python get_global_score.py ${output_dir}

# Compress local scores
echo "Compress local scores"
pushd ${output_dir}
tar -cvzf DeepAccNet.tar.gz DeepAccNet
rm -rf DeepAccNet
tar -cvzf DeepAccNet-Bert.tar.gz DeepAccNet-Bert
rm -rf DeepAccNet-Bert

echo Finished
