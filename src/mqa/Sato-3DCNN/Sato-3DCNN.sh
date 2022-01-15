#!/bin/bash
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=2:00:00
#$ -j y

# arguments
input_dir=${1}
output_csv_path=${2}
output_dir=${output_csv_path%.*}
target=${input_dir##*/}
echo "Target: ${target}"

# load execution environment
source /home/4/16B09097/.bashrc
. /etc/profile.d/modules.sh
module load intel/17.0.4.196
module load cuda/9.2.148 cudnn/7.6
conda activate mypython

# path to the P3CMQA code
p3cmqa_dir=/gs/hs0/tga-ishidalab/takei/P3CMQA

pushd ${p3cmqa_dir}/src

# run P3CMQA
echo "Running P3CMQA..."
python predict_atom_only.py -d ${input_dir} -o ${output_dir} -g 0 -s
popd

# post-process
if [ -e ${output_dir}/${target}.csv ]; then
    echo "Post-processing..."
    echo "python post_process.py ${output_dir} ${output_dir}/${target}.csv"
    python post_process.py ${output_dir} ${output_dir}/${target}.csv
    rm -f ${output_dir}/*.txt
    pushd ${output_dir}/..
    tar -cvzf ${target}.tar.gz ${target}
    rm -rf ${target}
    echo "Finished"
fi
