#!/bin/bash
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=3:00:00
#$ -j y

# arguments
input_dir=${1}
output_csv_path=${2}
output_dir=${output_csv_path%.*}
fasta=${3}
target=${4}
echo "Target: ${target}"

# load execution environment
source /home/4/16B09097/.bashrc
. /etc/profile.d/modules.sh
module load intel/17.0.4.196
module load cuda/9.2.148 cudnn/7.6
conda activate mypython

# path to the P3CMQA code
p3cmqa_dir=/gs/hs0/tga-ishidalab/takei/P3CMQA


# Set the path to the UniRef90 database
uniref90=/gs/hs0/tga-ishidalab/blastdb/uniref_2020_05/uniref90/uniref90
# preprocess
echo "Preprocessing..."
profile_dir=${output_dir}/profile
pushd ${p3cmqa_dir}/src
python preprocess.py -f ${fasta} -o ${profile_dir} -d ${uniref90} -n 14

# run P3CMQA
echo "Running P3CMQA..."
python predict.py -d ${input_dir} -f ${fasta} -o ${output_dir} -g 0 -p ${profile_dir}/${target} -s
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
