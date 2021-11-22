#!/bin/bash
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=3:00:00
#$ -j y

. /etc/profile.d/modules.sh
module load cuda/10.1.105 cudnn/7.6 openmpi r

# load environment for ProQ3D
source /home/4/16B09097/.bashrc
conda activate keras
python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"


# arguments
input_dir=${1}
output_csv_path=${2}
output_dir=${output_csv_path%.*}
fasta=${3}
target=${4}

if [ -e ${output_csv_path} ]; then
    echo "proq3 for ${target} have already been completed!"
    exit
fi

# make output directory
mkdir ${output_dir}

# make a list of pdb files
txtpath=${output_dir}/${target}.txt
if [ -e ${txtpath} ]; then
    rm ${txtpath}
fi

for i in ${input_dir}/*.pdb;do
    echo $i >> ${txtpath}
done


# Execute ProQ3D
echo "Running ProQ3D..."
run_proq3.sh -fasta ${fasta} -l ${txtpath} -outpath ${output_dir} -ncores 14 -repack no
# post-process
echo "Post-processing..."
python post_process.py ${output_dir}

# Compress output
cd ${output_dir}/..
tar -cvzf ${target}.tar.gz ${target}
rm -rf ${target}

echo "Finished"
