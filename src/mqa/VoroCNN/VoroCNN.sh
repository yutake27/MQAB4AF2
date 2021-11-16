#!/bin/bash
#$ -cwd
#$ -l s_core=1
#$ -l h_rt=2:00:00
#$ -j y

input_dir=${1}
output_csv_path=${2}
output_dir=${output_csv_path%.*}

# load python environment
. /etc/profile.d/modules.sh
module load python/3.9.2
source ../../../.venv/bin/activate
# Set vorocnn environment
vorocnn=/gs/hs0/tga-ishidalab/takei/bin/vorocnn/vorocnn # Set your vorocnn path
voronota_path=${vorocnn%/*}/voronota/linux/voronota

# Running VoroCNN
${vorocnn} -i ${input_dir} -o ${output_dir} -v ${voronota_path} -V

# Convert to local score to global score
python get_global_score.py ${output_dir}/vorocnn_output -o ${output_csv_path}

# Compress local score
cd ${output_dir}
tar -cvzf vorocnn_output.tar.gz vorocnn_output
rm -rf vorocnn_output

echo "Finished"
