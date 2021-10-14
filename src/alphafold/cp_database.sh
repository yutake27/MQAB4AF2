#!/bin/bash

. /etc/profile.d/modules.sh
module purge
module load cuda/11.0.3 alphafold/2.0.0
module li


echo Copy database to '$TMPDIR' ($TMPDIR)
echo =================================
echo `date` cp uniclust30 database
cp $ALPHAFOLD_DATA_DIR/uniclust30/uniclust30_2018_08_hhsuite.tar.gz $TMPDIR
echo `date` finish copy
pushd $TMPDIR
echo `date` compress uniclust30_2018_08_hhsuite.tar.gz
tar xvf uniclust30_2018_08_hhsuite.tar.gz --use-compress-prog=pigz
popd
# echo =================================
# echo `date` cp mgnify
# cp $ALPHAFOLD_DATA_DIR/mgnify/mgy_clusters_2018_12.fa $TMPDIR
# echo =================================
# echo `date` cp uniref90
# cp $ALPHAFOLD_DATA_DIR/uniref90/uniref90.fasta $TMPDIR
# echo =================================
# echo `date` cp pdb70
# mkdir $TMPDIR/pdb70
# cp $ALPHAFOLD_DATA_DIR/pdb70/pdb70_from_mmcif_200401.tar.gz $TMPDIR/pdb70
# pushd $TMPDIR/pdb70
# tar xvf pdb70_from_mmcif_200401.tar.gz --use-compress-prog=pigz
# popd
echo =================================
echo `date` finish
