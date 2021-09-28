#!/bin/bash

data_dir=../../data/raw
pushd ${data_dir}

# download bc40.out
bc40=bc-40.out
if [ ! -e ${bc40} ]; then
    wget https://cdn.rcsb.org/resources/sequence/clusters/bc-40.out
fi

# download pdb seqres
pdb_seqres=pdb_seqres.txt
if [ ! -e ${pdb_seqres} ]; then
    wget https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
    gzip -d pdb_seqres.txt.gz
fi

# Finish
popd