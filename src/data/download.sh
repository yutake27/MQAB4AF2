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

# download newest cath domain classification
mkdir -p cath
pushd cath
cath_newest=cath-b-newest-all.gz
if [ ! -e ${cath_newest} ]; then
    wget http://download.cathdb.info/cath/releases/daily-release/newest/cath-b-newest-all.gz
    gzip -d cath-b-newest-all.gz
fi
popd

# download ecod domain definition (version 282, 2021-10-04)
mkdir ecod
pushd ecod
ecod_domains=ecod.develop282.domains.txt
if [ ! -e ${ecod_domains} ]; then
    wget http://prodata.swmed.edu/ecod/distributions/${ecod_domains}
fi
popd

# Finish
popd