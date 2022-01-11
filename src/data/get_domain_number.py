"""
Get domain number from CATH and ECOD database
(download data by ./download.sh)
"""

import os
import subprocess
from pathlib import Path

import pandas as pd

data_dir = Path('../../data').absolute()
row_dir = data_dir / 'raw'
interim_dir = data_dir / 'interim'


def get_cath_domain():
    cath_domain_path = row_dir / 'cath' / 'cath-b-newest-all'
    cath_domain_df = pd.read_csv(cath_domain_path, sep=' ', header=None)
    cath_domain_df.columns = ['id', '_', '_', 'res']
    cath_domain_df['PDB_ID'] = cath_domain_df['id'].str[:4].str.upper()
    cath_domain_df['Chain'] = cath_domain_df['id'].str[4]
    cath_domain_df['id'] = cath_domain_df['PDB_ID'] + '_' + cath_domain_df['Chain']
    domain_num_df = pd.DataFrame(cath_domain_df.groupby(['id']).size(), columns=['num_domain']).reset_index()
    domain_num_df.to_csv(interim_dir / 'cath_domain_num.csv', index=False)


def get_ecod_domain():
    ecod_domain_path = row_dir / 'ecod' / 'ecod.develop282.domains.txt'
    ecod_df = pd.read_csv(ecod_domain_path, skiprows=4, sep='\t')
    ecod_df['id'] = ecod_df['pdb'] + '_' + ecod_df['chain']
    ecod_df['id'] = ecod_df['id'].str.upper()
    domain_num_df = pd.DataFrame(ecod_df.groupby(['id']).size(), columns=['num_domain']).reset_index()
    domain_num_df.to_csv(interim_dir / 'ecod_domain_num.csv', index=False)


def main():
    get_cath_domain()
    get_ecod_domain()


if __name__ == '__main__':
    main()
