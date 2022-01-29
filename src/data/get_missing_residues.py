import json
from operator import index
from pathlib import Path

import numpy as np
import pandas as pd
from prody import parsePDB
from tqdm import tqdm

native_dir = Path('../../data/out/dataset/native_pdb')


def get_missing_residues(target: str, length: int) -> list:
    native_path = native_dir / f'{target}.pdb'
    mol = parsePDB(str(native_path)).select('name CA')
    resolved_indices = mol.getResindices()
    all_indices = np.arange(length)
    missing_indices = list(np.setdiff1d(all_indices, resolved_indices))
    return missing_indices


def main():
    df = pd.read_csv('../../data/interim/target_list.csv', index_col=0)
    results_dict = {}
    for i, rows in tqdm(df.iterrows()):
        target, length = rows['id'], rows['length']
        missing_indices = get_missing_residues(target, length)
        results_dict[target] = missing_indices
    np.save('../../data/interim/missing_residues.npy', results_dict)

if __name__ == '__main__':
    main()
