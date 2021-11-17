"""
Specify the MQA method and throw jobs of that method
Also specify the target list.

When throwing a job, check the targets for which AlphaFold could not be executed,
and do not execute MQA on these targets.

If all the targets that should be predicted have been executed, combine the scores.
Save in subsets/subset_name/mqa_method.csv
"""
import argparse
import os
import subprocess
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from prody import parsePDB, confProDy

confProDy(verbosity='none')

dataset_dir = Path('../../data/out/dataset').resolve()


def not_be_modeled(target_pdb_dir: Path) -> bool:
    """
    Check if the AlphaFold execution was done successfully.
    (the coordinates of the output pdb are not nan)
    """
    sample_model = next(Path(target_pdb_dir).glob('*.pdb'))
    mol = parsePDB(sample_model)
    coords = mol.getCoords()
    return np.isnan(coords).any()


def throw_job(target_pdb_dir: Path, method: str, output_score_path: Path, qsub: bool, execute: bool):
    """
    Throw a job of the specified MQA method.
    """
    target = target_pdb_dir.stem
    cmd = [f'./{method}.sh', str(target_pdb_dir), str(output_score_path)]
    if method in ['ProQ3D', 'P3CMQA']:  # If fasta path needed
        fasta_path = dataset_dir / 'fasta' / f'{target}.fasta'
        cmd += [str(fasta_path), str(target)]
    if qsub:
        cmd = ['qsub', '-g', 'tga-ishidalab', '-N', f'{method}_{target}'] + cmd
    print(' '.join(cmd))
    if execute:
        subprocess.run(cmd)


def concatenate_scores(mqa_output_dir: Path, target_list: List, subset_score_path: Path):
    """
    Concatenate mqa scores for targets in the target list.
    """
    df_list = []
    for target in target_list:
        target_df = pd.read_csv(mqa_output_dir / f'{target}.csv', index_col=0)
        target_df['Target'] = target
        df_list.append(target_df)
    df = pd.concat(df_list)
    df.to_csv(subset_score_path)


def main():
    parser = argparse.ArgumentParser(description='Throw MQA job')
    parser.add_argument('-m', '--method', type=str, required=True, help='MQA method name',
                        choices=['DOPE', 'ProQ3D', 'SBROD', 'P3CMQA', 'DeepAccNet', 'VoroCNN'])
    parser.add_argument('-t', '--target_list_csv', type=str, required=True, help='target list csv')
    parser.add_argument('-q', '--qsub', action='store_true', help='throw job using qsub')
    parser.add_argument('-e', '--execute', action='store_true', help='execute the job')
    args = parser.parse_args()

    method = args.method
    target_list_csv = Path(args.target_list_csv)
    subset_name = target_list_csv.stem
    target_df = pd.read_csv(target_list_csv, index_col=0)
    target_list = target_df['id'].tolist()

    """Output score directory structure
    └── score
        ├── label
        ├── mqa
        │   └── method
        └── subsets
            └── subset
    """
    score_dir = dataset_dir / 'score'
    alphafold_output_dir = dataset_dir / 'alphafold_output'
    mqa_output_dir = score_dir / 'mqa' / method
    mqa_output_dir.mkdir(parents=True, exist_ok=True)
    subset_score_dir = score_dir / 'subsets' / subset_name
    subset_score_dir.mkdir(parents=True, exist_ok=True)

    # Throw MQA job
    uncompleted_jobs: int = 0
    os.chdir(method)
    print('pwd:', os.getcwd())
    final_target_list = target_list.copy()
    for target in target_list:
        # pdb directory
        target_pdb_dir = alphafold_output_dir / target
        # check if the AlphaFold execution was done successfully
        if not_be_modeled(target_pdb_dir):
            final_target_list.remove(target)
            continue
        # check if the MQA execution was already done
        output_score_path = mqa_output_dir / f'{target}.csv'
        if output_score_path.exists():
            continue
        # throw job
        uncompleted_jobs += 1
        throw_job(target_pdb_dir, method, output_score_path, args.qsub, args.execute)

    # If all jobs are completed, concatenate scores
    if uncompleted_jobs == 0:
        print('All jobs are completed. Concatenating scores...')
        subset_score_path = subset_score_dir / f'{method}.csv'
        concatenate_scores(mqa_output_dir, final_target_list, subset_score_path)


if __name__ == '__main__':
    main()
