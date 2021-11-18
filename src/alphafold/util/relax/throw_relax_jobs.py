"""
Read the csv of the target list and submit jobs to relax predicted structures.
"""

import argparse
import datetime
import math
import subprocess
from pathlib import Path
from typing import List, Tuple

import pandas as pd

# Directory of the data
data_dir = Path('../../../../data/').resolve()
dataset_name = 'dataset'
dataset_dir = data_dir / 'out' / dataset_name
alphafold_output_dir = dataset_dir / 'alphafold_output'
native_pdb_dir = dataset_dir / 'native_pdb'
score_dir = dataset_dir / 'score'
relaxed_label_dir = score_dir / 'label_relaxed'
relaxed_label_dir.mkdir(exist_ok=True)


class ThrowJob:
    @staticmethod
    def _throw_job_from_cmd(cmd: List, qsub=False):
        if not qsub:
            print(' '.join(cmd))
        if qsub:
            subprocess.run(cmd)

    @classmethod
    def throw_job(cls, row: pd.Series, qsub: bool):
        """throw a job to relax structures

        Args:
            row (pd.Series): a row of dataframe of the target list
            qsub (bool): whether to submit the job
        """
        target = row['id']
        output_target_dir = alphafold_output_dir / target
        # Check if alphafold execution for the target is already done
        if not (output_target_dir / 'model_pickle.tar.gz').exists():
            print(f'AlphaFold prediction for {target} has not been completed!')
            return
        native_pdb_path = native_pdb_dir / (target + '.pdb')
        output_score_path = relaxed_label_dir / (target + '.csv')
        # Check if relax for the target is already done
        if output_score_path.exists():
            print(f'Already finished: {target}')
            return
        qsub_header = ['qsub', '-g', 'tga-ishidalab', '-N', f'relax_{target}']
        main_cmd = ['./relax.sh', str(output_target_dir), str(native_pdb_path), str(output_score_path)]
        cmd = qsub_header + main_cmd
        cls._throw_job_from_cmd(cmd, qsub)


def main():
    parser = argparse.ArgumentParser(description='Throw relax jobs for Alphafold structures on the target list.')
    parser.add_argument('-t', '--target_list_path', type=str, help='Path to the target list.')
    parser.add_argument('-n', '--num_targets', type=int, default=-1,
                        help='Number of targets to run. -1 means all.')
    parser.add_argument('-q', '--qsub', action='store_true', help='Submit the job to the queue.')
    args = parser.parse_args()

    num_targets = args.num_targets
    qsub = args.qsub

    if args.target_list_path is None:
        target_list_path = data_dir / 'interim' / 'target_list.csv'
    else:
        target_list_path = args.target_list_path
    df = pd.read_csv(target_list_path, index_col=0)
    num_targets = num_targets if num_targets > 0 else len(df)
    df = df[: num_targets]
    for index, row in df.iterrows():
        ThrowJob.throw_job(row, qsub)


if __name__ == '__main__':
    main()
