"""
Read the csv of the target list and submit jobs to run Alphafold on each target.
"""

import argparse
import datetime
import math
import subprocess
import time
from pathlib import Path
from typing import List, Tuple

import pandas as pd

# Directory of the data
data_dir = Path('../../data/').resolve()


class RunMSAGenerator:
    @classmethod
    def run(cls, row: pd.Series, fasta_dir: Path, output_dir: Path):
        """throw a job to run alphafold or colabfold on a target

        Args:
            row (pd.Series): a row of dataframe of the target list
            fasta_dir (Path): directory of the fasta files to save
            output_dir (Path): directory of the alphafold output files
            method (str): 'alphafold' or 'colabfold'
            qsub (bool): whether to submit the job
        """
        entry_id = row['id']
        header, seq = row['header'], row['sequence']
        output_target_dir = output_dir / entry_id
        output_msa_pickle = output_target_dir / 'msa.pickle'
        if output_msa_pickle.exists():
            return
        fasta_path = make_fasta(entry_id, header, seq, fasta_dir)
        assert fasta_path.exists()
        cmd = ['./msa_only.sh', '-f', str(fasta_path), '-o', str(output_target_dir)]
        print(' '.join(cmd))
        subprocess.run(cmd)
        time.sleep(5 * 60)  # Wait some minutes after job submission to avoid overloading the MMseqs2 server


def make_fasta(entry_id: str, header: str, seq: str, fasta_dir: Path) -> Path:
    fasta_path = fasta_dir / (entry_id + '.fasta')
    if fasta_path.exists():
        return fasta_path
    with open(fasta_path, 'w') as f:
        f.write(header + '\n')
        len_line = 100
        for i in range((len(seq) - 1) // len_line + 1):
            f.write(seq[i * len_line: (i + 1) * len_line] + '\n')
        return fasta_path


def get_max_template_date_from_releasedate(releasedate: str):
    releasedate_str = releasedate.split('T')[0]
    releasedate_date = datetime.datetime.strptime(releasedate_str, '%Y-%m-%d')
    max_template_date: datetime.date = releasedate_date.date() - datetime.timedelta(days=1)
    return str(max_template_date)


def make_output_dir(output_dir_name: str) -> Tuple[Path, Path]:
    # Structure of output directory:
    # data/out
    # └── output_dir_name
    #     ├── fasta
    #     └── pdb
    #     └── native_pdb
    #     └── score
    #       └── label
    #     └── alphafold_output
    output_dir = data_dir / 'out' / output_dir_name
    output_dir.mkdir(parents=True, exist_ok=True)
    output_fasta_dir = output_dir / 'fasta'
    output_fasta_dir.mkdir(exist_ok=True)
    output_pdb_dir = output_dir / 'pdb'
    output_pdb_dir.mkdir(exist_ok=True)
    output_native_pdb_dir = output_dir / 'native_pdb'
    output_native_pdb_dir.mkdir(exist_ok=True)
    output_score_dir = output_dir / 'score'
    output_score_dir.mkdir(exist_ok=True)
    output_label_dir = output_score_dir / 'label'
    output_label_dir.mkdir(exist_ok=True)
    output_alphafold_dir = output_dir / 'alphafold_output'
    output_alphafold_dir.mkdir(exist_ok=True)
    return output_fasta_dir, output_alphafold_dir


def main():
    parser = argparse.ArgumentParser(description='Run Alphafold on the target list.')
    parser.add_argument('dataset_name', type=str, help='Name of the dataset.')
    parser.add_argument('-t', '--target_list_path', type=str, help='Path to the target list.')
    parser.add_argument('-n', '--num_targets', type=int, default=-1,
                        help='Number of targets to run. -1 means all.')
    args = parser.parse_args()

    num_targets = args.num_targets

    output_dir_name = args.dataset_name
    fasta_dir, output_alphafold_dir = make_output_dir(output_dir_name)
    if args.target_list_path is None:
        target_list_path = data_dir / 'interim' / 'target_list.csv'
    else:
        target_list_path = args.target_list_path
    df = pd.read_csv(target_list_path, index_col=0)
    num_targets = num_targets if num_targets > 0 else len(df)
    df = df[: num_targets]
    for index, row in df.iterrows():
        RunMSAGenerator.run(row, fasta_dir, output_alphafold_dir)


if __name__ == '__main__':
    main()
