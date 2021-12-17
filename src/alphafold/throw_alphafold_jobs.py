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


class DetermineJobParams:
    @staticmethod
    def estimate_time_from_length(length: int, ensemble: str = 'full',
                                  resume=False, resume_divider: int = 4) -> str:
        """
        Estimate the time required to run a job based on the length of the target.

        Args:
        length (int): Length of the target.
        ensemble (str): Whether to run ensemble or not. Choices=['full', 'ens', 'noens']. Default is full.
        resume (bool): Whether to resume or not. Default is False.
        resume_divider (int): Number to divide the estimated time when resuming. Default is 4.

        Returns:
            str: Estimated time required to run the job.
        """
        assert ensemble in ['full', 'ens', 'noens']
        fitting_params = {
            # Slope and intercept obtained from fitting curve in sample running.
            # The 24-hour rate is calculated. For example, 0.5 means 12 hours.
            'full': {
                's': 0.001410034,
                'i': -0.01651
            },
            'ens': {
                's': 0.001602107,
                'i': -0.082827655
            },
            'noens': {  # calculate only from sequence with over 400 residues
                's': 0.000424883,
                'i': -0.094979556
            },
        }
        estimated_hour = (fitting_params[ensemble]['s'] * length + fitting_params[ensemble]['i']) * 24
        spare_hour = 1.5
        total_time = estimated_hour + spare_hour
        if resume:
            total_time /= resume_divider
        _minute, hour = math.modf(total_time)
        assert _minute < 1
        minute = _minute * 60
        return f'{int(hour)}:{int(minute):02d}:00'

    @staticmethod
    def determine_run_ensemble_from_length(length: int) -> bool:
        """
        Determine whether to run the ensemble or not based on the length of the target.
        """
        if length < 400:
            return True
        else:
            return False


class checkSubmittedJob:
    def __init__(self):
        self.submitted_job_targets = self._get_submitted_job_targets()

    @staticmethod
    def _qstat() -> List[str]:
        """Return job id list"""
        cmd = ['qstat']
        qstat_result = subprocess.run(cmd, capture_output=True, text=True).stdout
        qstat_result_lines = qstat_result.split('\n')[2: -1]
        qstat_result_lines_filter = list(filter(lambda line: 'af_' in line, qstat_result_lines))
        job_ids = [line.split()[0] for line in qstat_result_lines_filter]
        return job_ids

    @staticmethod
    def _stat_j(job_id: str) -> str:
        """Return target name"""
        cmd = ['qstat', '-j', job_id]
        qstat_j_result = subprocess.run(cmd, capture_output=True, text=True).stdout
        qstat_j_result_line = qstat_j_result.split('\n')
        job_name_line = list(filter(lambda x: x.startswith('job_name'), qstat_j_result_line))[0]
        ar_id_line = list(filter(lambda x: x.startswith('ar_id'), qstat_j_result_line))
        print(job_name_line, ar_id_line[0] if ar_id_line else '')
        job_name = job_name_line.split()[1]
        target_name = job_name.split('_', 1)[1].rsplit('_', 1)[0]
        return target_name

    @classmethod
    def _get_submitted_job_targets(cls) -> List[str]:
        """Return target name list"""
        job_ids = cls._qstat()
        target_names = [cls._stat_j(job_id) for job_id in job_ids]
        return target_names

    def check_submitted_jobs(self, target_name: str) -> bool:
        """
        Check whether the job of the target is submitted or not.

        Args:
        target_name (str): Name of the target.

        Returns:
            bool: Whether the job is submitted or not.
        """
        return target_name in self.submitted_job_targets


class ThrowJob:
    @staticmethod
    def _throw_job_from_cmd(cmd: List, qsub=False):
        if not qsub:
            print(' '.join(cmd))
        if qsub:
            subprocess.run(cmd)

    @classmethod
    def throw_job(cls, row: pd.Series, fasta_dir: Path, output_dir: Path,
                  method: str, qsub: bool, runall: bool = False):
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
        release_date = row['releasedate']
        length = row['length']
        output_target_dir = output_dir / entry_id
        output_score_path = output_target_dir / 'scores.csv'
        if output_score_path.exists() and not runall:  # check if the output file already exists
            tar_pkl_path = output_target_dir / 'model_pickle.tar.gz'
            if not tar_pkl_path.exists():  # check if post processing is done
                cmd = ['python', 'post_process.py', f'{output_target_dir}']
                print(' '.join(cmd))
            return
        fasta_path = make_fasta(entry_id, header, seq, fasta_dir)
        assert fasta_path.exists()
        # Determine run ensemble or not
        run_ensemble = DetermineJobParams.determine_run_ensemble_from_length(length)
        if run_ensemble:
            ensemble = 'full'
        else:
            ensemble = 'noens'
        # Determine the time required to run the job
        resume: bool = len(list(output_target_dir.glob('*.pdb'))) > 100
        qsub_time = DetermineJobParams.estimate_time_from_length(length, ensemble, resume)
        qsub_header = ['qsub', '-g', 'tga-ishidalab', '-l', f'h_rt={qsub_time}', '-N', f'af_{entry_id}_{length}']
        if method == 'alphafold':
            max_template_date = get_max_template_date_from_releasedate(release_date)
            cmd = ['./alphafold.sh', str(fasta_path), str(output_target_dir), str(max_template_date)]
        elif method == 'colabfold':
            cmd = ['./colabfold.sh', '-f', str(fasta_path), '-o', str(output_target_dir)]
            if not run_ensemble:  # only use ensemble 1
                cmd += ['-e', '1']
        else:
            return
        cmd = qsub_header + cmd
        cls._throw_job_from_cmd(cmd, qsub)
        if qsub and not resume:
            time.sleep(10 * 60)  # Wait some minutes after job submission to avoid overloading the MMseqs2 server

    @classmethod
    def throw_job_for_target_df(cls, df: pd.DataFrame, **kwargs) -> None:
        crj = checkSubmittedJob()
        for index, row in df.iterrows():
            target_name = str(row['id'])
            if not crj.check_submitted_jobs(target_name):
                cls.throw_job(row, **kwargs)


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
    parser.add_argument('--method', default='colabfold', choices=['alphafold', 'colabfold'],
                        type=str, help='Alphafold or Colabfold.')
    parser.add_argument('-q', '--qsub', action='store_true', help='Submit the job to the queue.')
    parser.add_argument('--runall', action='store_true',
                        help='Run command even for targets that have already been completed')
    args = parser.parse_args()

    num_targets = args.num_targets
    method = args.method
    qsub = args.qsub
    runall = args.runall

    output_dir_name = args.dataset_name
    fasta_dir, output_alphafold_dir = make_output_dir(output_dir_name)
    if args.target_list_path is None:
        target_list_path = data_dir / 'interim' / 'target_list.csv'
    else:
        target_list_path = args.target_list_path
    df = pd.read_csv(target_list_path, index_col=0)
    num_targets = num_targets if num_targets > 0 else len(df)
    df = df[: num_targets]
    ThrowJob.throw_job_for_target_df(df, fasta_dir=fasta_dir, output_dir=output_alphafold_dir,
                                     method=method, qsub=qsub, runall=runall)


if __name__ == '__main__':
    main()
