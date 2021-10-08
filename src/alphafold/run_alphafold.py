"""Read the csv of the list of targets and submit a job to run Alphafold on each target.
"""

import subprocess
import pandas as pd
from pathlib import Path
import datetime

data_dir = Path('../../data/')


def throw_job(fasta_path: Path, output_dir: Path, max_template_date: str, time: str, qsub=False):
    cmd = ['./alphafold.sh', str(fasta_path), str(output_dir), str(max_template_date)]
    if not qsub:
        print(' '.join(cmd))
    if qsub:
        cmd = ['qsub', '-g', 'tga-ishidalab',
               '-l', 'h_rt={}'.format(time), '-N', 'af_' + fasta_path.stem] + cmd
        subprocess.run(cmd)


def make_fasta(entry_id: str, header: str, seq: str, fasta_dir: Path) -> str:
    fasta_path = fasta_dir / (entry_id + '.fasta')
    with open(fasta_path, 'w') as f:
        f.write(header + '\n')
        len_line = 100
        for i in range((len(seq) - 1) // len_line + 1):
            f.write(seq[i * len_line: (i + 1) * len_line] + '\n')
        return fasta_path
    return None


def get_max_template_date_from_releasedate(releasedate: str):
    releasedate_str = releasedate.split('T')[0]
    releasedate_date = datetime.datetime.strptime(releasedate_str, '%Y-%m-%d')
    max_template_date: datetime.date = releasedate_date.date() - datetime.timedelta(days=1)
    return str(max_template_date)


def run_alphafold(row: pd.Series, fasta_dir: Path, pdb_dir: Path):
    entry_id = row['id']
    header, seq = row['header'], row['sequence']
    release_date = row['releasedate']
    max_template_date = get_max_template_date_from_releasedate(release_date)
    # length = row['length']
    # time = estimate_time_from_length(length)
    fasta_path = make_fasta(entry_id, header, seq, fasta_dir)
    assert fasta_path is not None
    time = '12:00:00'
    throw_job(fasta_path, pdb_dir, max_template_date, time)


def make_output_dir(output_dir_name: str) -> (Path, Path):
    # Structure of output dir.
    # out
    # └── output_dir_name
    #     ├── fasta
    #     └── pdb
    output_dir = data_dir / 'out' / output_dir_name
    output_dir.mkdir(parents=True, exist_ok=True)
    output_fasta_dir = output_dir / 'fasta'
    output_fasta_dir.mkdir(exist_ok=True)
    output_pdb_dir = output_dir / 'pdb'
    output_pdb_dir.mkdir(exist_ok=True)
    return output_fasta_dir, output_pdb_dir


def main():
    output_dir_name = 'test'
    fasta_dir, pdb_dir = make_output_dir(output_dir_name)
    target_list_path = data_dir / 'interim' / 'target_list.csv'
    df = pd.read_csv(target_list_path, index_col=0)
    df = df[: 1]
    for index, row in df.iterrows():
        run_alphafold(row, fasta_dir, pdb_dir)


if __name__ == '__main__':
    main()
