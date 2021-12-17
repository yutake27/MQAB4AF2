"""
Calculate Neff (Number of effective sequences) in MSA using CD-HIT

The calculation method was adapted from P. Bryant, et al.,
"Improved prediction of protein-protein interactions using AlphaFold2", bioRxiv, 2021.
(https://www.biorxiv.org/content/10.1101/2021.09.15.460468v2)

1. Read two MSA files (bfd.a3m, uniref.a3m) and merge them into one MSA file
2. Remove gaps in MSA
3. Clustering MSA using CD-HIT
4. Set the number of clusters to Neff
"""

import argparse
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Union

import pandas as pd
import tqdm

data_dir = Path('../../../data/out/dataset/')
STRLIKE = Union[str, Path]


class CalculateNeff:
    @staticmethod
    def clstr_msa_by_cdhit(input_file: STRLIKE, output_file: STRLIKE,
                           verbose: bool = False, thread: int = 7) -> None:
        """
        Clustering MSA by cd-hit
        (cd-hit -i input_file -o output_file -c 0.62 -G 0 -n 3 -aS 0.9)
        """
        cmd = ['cd-hit', '-i', str(input_file), '-o', str(output_file),
               '-c', '0.62', '-G', '0', '-n', '3', '-aS', '0.9', '-T', str(thread)]
        output = subprocess.run(cmd, capture_output=True, text=True).stdout
        if verbose:
            print(output)

    @staticmethod
    def count_sequence_in_fasta(fasta_file: STRLIKE) -> int:
        """
        Count number of sequences in fasta file
        """
        with open(fasta_file, 'r') as f:
            lines = f.readlines()
            lines = list(filter(lambda x: x.startswith('>'), lines))
        return len(lines)

    @classmethod
    def calculate_neff(cls, msa_file: STRLIKE) -> Dict[str, int]:
        """
        Calculate Neff by cd-hit

        Return:
            Dict[str, int]: {'Neff': Neff, 'NumMSA': number of sequences in MSA}
        """
        lines = cls.read_MSA(msa_file)
        assert len(lines) > 0
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_file = Path(tmpdir) / 'tmp.fasta'
            tmp_file.write_text('\n'.join(lines))
            clstr_output = tmp_file.with_suffix('.clstr')
            cls.clstr_msa_by_cdhit(tmp_file, clstr_output)
            neff = cls.count_sequence_in_fasta(clstr_output)
            num_msa = cls.count_sequence_in_fasta(msa_file)
            assert neff <= num_msa
            return {'Neff': neff, 'NumMSA': num_msa}

    @staticmethod
    def remove_gap(seq: str) -> str:
        """
        Remove gap in sequence
        """
        return seq.replace('-', '')

    @classmethod
    def read_MSA(cls, msa_file: STRLIKE) -> list:
        """
        Read MSA file
        """
        with open(msa_file, 'r') as f:
            lines = f.readlines()
        lines = [cls.remove_gap(line.strip()) if line[0] != '>' else line for line in lines]
        return lines


def concat_MSA(msa_files) -> List[str]:
    """
    Concatenate MSA files and return lines
    """
    lines = []
    for msa_file in msa_files:
        with open(msa_file, 'r') as f:
            lines += [line.strip() for line in f.readlines() if line[0] != '\x00']
    return lines


def calc_neff_for_msa_files(msa_files) -> Dict[str, int]:
    """
    Calculate Neff for msa files (Concatenate MSA files and calculate Neff)
    """
    lines = concat_MSA(msa_files)
    with tempfile.NamedTemporaryFile('w') as f:
        f.write('\n'.join(lines))
        f.seek(0)
        neff = CalculateNeff.calculate_neff(f.name)
    num_seqs = sum([CalculateNeff.count_sequence_in_fasta(msa_file) for msa_file in msa_files])
    assert neff['NumMSA'] == num_seqs
    neff['NumMSA'] -= 1
    return neff


def calc_neff_for_target(target: str) -> Dict[str, int]:
    target_dir = data_dir / 'alphafold_output' / target
    assert target_dir.exists()
    msa_dir = next(target_dir.glob('*_env'))
    msa_files = list(msa_dir.glob('*.a3m'))
    return calc_neff_for_msa_files(msa_files)


def main():
    parser = argparse.ArgumentParser(description='Calculate Neff for each target in the subset')
    parser.add_argument('target_list', type=str, help='Target list file')
    args = parser.parse_args()

    target_list = args.target_list
    subset_name = Path(target_list).stem
    # Read target list
    target_df = pd.read_csv(args.target_list, index_col=0)
    neffs = []
    for target in tqdm.tqdm(target_df['id']):
        neff = calc_neff_for_target(target)
        neff['Target'] = target
        neffs.append(neff)
        print(target, neff)
    result = pd.DataFrame(neffs)
    print(result)
    output_path = data_dir / 'score' / 'subsets' / subset_name / 'neff.csv'
    result.to_csv(output_path)

def test_sample_target():
    target = '6AN4_A'
    # target = '6EXU_A'
    calc_neff_for_target(target)


def test_calc_neff():
    bfd_a3m = 'sample_msa/bfd.mgnify30.metaeuk30.smag30.a3m'
    uniref_a3m = 'sample_msa/uniref.a3m'
    msa_files = [bfd_a3m, uniref_a3m]
    neff = calc_neff_for_msa_files(msa_files)
    print(neff)


if __name__ == '__main__':
    main()
