"""
Calculate MQA scores only for the resolved region from local score.

MQA methods:
    - DeepAccNet
    - P3CMQA
    - ProQ3D
    - VoroCNN
"""

import argparse
import os
import subprocess
import tarfile
from pathlib import Path
from typing import Any, List, Union

import numpy as np
import pandas as pd
from prody import parsePDB, writePDB
from tqdm import tqdm

data_dir = Path('../../data')
interim_path = data_dir / 'interim'
score_dir = data_dir / 'out/dataset/score/mqa'


def open_tar(tar_file: Union[str, Path]) -> tarfile.TarFile:
    return tarfile.open(tar_file, 'r:gz')


def get_resolved_pdb(target: str, resolved_indices: List[int]) -> Path:
    target_pdb_dir = data_dir / 'out/dataset/alphafold_output' / target
    pdb_resolved_dir = data_dir / 'out/dataset/pdb/pdb_resolved'
    pdb_resolved_target_dir = pdb_resolved_dir / target
    pdb_resolved_target_dir.mkdir(parents=True, exist_ok=True)
    for pdb in target_pdb_dir.glob('*.pdb'):
        pdb_name = pdb.stem
        output_pdb_path = pdb_resolved_target_dir / f'{pdb_name}.pdb'
        if output_pdb_path.exists():
            continue
        mol = parsePDB(pdb)
        resindices = mol.getResnums() - 1
        resolved_atom_indices = np.where(np.isin(resindices, resolved_indices))[0]
        mol_resolved = mol[resolved_atom_indices]
        writePDB(str(output_pdb_path), mol_resolved)
    return pdb_resolved_target_dir

class CalcResolvedConfidence:
    missing_dict = np.load(interim_path / 'missing_residues.npy', allow_pickle=True).item()

    def __init__(self, method: str, target_csv: Union[str, Path]):
        self.method = method
        self.target_df = pd.read_csv(target_csv, index_col=0)

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        results = []
        with tqdm(self.target_df.iterrows(), total=len(self.target_df)) as pbar:
            for _, row in pbar:
                target = row['id']
                pbar.set_description(f'Target = {target}')
                length = row['length']
                result = self.for_target(target, length)
                results.append(result)
        if sum([1 if result is None else 0 for result in results]) > 0:
            print(f'{self.method} calculation not yet finished')
            exit()
        return pd.concat(results)

    def for_target(self, target: str, length: int) -> Union[pd.DataFrame, None]:
        resolved_indices = self.get_resolved_indices(target, length)
        if self.method == 'DeepAccNet' or self.method == 'DeepAccNet-Bert':
            result = self.DeepAccNet(target, length)
        elif self.method == 'P3CMQA' or self.method == 'Sato-3DCNN':
            result = self.P3CMQA(target, resolved_indices)
        elif self.method == 'ProQ3D':
            result = self.ProQ3D(target, resolved_indices)
        elif self.method == 'VoroCNN':
            result = self.VoroCNN(target, resolved_indices)
        elif self.method == 'DOPE':
            result = self.DOPE(target, resolved_indices)
        elif self.method == 'SBROD':
            result = self.SBROD(target, resolved_indices)
        else:
            raise ValueError(f'Unknown method: {self.method}')
        return result

    @classmethod
    def get_resolved_indices(cls, target: str, length: int) -> List[int]:
        return np.setdiff1d(np.arange(length), cls.missing_dict[target])

    def DeepAccNet(self, target: str, length: int) -> Union[pd.DataFrame, None]:
        deepaccnet_path = score_dir / 'DeepAccNet'
        result_path = deepaccnet_path / f'{target}_resolved.csv'
        # if calculation already finished
        if result_path.exists():
            result_df = pd.read_csv(result_path, index_col=0)
            return result_df
        # if calculation not yet finished
        os.chdir('DeepAccNet')
        cmd = ['qsub', '-g', 'tga-ishidalab', './get_score_resolved.sh', target, str(length)]
        subprocess.run(cmd)
        os.chdir('..')
        return None

    def P3CMQA(self, target: str, resolved_indices: List[int]) -> pd.DataFrame:
        p3cmqa_path = score_dir / self.method
        tar_path = p3cmqa_path / f'{target}.tar.gz'
        tar = open_tar(tar_path)
        results = []
        for tarinfo in tar:
            if tarinfo.name.endswith('.csv'):
                if Path(tarinfo.name).stem == target:
                    continue
                f = tar.extractfile(tarinfo.name)
                local_df = pd.read_csv(f, index_col=0)
                resolved_score = np.mean(local_df['Score'][resolved_indices])
                results.append([Path(tarinfo.name).stem, resolved_score])
        result_df = pd.DataFrame(results, columns=['Model', f'{self.method}_resolved'])
        result_df['Target'] = target
        return result_df

    def ProQ3D(self, target: str, resolved_indices: List[int]) -> pd.DataFrame:
        proq3d_path = score_dir / self.method
        tar_path = proq3d_path / f'{target}.tar.gz'
        tar = open_tar(tar_path)
        results = []
        for tarinfo in tar:
            if tarinfo.name.endswith('.local'):
                f = tar.extractfile(tarinfo.name)
                local_df = pd.read_csv(f, sep=' ')
                resolved_score_dict = local_df.iloc[resolved_indices].mean().to_dict()
                resolved_score_dict['Model'] = Path(tarinfo.name).stem.split('.')[0]
                results.append(resolved_score_dict)
        result_df = pd.DataFrame(results)
        result_df['Target'] = target
        return result_df

    def VoroCNN(self, target: str, resolved_indices: List[int]) -> pd.DataFrame:
        vorocnn_path = score_dir / self.method
        tar_path = vorocnn_path / target / 'vorocnn_output.tar.gz'
        tar = open_tar(tar_path)
        results = []
        for tarinfo in tar:
            if tarinfo.name.endswith('.scores'):
                f = tar.extractfile(tarinfo.name)
                local_df = pd.read_csv(f, sep='\t')
                resolved_score = np.mean(local_df['prediction'][resolved_indices])
                results.append([Path(tarinfo.name).stem, resolved_score])
        result_df = pd.DataFrame(results, columns=['Model', f'{self.method}_resolved'])
        result_df['Target'] = target
        return result_df

    def DOPE(self, target: str, resolved_indices: List[int]) -> Union[pd.DataFrame, None]:
        dope_path = score_dir / 'DOPE'
        result_path = (dope_path / f'{target}_resolved.csv').resolve()
        # if calculation already finished
        if result_path.exists():
            result_df = pd.read_csv(result_path, index_col=0)
            result_df['Target'] = target
            return result_df
        # if calculation not yet finished
        resolved_pdb_dir = get_resolved_pdb(target, resolved_indices).resolve()
        os.chdir('DOPE')
        cmd = ['qsub', '-g', 'tga-ishidalab', './DOPE.sh', str(resolved_pdb_dir), str(result_path)]
        subprocess.run(cmd)
        os.chdir('..')
        return None

    def SBROD(self, target: str, resolved_indices: List[int]) -> Union[pd.DataFrame, None]:
        sbrod_path = score_dir / 'SBROD'
        result_path = (sbrod_path / f'{target}_resolved.csv').resolve()
        # if calculation already finished
        if result_path.exists():
            result_df = pd.read_csv(result_path, index_col=0)
            result_df['Target'] = target
            return result_df
        # if calculation not yet finished
        resolved_pdb_dir = get_resolved_pdb(target, resolved_indices).resolve()
        os.chdir('SBROD')
        cmd = ['qsub', '-g', 'tga-ishidalab', './SBROD.sh', str(resolved_pdb_dir), str(result_path)]
        subprocess.run(cmd)
        os.chdir('..')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('target_csv', type=str, help='Path to target csv file')
    parser.add_argument('method', type=str, help='Method name')
    args = parser.parse_args()
    method = args.method
    target_csv_path = Path(args.target_csv)
    dataset_name = target_csv_path.stem
    crc = CalcResolvedConfidence(method, target_csv_path)
    result_df = crc()
    print(result_df)
    output_score_path = data_dir / 'out/dataset/score/subsets' / dataset_name / f'{method}_resolved.csv.gz'
    result_df.to_csv(output_score_path)


if __name__ == '__main__':
    main()
