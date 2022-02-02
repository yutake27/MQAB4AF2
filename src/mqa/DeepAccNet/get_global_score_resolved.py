import argparse
import tarfile
from pathlib import Path
from typing import Any, List, Union

import numpy as np
import pandas as pd
from tqdm import tqdm

data_dir = Path('../../../data')
interim_path = data_dir / 'interim'
score_dir = data_dir / 'out/dataset/score/mqa'


def open_tar(tar_file: Union[str, Path]) -> tarfile.TarFile:
    return tarfile.open(tar_file, 'r:gz')


class CalcResolvedConfidence:
    missing_dict = np.load(interim_path / 'missing_residues.npy', allow_pickle=True).item()

    def __init__(self, method: str):
        self.method = method

    def for_target(self, target: str, length: int) -> pd.DataFrame:
        resolved_indices = self.get_resolved_indices(target, length)
        if self.method == 'DeepAccNet' or self.method == 'DeepAccNet-Bert':
            result = self.DeepAccNet(target, resolved_indices)
        else:
            raise ValueError(f'Unknown method: {self.method}')
        return result

    @classmethod
    def get_resolved_indices(cls, target: str, length: int) -> List[int]:
        return np.setdiff1d(np.arange(length), cls.missing_dict[target])

    def DeepAccNet(self, target: str, resolved_indices: List[int]) -> pd.DataFrame:
        deepaccnet_path = score_dir / 'DeepAccNet'
        tar_path = deepaccnet_path / target / f'{self.method}.tar.gz'
        tar = open_tar(tar_path)
        results = []
        for tarinfo in tqdm(tar):
            if tarinfo.name.endswith('.npz'):
                f = tar.extractfile(tarinfo.name)
                npz = np.load(f)
                resolved_score = np.mean(npz['lddt'][resolved_indices])
                results.append([Path(tarinfo.name).stem, resolved_score])
        result_df = pd.DataFrame(results, columns=['Model', f'{self.method}_resolved'])
        result_df['Target'] = target
        return result_df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('target', type=str, help='Target name')
    parser.add_argument('length', type=int, help='Target length')
    args = parser.parse_args()
    target, length = args.target, args.length
    standard_result = CalcResolvedConfidence('DeepAccNet').for_target(target, length)
    bert_result = CalcResolvedConfidence('DeepAccNet-Bert').for_target(target, length)
    result = pd.merge(standard_result, bert_result, on=['Target', 'Model'])
    output_path = score_dir / 'DeepAccNet' / f'{target}_resolved.csv'
    result.to_csv(output_path)


if __name__ == '__main__':
    main()
