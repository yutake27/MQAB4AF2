"""
Calculate confidence scores only for the resolved region.
"""

import argparse
import pickle
import sys
import tarfile
from ctypes import Union
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.append('../alphafold/alphafold/common/')
import confidence

sys.path.append('../alphafold')

data_dir = Path('../../../data')
score_dir = data_dir / 'out/dataset/alphafold_output'


def calc_plddt_resolved_region(prediction_result: Dict, resolved_resindices: np.ndarray) -> float:
    plddt = prediction_result['plddt']
    plddt_resolved = plddt[resolved_resindices]  # get plddt for resolved residues
    plddt_resolved_mean = np.mean(plddt_resolved) / 100
    assert 0 <= plddt_resolved_mean <= 1
    return float(plddt_resolved_mean)


def calc_ptm_from_pae(pae: np.ndarray, residue_weights: np.ndarray) -> float:
    num_res = np.sum(residue_weights)
    # Clip num_res to avoid negative/undefined d0.
    clipped_num_res = max(num_res, 19)

    # Compute d_0(num_res) as defined by TM-score, eqn. (5) in
    # http://zhanglab.ccmb.med.umich.edu/papers/2004_3.pdf
    # Yang & Skolnick "Scoring function for automated
    # assessment of protein structure template quality" 2004
    d0 = 1.24 * (clipped_num_res - 15) ** (1./3) - 1.8
    tm_per = 1. / (1 + np.square(pae) / np.square(d0))
    tm_per_masked = tm_per * residue_weights  # mask unresolved residues
    predicted_tm_term = np.sum(tm_per_masked, axis=-1)
    normed_residue_mask = residue_weights / (1e-8 + residue_weights.sum())
    ptm_per_res = predicted_tm_term * normed_residue_mask
    ptm = np.max(ptm_per_res, axis=-1)  # find best ptm
    assert 0 <= ptm <= 1
    return float(ptm)


def calc_ptm_resolved_region(prediction_result: Dict, residue_weights: np.ndarray) -> float:
    pae = prediction_result['pae']
    ptm_resolved = calc_ptm_from_pae(pae, residue_weights=residue_weights)
    return ptm_resolved


def calc_confidence_resolved_region(prediction_result: Dict, missing_indices: List, length: int
                                    ) -> Tuple[float, float]:
    residue_weights = np.zeros(length)
    resolved_resindices = np.setdiff1d(np.arange(length), missing_indices)
    residue_weights[resolved_resindices] = 1
    plddt_resolved_mean = calc_plddt_resolved_region(prediction_result, resolved_resindices)
    if 'pae' in prediction_result:
        ptm_resolved = calc_ptm_resolved_region(prediction_result, residue_weights)
    else:
        ptm_resolved = None
    return plddt_resolved_mean, ptm_resolved


def calc_confidence_tar(tar_file, missing_indices, length):
    tar = tarfile.open(tar_file, 'r:gz')
    confidences = []
    for tarinfo in tar:
        if tarinfo.name.endswith('.pickle'):
            f = tar.extractfile(tarinfo.name)
            prediction_result = pickle.load(f)
            confidence = calc_confidence_resolved_region(prediction_result, missing_indices, length)
            confidences.append([tarinfo.name.split('.')[0], *confidence])
    return confidences


def main_target(target_name: str, missing_indices: List, length: int) -> pd.DataFrame:
    tar_file = score_dir / target_name / 'model_pickle.tar.gz'
    confidences = calc_confidence_tar(tar_file, missing_indices, length)
    df = pd.DataFrame(confidences, columns=['Model', 'pLDDT_resolved', 'pTM_resolved'])
    df['Target'] = target_name
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('target_csv', type=str, help='Path to target csv file')
    args = parser.parse_args()
    target_csv_path = Path(args.target_csv)
    dataset_name = target_csv_path.stem
    target_df = pd.read_csv(target_csv_path, index_col=0)
    interim_path = data_dir / 'interim'
    missing_dict = np.load(interim_path / 'missing_residues.npy', allow_pickle=True).item()
    results = []
    with tqdm(target_df.iterrows(), total=len(target_df)) as pbar:
        for i, row in pbar:
            target = row['id']
            pbar.set_description(f'Target = {target}')
            length = row['length']
            missing_indices = missing_dict[target]
            target_result_df = main_target(target, missing_indices, length)
            results.append(target_result_df)
    result_df = pd.concat(results)
    output_score_path = data_dir / 'out' / 'score' / 'subsets' / dataset_name / 'af2_confidence_resolved.csv.gz'
    result_df.to_csv(output_score_path)


if __name__ == '__main__':
    main()
