"""
Detect if there are low interaction region in the protein from contact map

Calculate the percentage of contacts between the residues before the i+1-th residue and
the residues after the i-th residue, and determine that there is a low interaction region
if the percentage is less than n%.
"""


import argparse
from logging import DEBUG, getLogger
from pathlib import Path

import numpy as np
import pandas as pd
import prody
import pytest

logger = getLogger(__name__)


def calc_contact_rate(dist_mat: np.ndarray, index: int,
                      contact_threshold: int, resnums: np.ndarray) -> float:
    """
    Calculate the contact rate between the residues before index + 1 (start <= i <= index)
    and the residues after index (index + 1 <= i <= end)

    Args:
        dist_mat (np.ndarray): distance matrix (CA atom)
        index (int): index of the residue
        contact_threshold (int): threshold for contact
        resnums (np.ndarray): residue numbers

    Returns:
        float: contact rate
    """
    dist_mat_after_i = dist_mat[:index+1, index + 1:]
    contact_num = np.sum(dist_mat_after_i < contact_threshold)
    if dist_mat_after_i.shape[1] == 0:
        print(index, dist_mat.shape, dist_mat_after_i.shape)
    assert dist_mat_after_i.shape[0] != 0
    assert dist_mat_after_i.shape[1] != 0
    contact_rate = contact_num / (dist_mat_after_i.shape[0] * dist_mat_after_i.shape[1])
    logger.debug(f'{index=}, {resnums[index]=}, {contact_num=}, {contact_rate=:.3f}')
    return contact_rate


def has_low_interaction_region(pdb_file: str, contact_threshold: int = 12,
                               contact_rate_threshold: float = 0.01) -> bool:
    """Judge if the protein has low interaction region

    Args:
        pdb_file (str): path to the pdb file.
        contact_threshold (int, optional): threshold for contact between domain. Defaults to 15.
        contact_rate_threshold (float, optional): threshold for contact rate. Defaults to 0.01.

    Returns:
        bool: True if the protein has low interaction region, False otherwise.
    """
    ca_mol = prody.parsePDB(pdb_file).select('name CA')
    dist_mat = prody.buildDistMatrix(ca_mol)  # build distance matrix
    resindices = np.arange(len(ca_mol))
    resnums = ca_mol.getResnums()
    for i in resindices[: -1]:
        if calc_contact_rate(dist_mat, i, contact_threshold, resnums=resnums) < contact_rate_threshold:
            return True
    return False


@pytest.mark.parametrize("input_pdb, expected", [
    ('../../data/out/dataset/native_pdb/6FF6_A.pdb', True),
    ('../../data/out/dataset/native_pdb/6E3A_A.pdb', False),
    ('../../data/out/dataset/native_pdb/6RUM_A.pdb', False),
    ('../../data/out/dataset/native_pdb/6SC5_A.pdb', True),
    ('../../data/out/dataset/native_pdb/5XSO_A.pdb', False),  # Either True or False
    ('../../data/out/dataset/native_pdb/6RGV_A.pdb', False),  # Either True or False
    ('../../data/out/dataset/native_pdb/6BHF_A.pdb', True),
    ('../../data/out/dataset/native_pdb/6JQ1_A.pdb', True),
    ('../../data/out/dataset/native_pdb/6JDR_A.pdb', True),
    ('../../data/out/dataset/native_pdb/6TJF_C.pdb', True),
    ('../../data/out/dataset/native_pdb/6XW2_A.pdb', False),  # Either True or False
    ('../../data/out/dataset/native_pdb/6HYG_A.pdb', False),  # Either True or False
    ('../../data/out/dataset/native_pdb/6IA5_D.pdb', True),
    ('../../data/out/dataset/native_pdb/6XHV_1P.pdb', True),
    ('../../data/out/dataset/native_pdb/7C2G_G.pdb', True),
    ('../../data/out/dataset/native_pdb/7CF7_A.pdb', True),
    ('../../data/out/dataset/native_pdb/7CHU_A.pdb', True),
    ('../../data/out/dataset/native_pdb/6UVQ_A.pdb', False),
    ('../../data/out/dataset/native_pdb/6VO5_A.pdb', False),
    ('../../data/out/dataset/native_pdb/6FHV_A.pdb', True),
    ('../../data/out/dataset/native_pdb/5XLL_A.pdb', True),
    ('../../data/out/dataset/native_pdb/6QP4_A.pdb', True),
    ('../../data/out/dataset/native_pdb/6WWX_A.pdb', True),
])
def test_judge_low_domain_interaction(input_pdb, expected, caplog):
    caplog.set_level(DEBUG)
    assert has_low_interaction_region(input_pdb) == expected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--target_list', type=str, required=True)
    args = parser.parse_args()
    target_list_path = Path(args.target_list)
    df = pd.read_csv(target_list_path, index_col=0)
    native_dir = Path('../../data/out/dataset/native_pdb')
    count = 0
    for target_id in df['id']:
        native_path = native_dir / f'{target_id}.pdb'
        if has_low_interaction_region(str(native_path)):
            print(target_id)
            count += 1
    print(count)


if __name__ == '__main__':
    main()
