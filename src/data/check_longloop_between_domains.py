"""
Exclude proteins that are connected in long loops between domains from target candidates
(Exclude proteins with long loops and few contact in a region)

Requirement:
- dssp version 3.0.0 (conda install -c salilab dssp)
"""

import subprocess
from functools import reduce
from logging import DEBUG, getLogger
from typing import Any, List, Tuple, Union

import numpy as np
import pandas as pd
import prody
import pytest
import tqdm

logger = getLogger(__name__)


def dssp(pdb_file: str) -> str:
    """run dssp on pdb and return output

    Args:
        pdb_file (str): pdb file

    Returns:
        str: dssp output
    """
    cmd = ['mkdssp', '-i', pdb_file]
    dssp_out = subprocess.run(cmd, capture_output=True, text=True).stdout
    return dssp_out


def parse_dssp_output(dssp_out: str) -> Tuple[List[str], List[str]]:
    """parse dssp output and return a numpy array of secondary structure

    Args:
        dssp_out (str): dssp output

    Returns:
        List[int]: residue numbers
        List[str]: Secondary structure array
    """
    lines = dssp_out.split('\n')
    lines = lines[28: -1]
    judge_misres = lambda line: line[13] == '!'
    resnums = [line[5: 10].strip() for line in lines if not judge_misres(line)]
    ss_array = [line[16] for line in lines if not judge_misres(line)]
    return resnums, ss_array


def get_ss_around_each_residue(ss_array: List[str], padding: int) -> np.ndarray:
    """get secondary structure around each residue

    Args:
        ss_array (List[str]): secondary structure array

    Returns:
        List[str]: secondary structure around each residue
    """
    not_loop = 'X'
    pad_ss_array = [not_loop] * padding + ss_array + [not_loop] * padding
    res_neighbor_ss = np.array([pad_ss_array[i: i + 2 * padding + 1] for i in range(len(ss_array))])
    return res_neighbor_ss


def detect_long_loop(ss_array: List[str], padding: int = 10, loop_ratio: float = 0.7) -> np.ndarray:
    """detect long loop in secondary structure array

    Args:
        ss_array (List[str]): secondary structure array

    Returns:
        List[int]: residue indices of long loop
    """
    ss_around_res = get_ss_around_each_residue(ss_array, padding=padding)
    num_loop_around_res = np.count_nonzero(ss_around_res == ' ', axis=1)
    len_around_res = 2 * padding
    longloop_resindices = np.where(num_loop_around_res >= len_around_res * loop_ratio)[0]
    return longloop_resindices


def calc_num_contact_of_longloop_residues(mol_ca: prody.atomic, ll_resindices: np.ndarray,
                                          threshold_dist: float) -> Union[np.ndarray, Any]:
    """Calculate number of contact of long loop residues

    Args:
        mol_ca (prody.atomic): prody atomic object
        ll_resindices (np.ndarray): long loop residue indices
        threshold_dist (float): threshold distance to be considered as contact

    Returns:
        np.ndarray: number of contact of long loop residues
    """
    dist_mat = prody.buildDistMatrix(mol_ca)
    dist_ll_res = dist_mat[ll_resindices]
    num_res_around_ll_residues = np.count_nonzero(dist_ll_res < threshold_dist, axis=1) - 1  # except itself
    return num_res_around_ll_residues


def has_longloop_between_domains(pdb_file: str, pad_residues: int = 10, threshold_loop_ratio: float = 0.7,
                                 threshold_dist: float = 10, threshold_num_contact: int = 8,
                                 threshold_num_res_low_contact: int = 2) -> np.bool_:
    """check if there is long loop between domains

    Procedure:
    1. run dssp on pdb
    2. detect long loop
    3. calculate contact number of long loop residues and return True
    if contact number is lower than threshold for each long loop residues

    Args:
        pdb_file (str): pdb file
        pad_residues (int): number of padding residues around central residue to be considered
        threshold_loop_ratio (float): percentage of loop to judge as long loop
        threshold_dist (float): threshold distance to be considered as contact
        threshold_num_contact (int): threshold number of contact to judge long loop between domains
        threshold_num_res_low_contact (int):
        threshold number of residues to be considered as long loop with low contact

    Returns:
        bool: True if there is long loop between domains
    """
    logger.info(f'Checking long loop between domains in {pdb_file}')
    dssp_out = dssp(pdb_file)
    logger.debug(dssp_out)
    resnums, ss_array = parse_dssp_output(dssp_out)
    logger.debug(f'{resnums=}')
    logger.debug(f'{ss_array=}')
    mol = prody.parsePDB(pdb_file)
    mol_ca = mol.select('name CA resnum {}'.format(reduce(lambda a, b: a + ' ' + b, resnums)))
    logger.debug(f'{mol_ca.getResnums()=}')
    assert len(mol_ca) == len(ss_array)
    longloop_resindices = detect_long_loop(ss_array, padding=pad_residues, loop_ratio=threshold_loop_ratio)
    logger.debug(f'{longloop_resindices=}')
    longloop_resnums = np.array(resnums)[longloop_resindices]
    logger.info(f'{longloop_resnums=}')
    num_contact_of_ll_residues = calc_num_contact_of_longloop_residues(
        mol_ca, longloop_resindices, threshold_dist
    )
    logger.info(f'{num_contact_of_ll_residues=}')
    assert len(longloop_resindices) == len(num_contact_of_ll_residues)
    low_contact_indices = np.where(num_contact_of_ll_residues <= threshold_num_contact)[0]
    resnums_low_contact = np.array(longloop_resnums)[low_contact_indices]
    logger.info(f'{resnums_low_contact=}')
    has_ll_with_low_contact = np.sum(
        num_contact_of_ll_residues <= threshold_num_contact
    ) >= threshold_num_res_low_contact
    return has_ll_with_low_contact


@pytest.mark.parametrize(('pdb_file', 'excluded'), [
    ('../../data/out/dataset/native_pdb/6IA5_D.pdb', True),
    ('../../data/out/dataset/native_pdb/6XHV_1P.pdb', True),
    ('../../data/out/dataset/native_pdb/7C2G_G.pdb', True),
    ('../../data/out/dataset/native_pdb/7CF7_A.pdb', False),
    ('../../data/out/dataset/native_pdb/6UVQ_A.pdb', False),
    ('../../data/out/dataset/native_pdb/6VO5_A.pdb', False),
    ('../../data/out/dataset/native_pdb/6JQ1_A.pdb', True),
    ('../../data/out/dataset/native_pdb/6TJF_C.pdb', True),
    ('../../data/out/dataset/native_pdb/6FF6_A.pdb', True)
])
def test_exclude_target_with_longloop_between_domain(pdb_file, excluded, caplog):
    caplog.set_level(DEBUG)
    """test long loop detection"""
    be_excluded = has_longloop_between_domains(pdb_file)
    assert be_excluded == excluded


if __name__ == '__main__':
    target_csv = '../../data/interim/target_subset_how_eq_random_num_300_seed_0.csv'
    df = pd.read_csv(target_csv, index_col=0)
    excluded_targets = []
    for i, row in tqdm.tqdm(df.iterrows(), total=len(df)):
        target_id = row['id']
        num_entry_in_cluster_af2_not_include = row['num_entry_in_cluster_AF2_notInclude']
        pdb_file = f'../../data/out/dataset/native_pdb/{target_id}.pdb'
        try:
            has_ll_with_low_contact = has_longloop_between_domains(pdb_file)
        except AssertionError:
            print('AssertionError', pdb_file)
            break
        if has_ll_with_low_contact:
            excluded_targets.append((target_id, num_entry_in_cluster_af2_not_include))
    print(len(excluded_targets))
    print(f'{excluded_targets=}')
