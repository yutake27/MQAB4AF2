import gzip
import os
import tempfile
from typing import Tuple
import urllib.request
import warnings
from functools import reduce
from pathlib import Path

import numpy as np
import prody.atomic.atomgroup
from Bio.PDB import PDBIO, MMCIFParser, PDBExceptions
from prody import parsePDB, writePDB

try:
    from .seq import AlignSeq, ReadSeq
except ImportError:
    from seq import AlignSeq, ReadSeq

warnings.simplefilter('ignore', PDBExceptions.PDBConstructionWarning)


class DownloadProtein:
    rcsb_download_url = 'https://files.rcsb.org/download/'

    @classmethod
    def download_native_pdb(cls, pdb_id: str, chain: str, output_path: str):
        """Download the pdb file. If the pdb format is not available, download the mmcif format,
        select only the specified chain, and convert it to pdb.

        Args:
            pdb_id (str): PDB ID.
            chain (str): Chain name.
            output_path (str): Path to the output file.
        """
        try:
            cls.download_pdb(pdb_id, output_path)
        except urllib.request.HTTPError:
            cls.download_mmcif(pdb_id, output_path)
            cls.mmcif2pdb(str(output_path), chain, str(output_path))

    @staticmethod
    def _download_and_decompress_file_from_url(url: str, output_path: str):
        urllib.request.urlretrieve(url, output_path)
        with gzip.open(output_path, 'rt') as f:
            lines = f.readlines()
        with open(output_path, 'w') as f:
            f.writelines(lines)

    @classmethod
    def download_pdb(cls, pdb_id: str, output_path: str = None):
        pdb_url = cls.rcsb_download_url + pdb_id.upper() + '.pdb.gz'
        if output_path is None:
            output_path = pdb_id.upper() + '.pdb'
        cls._download_and_decompress_file_from_url(pdb_url, output_path)

    @classmethod
    def download_mmcif(cls, pdb_id: str, output_path: str = None):
        mmcif_url = cls.rcsb_download_url + pdb_id.upper() + '.cif.gz'
        if output_path is None:
            output_path = pdb_id.upper() + '.cif'
        cls._download_and_decompress_file_from_url(mmcif_url, output_path)

    @classmethod
    def mmcif2pdb(cls, input_mmcif_path: str, chain: str, output_pdb_path: str) -> str:
        parser = MMCIFParser()
        structure = parser.get_structure('', input_mmcif_path)
        try:
            sel_structure = structure[0][chain]
        except KeyError:
            print('chain list:', sorted([chain.id for chain in structure.get_chains()]))
            raise KeyError('The specified chain "{}" does not exist in template mmcif'.format(chain))

        if len(chain) > 1:
            if chain[0] in [chain.id for chain in structure.get_chains()]:
                structure[0][chain[0]].id = None
            sel_structure.id = chain[0]
        io = PDBIO()
        io.set_structure(sel_structure)
        io.save(output_pdb_path)

    @staticmethod
    def correct_resnums(mol: prody.atomic.atomgroup.AtomGroup) -> np.ndarray:
        """correct the residue number if different residues are assigned the same residue number.

        Args:
            mol (prody.Atom.atomic): mol read by prody.

        Returns:
            np.ndarray: residue numbers after modification.
        """
        new_resnum_list = []
        plus = 0
        resnum_before = float('inf')
        resid_before = float('inf')
        for resnum, resid in zip(mol.getResnums(), mol.getResindices()):
            if resnum_before == resnum and resid_before != resid:
                plus += 1
            new_resnum_list.append(resnum + plus)
            resnum_before = resnum
            resid_before = resid
        return np.array(new_resnum_list)

    @staticmethod
    def _test_resnum_match(fasta_seq: str, mol: prody.atomic.atomgroup.AtomGroup) -> Tuple[int, int]:
        """test that the fasta sequence matches to the pdb sequence.

        Args:
            fasta_seq (str): Sequence of the fasta.
            mol (prody.atomic.atomgroup.AtomGroup): PDB object read by ProDy.
        """
        pdb_seq, pdb_resnum = ReadSeq.mol2seq(mol, insert_gap=False)
        fasta_seq_array = np.array(list(fasta_seq))
        pdb_seq_array = np.copy(fasta_seq_array)
        pdb_seq_array[pdb_resnum - 1] = list(pdb_seq)
        num_diff = np.count_nonzero(fasta_seq_array != pdb_seq_array)
        num_missing = len(fasta_seq) - len(pdb_seq)
        assert num_diff < len(fasta_seq) * 0.05
        print('length:', len(fasta_seq))
        print('num different residues between pdb and fasta:', num_diff)
        print('num missing residues:', num_missing)
        return num_diff, num_missing

    @staticmethod
    def test_match_pdb_and_fasta(pdb_path: str, fasta_path: str) -> Tuple[int, int, int]:
        """test that the pdb sequence matches to the fasta sequence.

        Args:
            pdb_path (str): Path to the pdb file.
            fasta_path (str): Path to the fasta file.

        Returns:
            Tuple[int, int, int]: Number of different residues between pdb and fasta,
            number of missing residues, and length of the target sequence.
        """
        mol = parsePDB(pdb_path)
        fasta_seq = ReadSeq.fasta2seq(fasta_path)
        num_diff, num_missing = DownloadProtein._test_resnum_match(fasta_seq, mol)
        return num_diff, num_missing, len(fasta_seq)

    @classmethod
    def download_pdb_specific_chain(cls, pdb_id: str, chain: str, fasta_path: str,
                                    output_pdb_path: str, alignment_percentage_threshold=0.3) -> None:
        """Download the specific chain of pdb and fix residue numbers.

        Args:
            pdb_id (str): PDB ID.
            chain (str): Chain name.
            fasta_path (str): Path to the fasta file of the pdb_id.
            output_pdb_path (str): Path to the output pdb file.
        """
        with tempfile.NamedTemporaryFile() as f:
            # tmp_pdb_path = pdb_id + '.pdb'
            tmp_pdb_path = f.name
            cls.download_native_pdb(pdb_id, chain, tmp_pdb_path)
            # fasta_seq = ReadSeq.fasta2seq(self.native_fasta_path)
            mol = parsePDB(tmp_pdb_path, chain=chain)
            # if the first residue number is negative, correct the residue numbers
            if mol.getResnums()[0] < 0:
                resnums = mol.getResnums()
                new_resnums = resnums - resnums[0] + 1
                mol.setResnums(new_resnums)
            new_resnums = cls.correct_resnums(mol)
            mol.setResnums(new_resnums)
            pdb_seq, pdb_resnums = ReadSeq.mol2seq(mol, insert_gap=True)
            fasta_seq = ReadSeq.fasta2seq(fasta_path)
            align_pseq, align_fseq, align_findices, align_pindices = AlignSeq.align_fasta_and_pdb(
                fasta_seq, pdb_seq, alignment_percentage_threshold=alignment_percentage_threshold)
            sel_mol_resnums = pdb_resnums[align_pindices]
            sel_mol = mol.select('resnum {}'.format(reduce(lambda a, b: str(a) + ' ' + str(b), sel_mol_resnums)))
            assert sel_mol is not None
            fasta_resnum = align_findices + 1
            convert_resnum_dict = dict(zip(sel_mol_resnums, fasta_resnum))
            new_resnum = [convert_resnum_dict[resnum] for resnum in sel_mol.getResnums()]
            sel_mol.setResnums(new_resnum)
            cls._test_resnum_match(fasta_seq, sel_mol)
            writePDB(str(output_pdb_path), sel_mol)
