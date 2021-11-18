import numpy as np
from Bio import SeqIO, pairwise2
from prody import parsePDB
import prody.atomic.atomgroup


class ReadSeq:

    @staticmethod
    def mol2seq(mol: prody.atomic.atomgroup.AtomGroup,
                insert_gap: bool = True) -> (str, np.ndarray):
        """Convert mol (prody.atomic.atomgroup.AtomGroup) to sequence.

        Args:
            mol (prody.atomic.atomgroup.AtomGroup): [description].
            insert_gap (bool, optional): Flag to insert gaps in areas of missing residue. Defaults to True.

        Returns:
            str: Sequece of mol.
            np.ndarray: Array of the residue numbers.
        """
        resindices = mol.getResindices()
        resindices_diff = np.diff(resindices, prepend=float('inf'))
        res_start_indices = np.where(resindices_diff != 0)[0]
        seq = ''.join(np.array(list(mol.getSequence()))[res_start_indices])
        resnums = mol.getResnums()[res_start_indices]
        assert len(resnums) == len(np.unique(resnums))
        if insert_gap:
            min_resnums = np.min(resnums)
            len_seq = np.max(resnums) - min_resnums + 1
            seq_array = np.array(['-'] * len_seq)
            seq_indices = resnums - min_resnums
            seq_array[seq_indices] = list(seq)
            seq = ''.join(seq_array)
            resnums = np.where(seq_array != '-')[0] + min_resnums
        return seq, resnums

    @classmethod
    def pdb2seq(cls,
                pdb_path: str,
                chain: str = None,
                read_HETATM: bool = False,
                insert_gap: bool = True) -> (str, np.ndarray):
        """Convert pdb to sequence.

        Args:
            pdb_path (str): Path to the pdb file.
            chain (str, optional): Chain name. Defaults to None.
            read_HETATM (bool, optional): Flag for parsing HETATM or not.
            If True, parse HETATM, else dismiss HETATM.
            Defaults to False.
            insert_gap (bool, optional): Flag to insert gaps in areas of missing residue. Defaults to True.

        Returns:
            str: Sequence of the pdb.
            np.ndarray: Array of the residue numbers.
        """
        if chain:
            mol = parsePDB(pdb_path, chain=chain)
        else:
            mol = parsePDB(pdb_path)
        if not read_HETATM:
            mol = mol.select('not hetatm')
        seq, resnums = cls.mol2seq(mol, insert_gap=insert_gap)
        return seq, resnums

    @staticmethod
    def fasta2seq(fasta_path: str) -> str:
        """Convert fasta to sequence.

        Args:
            fasta_path (str): Path to the fasta file.

        Returns:
            str: Sequence of the fasta.
        """
        seq_record = list(SeqIO.parse(fasta_path, 'fasta'))[0]
        seq = str(seq_record.seq)
        return seq


class AlignSeq:

    @staticmethod
    def _get_max_serial_position(np_array: np.ndarray) -> (int, int):
        """Get the longest consecutive residue number after alignment.

        Args:
            np_array (np.ndarray): Residue indices of template structure after alignment.

        Returns:
            int: index of the start residue.
            int: index of the end residue.
        """
        start, end = 0, 0
        max_serial_len, tmp_len = 0, 0
        b_num = -2
        for i, num in enumerate(np_array):
            if num == b_num + 1:
                tmp_len += 1
            else:
                tmp_len = 1
                tmp_start = i
            b_num = num
            tmp_end = i
            if max_serial_len < tmp_len:
                start = tmp_start
                end = tmp_end
                max_serial_len = tmp_len
        return start, end

    @classmethod
    def align_seq(cls,
                  seqA: str,
                  seqB: str,
                  match: int = 5,
                  mismatch: int = 1,
                  gap_start: int = -2,
                  gap_extend: int = -1) -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray):
        """align two sequence using biopython pairwise2.
        If there are multiple alignments, use the one with the smallest alignment variance.

        Args:
            seqA (str): Sequence A.
            seqB (str): Sequence B.

        Returns:
            np.ndarray: Sequence A after alignment.
            np.ndarray: Sequence B after alignment.
            np.ndarray: Indices of sequence A after alignment.
            np.ndarray: Indices of sequence B after alignment.
        """
        alignment_list = pairwise2.align.globalms(seqA, seqB, match, mismatch, gap_start, gap_extend)
        alignment_length = -1
        for alignment in alignment_list:
            align_seqA, align_seqB = np.array(list(alignment.seqA)), np.array(list(alignment.seqB))
            align_seqA, align_seqB = align_seqA[np.where(align_seqB != '-')], align_seqB[np.where(align_seqA != '-')]
            align_indicesA = np.where(align_seqB != '-')[0]
            align_indicesB = np.where(align_seqA != '-')[0]
            if len(align_indicesA) > alignment_length:
                alignment_length = len(align_indicesA)
                ret_align_seqA = align_seqA
                ret_align_seqB = align_seqB
                ret_align_indicesA = align_indicesA
                ret_align_indicesB = align_indicesB
        return ret_align_seqA, ret_align_seqB, ret_align_indicesA, ret_align_indicesB

    @classmethod
    def align_fasta_and_pdb(cls,
                            fasta_seq: str,
                            pdb_seq: str,
                            alignment_percentage_threshold: float = 0.7)\
            -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray):
        """align sequence of fasta and sequence of pdb using biopython pairwise2.
        If there are multiple alignments, use the one with the smallest alignment variance.
        Acquire only the positions where the alignment is continuous in order to prevent incorrect alignment.

        Args:
            fasta_seq (str): Sequence of fasta.
            pdb_seq (str): Sequence of pdb.
            alignment_percentage_threshold (float):
            Threshold value of the percentage of alignment length that causes alignment errors.
            If the percentage of the alignment length is lower than this threshold, a ValueError is raised.

        Returns:
            np.ndarray: Sequence of fasta after alignment.
            np.ndarray: Sequence of pdb after alignment.
            np.ndarray: Indices of fasta after alignment.
            np.ndarray: Indices of pdb after alignment.
        """
        align_fseq, align_pseq, align_findices, align_pindices = cls.align_seq(fasta_seq, pdb_seq)
        p_start, p_end = cls._get_max_serial_position(align_pindices)
        alignment_length = p_end - p_start + 1
        if alignment_length < len(fasta_seq) * alignment_percentage_threshold or alignment_length < 30:
            raise ValueError('Alignment length is not sufficient.' +
                             f'Only match {alignment_length} residues in total {len(fasta_seq)} residues.')
        p_not_serial = np.delete(np.arange(len(align_pindices)), np.arange(p_start, p_end + 1), 0)
        align_findices = np.delete(align_findices, p_not_serial, 0)
        align_pindices = np.delete(align_pindices, p_not_serial, 0)
        align_pseq[p_not_serial] = '-'
        return align_fseq, align_pseq, align_findices, align_pindices
