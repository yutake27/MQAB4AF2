"""
    Select targets from pdb entries
"""
import datetime
import inspect
import json
import sqlite3
import sys
import time
from pathlib import Path
from typing import Tuple, Union

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

sys.path.append('../alphafold/util/')
from download_native_pdb import DownloadProtein

data_dir = Path('../../data')


def now():
    return datetime.datetime.now()


class SearchAPI:
    """Search entries using SearchAPI of PDB
    """
    api_url = 'https://search.rcsb.org/rcsbsearch/v1/query?json='

    @classmethod
    def get(cls, query: json) -> dict:
        res = requests.get(cls.api_url + json.dumps(query))
        assert res.status_code == 200
        return res.json()

    @classmethod
    def get_entries(cls, query: json) -> list:
        res_json = cls.get(query)
        entry_list = res_json['result_set']
        return entry_list


def get_entry() -> list:
    with open('./specific_PDB_entry_query.json', 'r') as f:
        query = json.load(f)
        entry_list = SearchAPI.get_entries(query)
        return entry_list


class GraphQL:
    """Fetch entry information of PDB using data-api
    """
    api_url = 'https://data.rcsb.org/graphql?query='

    @classmethod
    def get(cls, query: str) -> dict:
        res = requests.get(cls.api_url + query)
        assert res.status_code == 200
        return res.json()

    @classmethod
    def post(cls, query: str, variables: dict) -> dict:
        res = requests.post(cls.api_url, json={'query': query, 'variables': variables})
        assert res.status_code == 200
        return res.json()

    @classmethod
    def fetch_resolution_and_releasedate(cls, pdb_ids: list) -> Tuple[list, list]:
        query = ('query($ids: [String!]!){entries(entry_ids: $ids){rcsb_entry_info{resolution_combined},' +
                 'rcsb_accession_info{initial_release_date}}}')
        variables = {'ids': pdb_ids}
        res_json = GraphQL.post(query, variables)
        resolutions = [res[0] if (res := ent['rcsb_entry_info']['resolution_combined']) is not None else None
                       for ent in res_json['data']['entries']]
        releasedates = [ent['rcsb_accession_info']['initial_release_date']
                        for ent in res_json['data']['entries']]
        assert len(pdb_ids) == len(resolutions)
        assert len(pdb_ids) == len(releasedates)
        return resolutions, releasedates


class SearchSequence:
    pdb_seqres_file = data_dir / 'raw/pdb_seqres.txt'

    def __init__(self, entries):
        entry_set = set(entries)
        self.headers, self.sequences = self._prefiltering(entry_set)

    def _prefiltering(self, entry_set: set) -> Tuple[list, list]:
        with open(self.pdb_seqres_file, 'r') as f:
            headers = []
            sequences = []
            while True:
                header = f.readline().strip()
                seq = f.readline().strip()
                if not header:
                    break
                if header[1: 5].upper() in entry_set:
                    headers.append(header)
                    sequences.append(seq)
            return headers, sequences

    def search_sequence(self, id: str) -> Tuple[str, str]:
        """Search sequence from pdb seqres.

        Args:
            id (str): 'pdb_id' + '_' + 'chain'. Upper case.

        Returns:
            [str]: fasta header.
            [str]: fasta sequence.
        """
        lines = []
        for header, sequence in zip(self.headers, self.sequences):
            # header example: '>6arp_D mol:protein length:222 ...'
            header_head = header.split()[0][1:]
            header_id = header_head[: 4].upper() + header_head[4:]
            if header_id == id:
                lines.append(header)
                lines.append(sequence)
        assert len(lines) == 2
        header, sequence = lines
        return header, sequence


class CheckSequence:
    @classmethod
    def check_sequence(cls, sequence: str) -> bool:
        for method_name, func in inspect.getmembers(cls):
            if method_name.startswith('_check') and not func(sequence):
                return False
        return True

    @staticmethod
    def _check_sequence_length(sequence) -> bool:
        if len(sequence) < 80:  # min (determined from AF2 paper)
            return False
        elif len(sequence) > 700:  # max (determined from the sequence length distribution)
            return False
        else:
            return True

    @staticmethod
    def _check_sequence_character(sequence) -> bool:
        set_sequence_character = set(list(sequence))
        # check character number
        if len(set_sequence_character) < 5:
            return False
        # check all amino acids are standard amino acid
        standard_amino_acids = ('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
                                'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
        for amino_acid in set_sequence_character:
            if amino_acid not in standard_amino_acids:
                return False
        return True


def make_fasta(fasta_path: Path, header: str, sequence: str, len_line=100):
    with open(fasta_path, 'w') as f:
        f.write(header + '\n')
        for i in range((len(sequence) - 1) // len_line + 1):
            f.write(sequence[i * len_line: (i + 1) * len_line] + '\n')


class CheckNative:
    @classmethod
    def check_native(cls, id: str, header: str, sequence: str) -> bool:
        pdb_path = data_dir / 'out' / 'dataset' / 'native_pdb' / f'{id}.pdb'
        pdb_id, chain = id.split('_')
        fasta_path = data_dir / 'out' / 'dataset' / 'fasta' / f'{id}.fasta'

        def _remove_files():
            pdb_path.unlink()
            fasta_path.unlink()

        if not fasta_path.exists():  # make fasta file
            make_fasta(fasta_path, header, sequence)
        if not pdb_path.exists():  # download pdb file
            try:
                time.sleep(1)
                DownloadProtein.download_pdb_specific_chain(pdb_id, chain, fasta_path, pdb_path)
            except (ValueError, AssertionError):
                fasta_path.unlink()
                return False
        # Check the ratio of missing residues
        try:
            num_diff, num_missing, length = DownloadProtein.test_match_pdb_and_fasta(pdb_path, fasta_path)
        except AssertionError:
            _remove_files()
            return False
        if (num_missing / length) > 0.2:
            _remove_files()
            return False
        return True


class Entry:
    """A class that stores information about an entry
    and determines if it is a suitable for a target.
    """

    def __init__(self, id: str, resolution: float, releasedate: str,
                 header: str, sequence: str, num_entry_in_cluster: int,
                 num_entry_in_cluster_AF2_notInclude: int):
        """Init

        Args:
            id (str): entry id. pdb_id(upper case) + '_' + chain_name.
            resolution (float): resolution.
            releasedate (str): Initial release date of entry.
            header (str): Header of fasta sequence.
            sequence (str): Amino acid sequence.
            num_entry_in_cluster (int): Number of entries in the cluster clustered by sequence identity 40%.
            num_entry_in_cluster_AF2_notInclude (int):
            Number of entries in the cluster that are not included in AF2 training (After 2018-05-01).
        """
        self.id = id
        self.resolution = resolution
        self.releasedate = releasedate
        self.header = header
        self.sequence = sequence
        self.length = len(sequence)
        self.num_entry_in_cluster = num_entry_in_cluster
        self.num_entry_in_cluster_AF2_notInclude = num_entry_in_cluster_AF2_notInclude
        # Whether similar sequences are included in the training data of AF2
        self.is_similar_AF2 = num_entry_in_cluster != num_entry_in_cluster_AF2_notInclude

    def _check_sequence(self) -> bool:
        return CheckSequence.check_sequence(self.sequence)

    def _check_resolution(self):
        pass

    def _check_native(self):
        return CheckNative.check_native(self.id, self.header, self.sequence)

    def check_entry(self) -> bool:
        return self._check_sequence() and self._check_native()

    def get_all_attribute(self) -> dict:
        return {key: value for key, value in self.__dict__.items()}


class Cluster:
    bc40 = data_dir / 'raw/bc-40.out'

    def __init__(self, entries: list):
        cluster_list = []
        with open(self.bc40, 'r') as f:
            for line in f.readlines():
                cluster_list.append(line.split())
        self.origin_clusters: list = cluster_list
        self.clusters: list = self._clustering(entries)
        self.ss = SearchSequence(entries)
        resolutions, releasedates = self._pre_fetch_resolution_and_releasedate(entries)
        self.resolution_dict = {ent: resolution for ent, resolution in zip(entries, resolutions)}
        self.releasedates_dict = {ent: releasedate for ent, releasedate in zip(entries, releasedates)}

    def _clustering(self, entries: list) -> list:
        """Split entries into cluster based on bc-40
        (sequence clusters published by rcsb).

        Args:
            entries (list): List of entries.

        Returns:
            list: List of clusters.
        """
        entry_set = set(entries)
        extracted_cluster_list = []
        for cluster in self.origin_clusters:
            extracted_cluster = []
            for ent in cluster:
                if ent[: 4] in entry_set:
                    extracted_cluster.append(ent)
            if len(extracted_cluster) > 0:
                extracted_cluster_list.append([extracted_cluster, len(cluster)])
        return extracted_cluster_list

    def _pre_fetch_resolution_and_releasedate(self, pdb_ids: list) -> Tuple[list, list]:
        """Fetch the resolution and release date of all entries before clustering in advance

        Args:
            pdb_ids (list): List of pdb id.

        Returns:
            list: List of resolution.
            list: List of release date.
        """
        print(now(), 'Fetch resolution and release date')
        resolutions, releasedates = GraphQL.fetch_resolution_and_releasedate(pdb_ids)
        print(now(), 'Finish fetching')
        return resolutions, releasedates

    def _get_resolution_and_releasedate_from_prefetched(self, pdb_ids: list) -> Tuple[list, list]:
        """Get resolution and release date from pre-fetched information.
        """
        resolutions = [self.resolution_dict[pdb_id] for pdb_id in pdb_ids]
        releasedates = [self.releasedates_dict[pdb_id] for pdb_id in pdb_ids]
        return resolutions, releasedates

    def _select_target_from_cluster(self, entries: list, num_entry_in_cluster: int) -> Union[dict, None]:
        """Select the target from the cluster.
        Retrieve the entries from the cluster in the order of best resolution
        and check if they are suitable for the target.

        Args:
            entries (list): List of entry.
            num_entry_in_cluster (int): Total number of entries in the cluster.

        Returns:
            dict: Dict of target information (entry_id, resolution, sequence, size of cluster.)
            If there is no suitable entry, return None.
        """
        pdb_ids = [ent[: 4] for ent in entries]
        resolutions, releasedates = self._get_resolution_and_releasedate_from_prefetched(pdb_ids)
        indices = np.argsort(resolutions)
        for index in indices:
            entry_id, resolution, releasedate = entries[index], resolutions[index], releasedates[index]
            header, sequence = self.ss.search_sequence(entry_id)
            entry = Entry(entry_id, resolution, releasedate, header, sequence, num_entry_in_cluster, len(entries))
            if entry.check_entry():
                return entry.get_all_attribute()
        return None

    def _get_target_info_by_name(self, target_name: str, num_entry_in_cluster: int,
                                 num_entry_in_cluster_AF2_notInclude: int) -> dict:
        """Get target information by name.
        """
        pdb_id = target_name[: 4]
        resolution = self.resolution_dict[pdb_id]
        releasedate = self.releasedates_dict[pdb_id]
        header, sequence = self.ss.search_sequence(target_name)
        entry = Entry(target_name, resolution, releasedate, header, sequence,
                      num_entry_in_cluster, num_entry_in_cluster_AF2_notInclude)
        return entry.get_all_attribute()

    def select_targets_from_each_cluster(self) -> list[dict]:
        """Select targets from each cluster.

        Returns:
            list: List of dict of target information.
        """
        targets: list[dict] = []
        db_path = data_dir / 'interim' / 'db' / 'targets.db'
        db_path.parent.mkdir(parents=True, exist_ok=True)
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        cur.execute('CREATE TABLE IF NOT EXISTS targets(cluster_index INT, id TEXT)')
        con.commit()
        for i, (entries, num_entry_in_cluster) in enumerate(tqdm(self.clusters, desc='Select target from cluster')):
            cache_exist = cur.execute('SELECT * FROM targets WHERE cluster_index = ?', (i,)).fetchone()
            if cache_exist:
                target_name = cache_exist[1]
                if target_name is None:
                    continue
                target = self._get_target_info_by_name(target_name, num_entry_in_cluster, len(entries))
                targets.append(target)
            else:
                target = self._select_target_from_cluster(entries, num_entry_in_cluster)
                if target is None:
                    cur.execute("""INSERT INTO targets VALUES(?, ?)""", (i, None))
                else:
                    targets.append(target)
                    cur.execute("""INSERT INTO targets VALUES(?, ?)""", (i, target['id']))
                con.commit()
        con.close()
        return targets


def main():
    print(now(), 'Fetch specific entries')
    entry_list = get_entry()
    print(now(), 'Clustering')
    cluster = Cluster(entry_list)
    print(now(), 'Select targets from each cluster')
    targets = cluster.select_targets_from_each_cluster()
    print(len(targets))
    target_df = pd.DataFrame(targets)
    output_csv = data_dir / 'interim' / 'target_list.csv'
    target_df.to_csv(output_csv)


if __name__ == '__main__':
    main()
