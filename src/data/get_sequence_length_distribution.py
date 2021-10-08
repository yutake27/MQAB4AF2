import json

import pandas as pd

from select_target import SearchAPI, SearchSequence


def get_entry():
    with open('./specific_PDB_entry_query.json', 'r') as f:
        query = json.load(f)
        entry_list = SearchAPI.get_entries(query)
        return entry_list


def get_sequence_lengths_from_header(headers: list) -> list:
    lengths = []
    for header in headers:
        header_split = header.split()
        if header_split[1] == 'mol:protein':
            length = int(header_split[2].split(':')[-1])
            lengths.append(length)
    return lengths


def main():
    entry_list = get_entry()
    ss = SearchSequence(entry_list)
    headers = ss.headers
    lengths = get_sequence_lengths_from_header(headers)
    print('Num sequences:', len(lengths))
    sorted_lengths = sorted(lengths)
    lower_limit = 0.025
    upper_limit = 0.975
    print('2.5%: {}, 97.5%: {}'.format(
        sorted_lengths[int(len(lengths) * lower_limit)],
        sorted_lengths[int(len(lengths) * upper_limit)]
    ))
    print('min: {}, max: {}'.format(
        sorted_lengths[0],
        sorted_lengths[-1]
    ))
    df = pd.DataFrame(sorted_lengths, columns=['Length'])
    df.to_csv('../../data/interim/sequence_length.csv')


if __name__ == '__main__':
    main()
