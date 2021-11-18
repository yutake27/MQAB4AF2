"""
Select subsets of the target from target list.
"""

import argparse
import pandas as pd
from pathlib import Path


def get_target_subset(csv_path, how='eq_random', target_num=100, random_state=0) -> pd.DataFrame:
    """
    Select subsets of the target from target list.

    Args:
        csv_path (str): path to the target list.
        how (str): how to select the subset.
            choices: ['random', 'eq_random', 'head', 'tail']
            random: select a random subset of the target list.
            eq_random: select a random subset of the target list with the same size of the target
            that have similar sequence to AF2 traingin dataset and that have not.
            head: select the first target_num targets.
            tail: select the last target_num targets.
        target_num (int): number of target to select.
        random_state (int): random seed.
    Returns:
        pd.DataFrame: subset of the target list.
    """
    df = pd.read_csv(csv_path, index_col=0)
    if how == 'random':
        df_sample = df.sample(n=target_num, random_state=random_state)
    elif how == 'eq_random':
        similar_df = df[df['is_similar_AF2'] == True]
        non_similar_df = df[df['is_similar_AF2'] == False]
        similar_sample = similar_df.sample(target_num // 2, random_state=random_state)
        non_similar_sample = non_similar_df.sample(target_num // 2, random_state=random_state)
        df_sample = pd.concat([similar_sample, non_similar_sample])
    elif how == 'head':
        df_sample = df.head(n=target_num)
    elif how == 'tail':
        df_sample = df.tail(n=target_num)
    else:
        raise ValueError('Invalid how: {}'.format(how))
    return df_sample


def main():
    target_list_path = Path('../../data/interim') / 'target_list.csv'
    parser = argparse.ArgumentParser(description="Select subsets of the target from target list")
    parser.add_argument('-i', '--input_csv_path', type=str, default=target_list_path,
                        help='path to the target list')
    parser.add_argument('--how', type=str, default='eq_random',
                        choices=['random', 'eq_random', 'head', 'tail'],
                        help="""
                        how to select the subset.
                        random: select a random subset of the target list.
                        eq_random: select a random subset of the target list
                        with the same size of the is_similar_AF2.
                        head: select the first target_num targets of the target list.
                        tail: select the last target_num targets of the target list.
                        """)
    parser.add_argument('--target_num', type=int, default=100, help='number of target to select')
    parser.add_argument('--random_state', type=int, default=0)
    parser.add_argument('-o', '--output_csv_path', type=str, default='',
                        help='output path of subset of target list.')
    args = parser.parse_args()

    input_csv_path = args.input_csv_path
    how = args.how
    target_num = args.target_num
    random_state = args.random_state

    if args.output_csv_path == '':
        output_csv_stem = f'target_subset_how_{how}_num_{target_num}'
        if how in ['random', 'eq_random']:
            output_csv_stem += f'_seed_{random_state}'
        output_csv_path = target_list_path.parent / (output_csv_stem + '.csv')
    else:
        output_csv_path = args.output_csv_path
    df_sample = get_target_subset(input_csv_path, how=how, target_num=target_num, random_state=random_state)
    df_sample.to_csv(output_csv_path)


if __name__ == '__main__':
    main()
