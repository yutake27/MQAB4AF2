import argparse
from pathlib import Path
from typing import List

import pandas as pd

data_dir = Path('../../../data/out/dataset/')
score_dir = data_dir / 'score'
unrelaxed_label_dir = score_dir / 'label'
relaxed_label_dir = score_dir / 'label_relaxed'
subset_score_dir = score_dir / 'subsets'


def concat_labels(label_dir: Path, target_list: List, skip: bool = False) -> pd.DataFrame:
    label_dfs = []
    for target in target_list:
        label_csv = label_dir / f'{target}.csv'
        if not label_csv.exists():
            if skip:
                print(f'Label for {target} does not exist. Skipping.')
                continue
            raise FileNotFoundError(f'{label_csv} does not exist')
        label_df = pd.read_csv(label_csv, index_col=0)
        label_df['Target'] = target
        label_df = label_df.sort_values('Model').reset_index(drop=True)
        label_dfs.append(label_df)
    return pd.concat(label_dfs).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser("Concatenate labels of targets in target list")
    parser.add_argument('target_list_csv', type=str, help='Target list')
    parser.add_argument('-r', '--relaxed', action='store_true', help='Concat labels for relaxed')
    parser.add_argument('-o', '--output_csv_path', type=str, help='Output file name')
    parser.add_argument('-s', '--skip', action='store_true', help='Skip not existing labels')
    args = parser.parse_args()

    target_list_csv = Path(args.target_list_csv)
    subset_name = target_list_csv.stem
    label_dir = unrelaxed_label_dir if not args.relaxed else relaxed_label_dir
    target_list = pd.read_csv(target_list_csv, index_col=0)['id'].tolist()
    label_df = concat_labels(label_dir, target_list, args.skip)
    if args.output_csv_path:
        output_csv_path = Path(args.output_csv_path)
    else:
        output_csv_name = 'label.csv' if not args.relaxed else 'label_relaxed.csv'
        output_csv_path = subset_score_dir / subset_name / output_csv_name
    output_csv_path.parent.mkdir(parents=True, exist_ok=True)
    label_df.to_csv(output_csv_path)


if __name__ == "__main__":
    main()
