"""
Get global scores from local scores.
"""

import argparse
from pathlib import Path

import pandas as pd


def get_global_score(csv_path: str):
    df = pd.read_csv(csv_path, sep='\t')
    global_score = df['prediction'].mean()
    return global_score


def get_df(dir_path: Path) -> pd.DataFrame:
    model_name_list = []
    global_score_list = []

    for csv_path in dir_path.glob('*.scores'):
        global_score = get_global_score(str(csv_path))
        model_name_list.append(csv_path.stem)
        global_score_list.append(global_score)
    df = pd.DataFrame({'Model': model_name_list, 'VoroCNN': global_score_list})
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('score_dir', type=str, help='Directory that contains *.scores.')
    parser.add_argument('-o', '--out', type=str,
                        help='Output csv path. If unspecified, score_dir.csv will be made.')
    args = parser.parse_args()

    score_dir = Path(args.score_dir)
    df = get_df(score_dir)
    df = df.sort_values('Model').reset_index(drop=True)
    if args.out is None:
        output_path = score_dir.with_suffix('.csv')
    else:
        output_path = args.out
    df.to_csv(output_path)
