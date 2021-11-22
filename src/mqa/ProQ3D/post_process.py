"""
Convert ProQ3D output to csv
"""

import argparse
from pathlib import Path
import pandas as pd


def parse_score(txt):
    with open(txt, 'r') as f:
        lines = f.readlines()
    score = lines[1].split()
    return score


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('score_dir', type=str, help='proq3 score dir for one target')
    args = parser.parse_args()

    proq3_score_dir = Path(args.score_dir)

    result_array = []
    for txt in proq3_score_dir.glob('*.global'):
        model_name = str(txt.name).split('.')[0]
        score = parse_score(txt)
        result_array.append([model_name] + score)
    df = pd.DataFrame(result_array, columns=['Model', 'ProQ2D', 'ProQRosCenD', 'ProQRosFAD', 'ProQ3D'])

    df.to_csv(proq3_score_dir.with_suffix('.csv'))
