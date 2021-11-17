"""
Format the output file of P3CMQA
"""


import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser('Format the output file of P3CMQA')
    parser.add_argument('output_score_dir', type=str, help='Output directory of P3CMQA')
    parser.add_argument('output_global_csv_path', type=str, help='Path to global score of P3CMQA.')
    args = parser.parse_args()

    # Convert txt of local score to csv
    output_score_dir = Path(args.output_score_dir)
    for txt in output_score_dir.glob('*.txt'):
        df = pd.read_csv(txt, sep='\t', header=2)
        output_csv_path = txt.with_suffix('.csv')
        df.to_csv(output_csv_path)

    # Format global score csv
    global_score_df = pd.read_csv(args.output_global_csv_path)
    global_score_df = global_score_df.rename(columns={'Model_name': 'Model', 'Score': 'P3CMQA'})
    global_score_df = global_score_df.sort_values('Model').reset_index(drop=True)
    output_global_score_path = output_score_dir.with_suffix('.csv')
    global_score_df.to_csv(output_global_score_path)


if __name__ == "__main__":
    main()
