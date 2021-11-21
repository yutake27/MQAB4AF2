import argparse
import subprocess
import pandas as pd
from pathlib import Path


def run_sbrod(model_dir):
    output = subprocess.check_output(['sbrod', model_dir])
    return output.decode()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_dir', type=str, help='input directory')
    parser.add_argument('output_csv_path', type=str, help='output csv path')
    args = parser.parse_args()

    input_dir = Path(args.input_dir)

    model_name_array = []
    sbrod_score_array = []
    output = run_sbrod(input_dir / '*.pdb')
    output_list = output.split()
    df = pd.DataFrame({'Model': output_list[::2], 'SBROD': output_list[1::2]})
    df["Model"] = df["Model"].apply(lambda x: Path(x).stem)
    df.to_csv(args.output_csv_path)
