'''
Calculate the accuracy of the relaxed structure.
'''

import argparse
import sys
sys.path.append('..')
from get_model_accuracy import ModelAccuracy


def main():
    parser = argparse.ArgumentParser(description='Calculate the accuracy of the relaxed structure.')
    parser.add_argument('-i', '--input_pdb_dir', help='Input directory', required=True)
    parser.add_argument('-n', '--native_pdb', help='Native pdb', required=True)
    parser.add_argument('-o', '--output_csv_path', help='Output file', required=True)
    args = parser.parse_args()

    # Calculate the accuracy of the relaxed structure.
    tmscore_df = ModelAccuracy.get_gdt_for_dir(args.native_pdb, args.input_pdb_dir)
    tmscore_df = tmscore_df.sort_values('Model').reset_index(drop=True)
    tmscore_df.to_csv(args.output_csv_path)

if __name__ == '__main__':
    main()
