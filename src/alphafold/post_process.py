import argparse
import tarfile
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from util.download_native_pdb import DownloadProtein
from util.get_model_accuracy import ModelAccuracy


data_dir = Path('../../data')


def main():
    parser = argparse.ArgumentParser(description='Post process the output of alphafold')
    parser.add_argument('alphafold_output_dir', type=str, help='alphafold output directory of a target.')
    args = parser.parse_args()

    alphafold_output_dir = Path(args.alphafold_output_dir)
    assert alphafold_output_dir.exists()
    target_name = alphafold_output_dir.stem
    dataset_name = alphafold_output_dir.parents[1].stem
    dataset_path = data_dir / 'out' / dataset_name
    assert dataset_path.exists()
    native_pdb_path = dataset_path / 'native_pdb' / (target_name + '.pdb')
    fasta_path = dataset_path / 'fasta' / (target_name + '.fasta')
    alphafold_score_path = alphafold_output_dir / 'scores.csv'
    output_label_path = dataset_path / 'score' / 'label' / (target_name + '.csv')

    # Download native pdb
    print('Downloading native pdb...')
    pdb_id, chain = target_name.split('_')
    print(f'PDB ID: {pdb_id}, Chain: {chain}')
    DownloadProtein.download_pdb_specific_chain(pdb_id, chain, str(fasta_path), str(native_pdb_path))
    assert native_pdb_path.exists()

    # Get model accuracy
    print('Get model accuracy...')
    label_df = ModelAccuracy.get_gdt_for_dir(native_pdb_path, alphafold_output_dir)
    score_df = pd.read_csv(alphafold_score_path, index_col=0)
    df = pd.merge(label_df, score_df, on='Model').sort_values('GDT_TS')
    df.to_csv(output_label_path)

    # Compress model_*.pickle
    tar_path = alphafold_output_dir / 'model_pickle.tar.gz'
    if not tar_path.exists():
        print('Compressing model.pickle...')
        with tarfile.open(str(tar_path), 'w:gz') as tar:
            for pkl_path in tqdm(list(alphafold_output_dir.glob('model_*.pickle'))):
                tar.add(str(pkl_path), arcname=pkl_path.name)
        # delete model_*.pickle
        for pkl_path in alphafold_output_dir.glob('model_*.pickle'):
            pkl_path.unlink()


if __name__ == "__main__":
    main()
