"""
Compress pdb files and scores.

Source
.
├── fasta
│   ├── target1.fasta
│   └── target2.fasta
├── native_pdb
│   ├── target1.pdb
│   └── target2.pdb
├── alphafold_output
    ├── target1
    └── target2

Output
.
├── fasta
│   ├── target1.fasta
│   └── target2.fasta
├── native_pdb
│   ├── target1.pdb
│   └── target2.pdb
├── pdb
│   ├── target1
│   └── target2
└── score

"""


import tarfile
from cProfile import label
from pathlib import Path

import pandas as pd
import tqdm
from chainer import dataset

data_dir = Path('../../data/out/dataset')

def compress_pdb(target_df: pd.DataFrame, output_dir: Path) -> Path:
    """Compress pdb files

    Args:
        target_df (pd.DataFrame): DataFrame of targets
        output_dir (Path): Output directory

    Returns:
        Path: output tar file path
    """
    print('Compress pdb')
    output_tar_pdb_path = output_dir / 'pdb.tar.gz'
    if output_tar_pdb_path.exists():
        return output_tar_pdb_path
    pdb_dir = data_dir / 'alphafold_output'
    with tarfile.open(output_tar_pdb_path, 'w:gz') as f:
        for target in tqdm.tqdm(target_df['id'].unique()):
            for pdb in pdb_dir.glob(f'{target}/*.pdb'):
                f.add(str(pdb), f'pdb/{target}/{pdb.name}')
    return output_tar_pdb_path


def compress_native_pdb(tf: tarfile.TarFile, target_df: pd.DataFrame, arcdir: str) -> None:
    # Compress native pdb files
    print('Compress native pdb')
    native_pdb_dir = data_dir / 'native_pdb'
    for target in target_df['id'].unique():
        native_pdb_path = native_pdb_dir / f'{target}.pdb'
        tf.add(str(native_pdb_path), arcname=f'{arcdir}/native_pdb/{target}.pdb')


def compress_fasta(tf: tarfile.TarFile, target_df: pd.DataFrame, arcdir: str) -> None:
    # Compress fasta files
    print('Compress fasta')
    fasta_dir = data_dir / 'fasta'
    for target in target_df['id'].unique():
        fasta_path = fasta_dir / f'{target}.fasta'
        tf.add(str(fasta_path), arcname=f'{arcdir}/fasta/{target}.fasta')


def compress_score(tf: tarfile.TarFile, target_df: pd.DataFrame, arcdir: str,
                   dataset_name: str, output_dir: Path) -> None:
    print('Compress score')
    # Compress score files
    score_dir = data_dir / 'score/subsets' / dataset_name
    # add target list
    domain_df = pd.read_csv(Path('../../data/interim') / 'ecod_domain_num.csv', index_col=0)
    neff_df = pd.read_csv(score_dir / 'neff.csv', index_col=0)
    target_df_ = pd.merge(target_df, domain_df, on='id', how='left')
    target_df_ = target_df_.rename(columns={'id': 'Target'})
    target_df_ = pd.merge(target_df_, neff_df, on='Target', how='left')
    target_df_['num_similar_seq_in_AF2_training'] =\
        target_df_['num_entry_in_cluster'] - target_df_['num_entry_in_cluster_AF2_notInclude']
    target_df_ = target_df_.drop(['num_entry_in_cluster_AF2_notInclude'], axis=1)
    target_list_path = output_dir / 'target_list.csv'
    target_df_.to_csv(target_list_path)
    tf.add(str(target_list_path), arcname=f'{arcdir}/score/target_list.csv')

    # read label
    label_df = pd.read_csv(score_dir / 'label.csv', index_col=0)
    label_df = label_df.drop(['Num_diff', 'Num_missing', 'Length', 'pLDDT', 'pTMscore'], axis=1)
    # read MQA methods
    methods = ['af2_confidence', 'DOPE', 'ProQ3D', 'SBROD', 'VoroCNN', 'P3CMQA', 'DeepAccNet']
    for i, method in enumerate(methods):
        method_path = score_dir / f'{method}_resolved.csv.gz'
        method_df = pd.read_csv(method_path, index_col=0)
        if i == 0:
            mqa_df = method_df
            continue
        mqa_df = pd.merge(mqa_df, method_df, on=['Model', 'Target'], how='inner')
    mqa_df = mqa_df.rename({m: m.split('_')[0] for m in mqa_df.columns.tolist() if '_resolved' in m}, axis=1)
    mqa_df = mqa_df.drop(['ProQ2D', 'ProQRosCenD', 'ProQRosFAD'], axis=1)
    all_score_df = pd.merge(label_df, mqa_df, on=['Target', 'Model'], how='inner')
    score_csv_path = output_dir / 'all_score.csv'
    all_score_df.to_csv(score_csv_path)
    tf.add(str(score_csv_path), arcname=f'{arcdir}/score/all_score.csv')


def main():
    dataset_name = 'target_subset_how_eq_random_num_500_seed_0'
    csv_path = Path('../../data/interim') / f'{dataset_name}.csv'
    target_df = pd.read_csv(csv_path, index_col=0)
    # output directory
    output_dir = Path('../../data/out/dataset/compress') / dataset_name
    output_dir.mkdir(parents=True, exist_ok=True)

    # compress pdb
    output_tar_pdb_path = compress_pdb(target_df, output_dir)

    # compress other files
    output_tar_filename = 'af2_500_targets.tar.gz'
    output_tar_path = output_dir / output_tar_filename
    with tarfile.open(output_tar_path, 'w:gz') as tf:
        arcdir = output_tar_filename
        # add pdb tar file
        tf.add(output_tar_pdb_path, arcname=f'{arcdir}/pdb.tar.gz')
        # add native_pdb
        compress_native_pdb(tf, target_df, arcdir)
        # add fasta
        compress_fasta(tf, target_df, arcdir)
        # add score
        compress_score(tf, target_df, arcdir, dataset_name, output_dir)
        # add README
        readme_path = data_dir / 'compress' / 'README.md'
        tf.add(str(readme_path), arcname=f'{arcdir}/README.md')


if __name__ == '__main__':
    main()
