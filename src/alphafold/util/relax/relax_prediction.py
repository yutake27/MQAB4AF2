"""
Relax predicted structures by amber.
Input: Directory of model_output.pickle
Output: relaxed_protein.pdb
"""

import argparse
import os
import pickle
import sys
from pathlib import Path

from tqdm import tqdm

sys.path.append('../../alphafold/colabfold-conda/lib/python3.7/site-packages')
alphafold_dir = Path('../../alphafold')
sys.path.append(str(alphafold_dir))
from alphafold.relax import relax, utils
os.chdir(str(alphafold_dir))


def relax_prediction(pickle_path: Path, output_dir: Path):
    """Relax predicted structures by amber."""
    with open(pickle_path, 'rb') as f:
        data = pickle.load(f)
        protein = data['unrelaxed_protein']
        pred_output_path = output_dir / f'{pickle_path.stem}_relaxed.pdb'
        if not pred_output_path.exists():
            amber_relaxer = relax.AmberRelaxation(
                max_iterations=0,
                tolerance=2.39,
                stiffness=10.0,
                exclude_residues=[],
                max_outer_iterations=20)
            relaxed_pdb_lines, _, _ = amber_relaxer.process(prot=protein)
            with open(pred_output_path, 'w') as f:
                f.write(relaxed_pdb_lines)


def main():
    parser = argparse.ArgumentParser('Relax predicted structures by amber and calculate accuracy.')
    parser.add_argument('-p', '--pickle_dir', type=str, required=True, help='Directory of model_output.pickle')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Output directory of relaxed_protein.pdb')
    args = parser.parse_args()

    pickle_dir = Path(args.pickle_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for pickle_path in tqdm(list(pickle_dir.glob('model_*.pickle'))):
        relax_prediction(pickle_path, output_dir)


if __name__ == '__main__':
    main()
