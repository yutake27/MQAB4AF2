"""
Mv local scores to the specified directory.
"""

import argparse
import re
import shutil
from pathlib import Path


def mv(src, dist):
    shutil.move(src, dist)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dst_dir', type=str, help='destination of local score.')
    args = parser.parse_args()

    dst_dir = Path(args.dst_dir)
    dst_dir.mkdir(exist_ok=True)
    src_dir = dst_dir.parent
    pattern = '.*.features.npz'
    repatter = re.compile(pattern)
    for npz_path in src_dir.glob('*.npz'):
        if repatter.match(str(npz_path)) is None:
            mv(npz_path, dst_dir / npz_path.name)
