import subprocess
from pathlib import Path
from typing import List, Tuple, Union
import pandas as pd

from tqdm import tqdm


class ModelAccuracy:
    @staticmethod
    def _parse_TMscore(result) -> Tuple[float, float, float]:
        lines = result.split('\n')
        for line in lines:
            line_split = line.split()
            if len(line_split) == 0:
                continue
            elif line_split[0] == 'TM-score':
                tmscore = float(line_split[2])
            elif line_split[0] == 'GDT-TS-score=':
                gdtts = float(line_split[1])
            elif line_split[0] == 'GDT-HA-score=':
                gdtha = float(line_split[1])
        return tmscore, gdtts, gdtha

    @staticmethod
    def _run_TMscore(model_pdb, native_pdb) -> str:
        cmd = ['TMscore', model_pdb, native_pdb, '-outfmt', '-1']
        result = subprocess.check_output(cmd)
        return result.decode('utf-8')

    @classmethod
    def get_gdt(cls, model_pdb, native_pdb) -> Tuple[float, float, float]:
        result = cls._run_TMscore(model_pdb, native_pdb)
        tmscore, gdtts, gdtha = cls._parse_TMscore(result)
        return tmscore, gdtts, gdtha

    @classmethod
    def get_gdt_for_dir(cls, native_pdb_path: str, model_pdb_dir: str) -> pd.DataFrame:
        results = []
        for model in tqdm(list(Path(model_pdb_dir).glob('*.pdb'))):
            tmscore, gdtts, gdtha = cls.get_gdt(model, native_pdb_path)
            results.append([model.stem, gdtts, gdtha, tmscore])

        df = pd.DataFrame(results, columns=['Model', 'GDT_TS', 'GDT_HA', 'TMscore'])
        df = df.sort_index()
        return df