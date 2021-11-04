import subprocess
import pandas as pd

from tqdm import tqdm


class ModelAccuracy:
    @staticmethod
    def _parse_TMscore(result):
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
    def _run_TMscore(native_pdb, model_pdb):
        cmd = ['TMscore', native_pdb, model_pdb, '-outfmt', '-1']
        result = subprocess.check_output(cmd)
        return result.decode('utf-8')

    @classmethod
    def get_gdt(cls, native_pdb, model_pdb):
        result = cls._run_TMscore(native_pdb, model_pdb)
        tmscore, gdtts, gdtha = cls._parse_TMscore(result)
        return tmscore, gdtts, gdtha

    @classmethod
    def get_gdt_for_dir(cls, native_pdb_path, model_pdb_dir) -> pd.DataFrame:
        results = []
        for model in tqdm(list(model_pdb_dir.glob('*.pdb'))):
            tmscore, gdtts, gdtha = cls.get_gdt(native_pdb_path, model)
            results.append([model.stem, gdtts, gdtha, tmscore])

        df = pd.DataFrame(results, columns=['Model', 'GDT_TS', 'GDT_HA', 'TMscore'])
        df = df.sort_index()
        return df