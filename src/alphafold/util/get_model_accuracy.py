import subprocess
import pandas as pd


class ModelAccuracy:
    @staticmethod
    def _parse_TMscore(result):
        lines = result.split('\n')
        for line in lines:
            line_split = line.split()
            if len(line_split) == 0:
                continue
            elif line_split[0] == 'GDT-TS-score=':
                gdtts = line_split[1]
            elif line_split[0] == 'GDT-HA-score=':
                gdtha = line_split[1]
        return gdtts, gdtha

    @staticmethod
    def _run_TMscore(native_pdb, model_pdb):
        cmd = ['TMscore', native_pdb, model_pdb, '-outfmt', '-1']
        result = subprocess.check_output(cmd)
        return result.decode('utf-8')

    @classmethod
    def get_gdt(cls, native_pdb, model_pdb):
        result = cls._run_TMscore(native_pdb, model_pdb)
        gdtts, gdtha = cls._parse_TMscore(result)
        return gdtts, gdtha

    @classmethod
    def get_gdt_for_target_df(cls, native_pdb_path, model_pdb_dir) -> pd.DataFrame:
        model_array = []
        gdtts_array = []
        gdtha_array = []

        for model in model_pdb_dir.iterdir():
            model_array.append(model.stem)
            gdtts, gdtha = cls.get_gdt(native_pdb_path, model)
            gdtts_array.append(gdtts)
            gdtha_array.append(gdtha)

        df = pd.DataFrame({'GDT_TS': gdtts_array, 'GDT_HA': gdtha_array}, index=model_array)
        df = df.astype('float')
        df = df.sort_index()
        return df
