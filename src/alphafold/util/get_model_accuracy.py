import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Union

import joblib
import numpy as np
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
    def get_gdt_for_dir(cls, native_pdb_path: Union[str, Path], model_pdb_dir: Union[str, Path]) -> pd.DataFrame:
        def get_gdt(model_pdb: Path):
            return model_pdb.stem, *cls.get_gdt(model_pdb, native_pdb_path)

        results = joblib.Parallel(n_jobs=-1)(
            joblib.delayed(get_gdt)(model_pdb) for model_pdb in tqdm(list(Path(model_pdb_dir).glob('*.pdb')))
        )

        df = pd.DataFrame(results, columns=['Model', 'GDT_TS', 'GDT_HA', 'TMscore'])
        df = df.sort_index()
        return df

    @staticmethod
    def _parse_lddt(result: str) -> Tuple[float, float, np.ndarray]:
        lines = result.split('\n')
        global_lddt = float(lines[7].split()[-1])
        local_lddts = [float(lddt) if (lddt := line.split()[4]) != '-' else np.nan for line in lines[11:] if len(line)]
        local_lddts_except_none = np.array(list(filter(lambda x: x, local_lddts)))
        mean_lddt = float(np.mean(local_lddts_except_none))
        return global_lddt, mean_lddt, np.array(local_lddts)

    @staticmethod
    def _run_lddt(model_pdb, native_pdb) -> str:
        cmd = ['lddt', model_pdb, native_pdb]
        result = subprocess.check_output(cmd)
        return result.decode('utf-8')

    @classmethod
    def get_lddt(cls, model_pdb, native_pdb) -> Tuple[float, float, np.ndarray]:
        result = cls._run_lddt(model_pdb, native_pdb)
        global_lddt, mean_lddt, local_lddts = cls._parse_lddt(result)
        return global_lddt, mean_lddt, local_lddts

    @classmethod
    def get_lddt_for_dir(cls, native_pdb_path: Union[str, Path],
            model_pdb_dir: Union[str, Path]) -> Tuple[pd.DataFrame, Dict[str, np.ndarray]]:
        global_results = []
        local_lddt_dict = {}

        def get_lddt(model_pdb: Path):
            return model_pdb.stem, *cls.get_lddt(model_pdb, native_pdb_path)

        results = joblib.Parallel(n_jobs=-1)(
            joblib.delayed(get_lddt)(model_pdb) for model_pdb in tqdm(list(Path(model_pdb_dir).glob('*.pdb')))
        )
        if results is None:
            raise ValueError(f'No models found in {model_pdb_dir}')
        for result in results:
            model_name, global_lddt, mean_lddt, local_lddts = result
            global_results.append([model_name, global_lddt, mean_lddt])
            local_lddt_dict[model_name] = local_lddts

        global_df = pd.DataFrame(global_results, columns=['Model', 'Global_LDDT', 'Mean_LDDT'])
        global_df = global_df.sort_index()
        res_indices = np.arange(1, len(local_lddts) + 1)
        local_lddt_dict['res_indices'] = res_indices
        return global_df, local_lddt_dict
