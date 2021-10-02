import pytest
import pandas as pd
from download_native_pdb import DownloadProtein


# Sample case
#     pdb_id = '6VI2'
#     chain = 'C'
#     fasta_seq = 'SDIQMTQSPSSLSASVGDRVTITCRASQSVSSAVAWYQQKPGKAPKLLIYSASSLYSGVPSRFS...'
#     output_pdb_path = tmpdir / pdb_id + '_' + chain + '.pdb'
#     DownloadProtein.download_pdb_specific_chain(pdb_id, chain, fasta_seq, output_pdb_path)


df = pd.read_csv('../../../data/interim/target_list.csv', index_col=0)
num_test = 10


@pytest.fixture(params=range(num_test))
def load_test_case(request, tmpdir):
    index = request.param
    row = df.iloc[index]
    entry_id = row['id']
    pdb_id, chain = entry_id.split('_')
    sequence = row['sequence']
    tmpfile = tmpdir.join(pdb_id + '_' + chain + '.pdb')
    yield (pdb_id, chain, sequence, tmpfile)
    tmpfile.remove()


def test_download_protein(load_test_case):
    DownloadProtein.download_pdb_specific_chain(*load_test_case)
