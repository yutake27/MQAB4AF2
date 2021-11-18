import pytest
import pandas as pd
from download_native_pdb import DownloadProtein


# Sample case
#     pdb_id = '6VI2'
#     chain = 'C'
#     fasta_seq = 'SDIQMTQSPSSLSASVGDRVTITCRASQSVSSAVAWYQQKPGKAPKLLIYSASSLYSGVPSRFS...'
#     output_pdb_path = tmpdir / pdb_id + '_' + chain + '.pdb'
#     DownloadProtein.download_pdb_specific_chain(pdb_id, chain, fasta_seq, output_pdb_path)


def make_fasta(header, sequence, fasta_file):
    with open(fasta_file, 'w') as f:
        f.write(f'{header}\n{sequence}\n')


df = pd.read_csv('../../../data/interim/target_list.csv', index_col=0)
num_test = 5


@pytest.fixture(params=range(num_test))
def load_test_case(request, tmpdir):
    index = request.param
    row = df.iloc[index]
    entry_id = row['id']
    pdb_id, chain = entry_id.split('_')
    header = row['header']
    sequence = row['sequence']
    tmp_fastafile = tmpdir.join(f'{entry_id}.fasta')
    make_fasta(header, sequence, tmp_fastafile)
    tmp_pdbfile = tmpdir.join(f'{entry_id}.pdb')
    yield (pdb_id, chain, tmp_fastafile, tmp_pdbfile)
    tmp_pdbfile.remove()


def test_download_protein(load_test_case):
    DownloadProtein.download_pdb_specific_chain(*load_test_case)


missing_caces = ['7OF7_8', '5ZX9_A', '5QU5_A', '6HCZ_A', '7DCO_W']


@pytest.fixture(params=missing_caces)
def load_test_missing_case(request, tmpdir):
    entry_id = request.param
    pdb_id, chain = entry_id.split('_')
    target_row = df.query('id == @entry_id').iloc[0]
    header = target_row['header']
    sequence = target_row['sequence']
    print(header, sequence)
    tmp_fastafile = tmpdir.join(f'{entry_id}.fasta')
    make_fasta(header, sequence, tmp_fastafile)
    tmp_pdbfile = tmpdir.join(f'{entry_id}.pdb')
    yield (pdb_id, chain, tmp_fastafile, tmp_pdbfile)
    tmp_pdbfile.remove()
    tmp_fastafile.remove()


def test_download_protein_missing(load_test_missing_case):
    DownloadProtein.download_pdb_specific_chain(*load_test_missing_case)
    _, _, tmp_fastafile, tmp_pdbfile = load_test_missing_case
    num_diff, num_missing = DownloadProtein.test_match_pdb_and_fasta(tmp_pdbfile, tmp_fastafile)
    assert num_diff < 10
