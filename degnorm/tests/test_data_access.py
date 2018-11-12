import pytest
import os
import shutil
import pickle as pkl
import numpy as np
from matplotlib.figure import Figure
from pandas import DataFrame
from random import choice
from degnorm.data_access import get_coverage_data, get_coverage_plots

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------- #
# define fixtures for BamReadsCoverageProcessor
# ----------------------------------------------------- #
@pytest.fixture
def da_setup(request):
    salt = ''.join([choice(['a', 'b', 'c', 'd', 'e', 'f', 'g', '1', '2', '3', '4', '5', '6', '7']) for _ in range(20)])
    outdir = os.path.join(THIS_DIR, 'data_access_test_' + salt)
    os.makedirs(outdir)

    def teardown():
        shutil.rmtree(outdir)

    request.addfinalizer(teardown)

    # create dummy gene_exon_metadata.csv file.
    exon_df = DataFrame({'chr': 'chr1'
                         , 'start': [323892, 324288, 1152288, 1153068]
                         , 'end': [324060, 324345, 1153838, 1154013]
                         , 'gene': ['GENE_1', 'GENE_1', 'GENE_2', 'GENE_2']
                         , 'exon_id': np.arange(4)
                         , 'gene_start': [323892, 323892, 1152288, 1152288]
                         , 'gene_end': [665731, 665731, 1167447, 1167447]})

    exon_df = exon_df[['chr', 'start', 'end', 'gene', 'exon_id', 'gene_end', 'gene_start']]
    exon_df.to_csv(os.path.join(outdir, 'gene_exon_metadata.csv')
                   , index=False)

    # create dummy DI scores file for gene_1, gene_2.
    di_df = DataFrame({'chr': ['chr1', 'chr1']
                       , 'gene': ['GENE_1', 'GENE_2']
                       , 'sample_1': [0.55, 0.11]
                       , 'sample_2': [0.46, 0.04]})
    di_df.to_csv(os.path.join(outdir, 'degradation_index_scores.csv')
                 , index=False)

    # dummy read counts to satisfy appearance of a DegNorm output directory.
    di_df.to_csv(os.path.join(outdir, 'read_counts.csv')
                 , index=False)

    # create chromosome 1 output dir.
    os.mkdir(os.path.join(outdir, 'chr1'))

    # create fake coverage matrices.
    cov_mat_dict = dict()
    cov_mat_dict['GENE_1'] = np.random.negative_binomial(n=400, p=0.8, size=[2, 225]).astype(float)
    cov_mat_dict['GENE_2'] = np.random.negative_binomial(n=400, p=0.8, size=[2, 2495]).astype(float)

    # save fake coverage matrices.
    with open(os.path.join(outdir, 'chr1', 'coverage_matrices_chr1.pkl'), 'wb') as f:
        pkl.dump(cov_mat_dict, f)

    with open(os.path.join(outdir, 'chr1', 'estimated_coverage_matrices_chr1.pkl'), 'wb') as f:
        pkl.dump(cov_mat_dict, f)

    return outdir


# ----------------------------------------------------- #
# degnorm.data_access.get_coverage_data tests
# ----------------------------------------------------- #
def test_get_coverage_data(da_setup):
    genes = ['GENE_1', 'GENE_2']
    cov_dat = get_coverage_data(genes=genes
                                , degnorm_dir=da_setup
                                , save_dir=da_setup)

    # check that coverage matrices were loaded.
    assert all([isinstance(cov_dat[gene]['raw'], DataFrame) for gene in genes])
    assert all([isinstance(cov_dat[gene]['estimate'], DataFrame) for gene in genes])

    # check that coverage matrices were saved.
    assert all([os.path.isfile(os.path.join(da_setup, 'chr1', x)) for x in
                ['GENE_1_raw_coverage.txt', 'GENE_1_estimated_coverage.txt']])


# ----------------------------------------------------- #
# degnorm.data_access.get_coverage_plots tests
# ----------------------------------------------------- #
def test_get_coverage_plots_no_save(da_setup):
    genes = ['GENE_1', 'GENE_2']
    out = get_coverage_plots(genes
                             , degnorm_dir=da_setup)

    assert len(out) == 2
    assert isinstance(out[0], Figure)

def test_get_coverage_plots_with_save(da_setup):
    genes = ['GENE_1', 'GENE_2']
    out = get_coverage_plots(genes
                             , degnorm_dir=da_setup
                             , save_dir=da_setup)

    assert len(out) == 2
    assert all([os.path.isfile(os.path.join(da_setup, 'chr1', x)) for x in out])