import pytest
import os
from pandas import DataFrame
from degnorm.gene_processing import GeneAnnotationProcessor

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------- #
# GeneAnnotationProcessor tests
# ----------------------------------------------------- #

@pytest.fixture
def gene_processor_setup():
    gtf_file = os.path.join(THIS_DIR, 'data', 'chr1_small.gtf')
    gtf_processor = GeneAnnotationProcessor(gtf_file
                                            , n_jobs=1)
    return gtf_processor


def test_gtf_processor_load(gene_processor_setup):
    exons_df = gene_processor_setup.load()
    assert isinstance(exons_df, DataFrame)
    assert not exons_df.empty


def test_gtf_processor_run(gene_processor_setup):
    reqd_cols = ['chr', 'gene', 'gene_start', 'gene_end', 'start', 'end']
    exons_df = gene_processor_setup.run()
    assert isinstance(exons_df, DataFrame)
    assert not exons_df.empty
    assert all([col in exons_df.columns.tolist() for col in reqd_cols])