import pytest
import os
from pandas import DataFrame
from degnorm.gene_processing import GeneAnnotationProcessor, get_gene_overlap_structure

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------- #
# GeneAnnotationProcessor tests
# ----------------------------------------------------- #

@pytest.fixture
def gene_processor_setup():
    gtf_file = os.path.join(THIS_DIR, 'data', 'chr1_small.gtf')
    gtf_processor = GeneAnnotationProcessor(gtf_file)
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


# ----------------------------------------------------- #
# Other degnorm.gene_processing function tests
# ----------------------------------------------------- #

def test_get_gene_overlap_structure():
    genes_df = DataFrame({'gene': ['A', 'B', 'C', 'D']
                          , 'gene_start': [100, 150, 215, 600]
                          , 'gene_end': [200, 230, 280, 822]
                          , 'chr': ['chr2']*4})
    reqd_keys = ['isolated_genes', 'overlap_genes']
    gene_overlap_dat = get_gene_overlap_structure(genes_df)
    assert all([gene_overlap_dat.get(key) is not None for key in reqd_keys])
    assert len(gene_overlap_dat['overlap_genes']) == 1
    assert len(set(gene_overlap_dat['overlap_genes'][0]) - {'A', 'B', 'C'}) == 0
    assert len(set(gene_overlap_dat['isolated_genes']) - {'D'}) == 0
