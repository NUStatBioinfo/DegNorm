import pytest
from pandas import DataFrame
from degnorm.utils import *
from numpy import ndarray

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------- #
# subset_to_chrom tests
# ----------------------------------------------------- #
@pytest.fixture
def sample_gene_data():
    df = DataFrame({'chr': ['chr1', 'chr1', 'chr2', 'chr2', 'chr3', 'chr3']
                    , 'gene': ['AZMP1', 'GENEE', 'ATP5F1', 'ATP5F1', 'GENE20', 'GENE50']
                    , 'exon_start': [105, 170, 200, 300, 550, 670]})
    return df


def test_positive_subset_to_chrom(sample_gene_data):
    chr2_df = subset_to_chrom(sample_gene_data, chrom='chr2', reindex=False)
    chr1_3_df = subset_to_chrom(sample_gene_data, chrom=['chr1', 'chr3'], reindex=True)

    assert isinstance(chr2_df, DataFrame)
    assert isinstance(chr1_3_df, DataFrame)

    assert len(chr2_df.chr.unique()) == 1
    assert chr2_df.chr.unique() == 'chr2'
    assert chr2_df.shape == (2, 3)

    assert len(chr1_3_df.chr.unique()) == 2
    assert chr1_3_df.shape == (4, 3)


def test_negative_subset_to_chrom(sample_gene_data):
    with pytest.raises(ValueError):
        subset_to_chrom(sample_gene_data, chrom=['not_a_chromosome'])


# ----------------------------------------------------- #
# list manipulation tests
# ----------------------------------------------------- #
def test_flatten_2d():
    l = [['gene' + str(i) for i in range(x)] for x in range(1, 5)]
    flat = flatten_2d(l)

    assert isinstance(flat, ndarray)
    assert flat.shape == (10,)
    assert len(flat) == 10


def test_split_into_chunks():
    l = ['gene' + str(i) for i in range(10)]

    assert len(split_into_chunks(l, n=3)) == 3
    assert len(split_into_chunks(l, n=7)) == 5
    assert len(split_into_chunks(l, n=10)) == 10
    assert len(split_into_chunks(l, n=20)) == 10