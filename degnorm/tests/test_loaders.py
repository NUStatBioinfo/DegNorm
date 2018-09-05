import pytest
import os
from pandas import DataFrame
from degnorm.loaders import Loader, SamLoader, GeneAnnotationLoader

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------- #
# Loader tests
# ----------------------------------------------------- #
def test_loader_file_dne_error():
    nonexistent_file = 'id_never_name_a_file_this.txt'
    l = Loader('.txt')
    with pytest.raises(IOError, message='to_load file {0} not found'.format(nonexistent_file)):
        l.get_file(nonexistent_file)


# ----------------------------------------------------- #
# SamLoader tests
# ----------------------------------------------------- #
@pytest.fixture
def sam_loader():
    sam_file = os.path.join(THIS_DIR, 'data', 'hg_small_1.sam')
    sam_loader = SamLoader(sam_file)
    return sam_loader.get_data()


def test_sam_header(sam_loader):
    header = sam_loader['header']
    assert isinstance(header, DataFrame)
    assert not header.empty
    assert all([col in header.columns.tolist() for col in ['chr', 'length']])


def test_sam_data(sam_loader):
    dat = sam_loader['data']
    assert isinstance(dat, DataFrame)
    assert not dat.empty
    assert all([col in dat.columns.tolist() for col in ['qname', 'chr', 'pos', 'cigar', 'rnext', 'qname_unpaired']])

# ----------------------------------------------------- #
# GeneAnnotationLoader tests
# ----------------------------------------------------- #
@pytest.fixture
def gene_loader():
    print('TESTING GeneAnnotationLoader')
    gtf_file = os.path.join(THIS_DIR, 'data', 'chr1_small.gtf')
    gtf_loader = GeneAnnotationLoader(gtf_file)
    return gtf_loader.get_data()


def test_gtf_data(gene_loader):
    assert isinstance(gene_loader, DataFrame)
    assert not gene_loader.empty
    assert all([col in gene_loader.columns.tolist() for col in ['chr', 'start', 'end', 'gene']])