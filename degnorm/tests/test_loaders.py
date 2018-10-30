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
def sam_loader_paired():
    sam_file = os.path.join(THIS_DIR, 'data', 'hg_small_1.sam')  # paired reads subsample dataset.
    sam_loader = SamLoader(sam_file)
    return sam_loader.get_data()


@pytest.fixture
def sam_loader_unpaired():
    sam_file = os.path.join(THIS_DIR, 'data', 'ff_small.sam')  # unpaired reads subsample dataset.
    sam_loader = SamLoader(sam_file)
    return sam_loader.get_data(unique_alignment=True)  # run with NH:i:1 only.


# check that paired and unpaired reads files are classified as such.
def test_sam_pairedness(sam_loader_paired, sam_loader_unpaired):
    assert sam_loader_paired['paired']
    assert not sam_loader_unpaired['paired']


# check that headers from .sam files are parsed properly for chromosome name and length metadata.
def test_sam_header(sam_loader_paired, sam_loader_unpaired):
    for header in [sam_loader_paired['header'], sam_loader_unpaired['header']]:
        assert isinstance(header, DataFrame)
        assert not header.empty
        assert all([col in header.columns.tolist() for col in ['chr', 'length']])

# check that the reads contents of .sam files are parsed properly dep on paired/unpaired status.
def test_sam_data(sam_loader_paired, sam_loader_unpaired):
    dat_paired = sam_loader_paired['data']
    dat_unpaired = sam_loader_unpaired['data']
    assert all([isinstance(df, DataFrame) for df in [dat_paired, dat_unpaired]])
    assert all([not df.empty for df in [dat_paired, dat_unpaired]])
    assert all([col in dat_paired.columns.tolist() for col in ['qname', 'chr', 'pos', 'cigar', 'rnext', 'qname_unpaired']])
    assert all([col in dat_unpaired.columns.tolist() for col in ['qname', 'chr', 'pos', 'cigar', 'rnext']])

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