import pytest
import os
from pandas import DataFrame
from pysam.libcalignmentfile import AlignmentFile
from degnorm.loaders import Loader, BamLoader, GeneAnnotationLoader

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------- #
# Loader tests
# ----------------------------------------------------- #
def test_loader_file_dne_error():
    nonexistent_file = 'id_never_name_a_file_this.txt'
    l = Loader('.txt')
    with pytest.raises(FileNotFoundError, message='to_load file {0} not found'.format(nonexistent_file)):
        l.get_file(nonexistent_file)


# ----------------------------------------------------- #
# SamLoader tests
# ----------------------------------------------------- #
def test_bam_loader():
    bam_file = os.path.join(THIS_DIR, 'data', 'hg_small_1.bam')  # paired reads subsample dataset.
    bai_file = os.path.join(THIS_DIR, 'data', 'hg_small_1.bai')  # corresponding .bam index file
    bam_loader = BamLoader(bam_file, index_file=bai_file)
    bam_file = bam_loader.get_data()
    assert isinstance(bam_file, AlignmentFile)
    bam_file.close()


def test_bam_loader_index_file_dne_error():
    bam_file = os.path.join(THIS_DIR, 'data', 'ff_small.bam')
    nonexistent_bai_file = 'id_also_never_name_a_file_this.txt'
    with pytest.raises(FileNotFoundError, message='{0} .bam index file not found'.format(nonexistent_bai_file)):
        bam_loader = BamLoader(bam_file, index_file=nonexistent_bai_file)


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
    reqd_cols = ['chr', 'start', 'end', 'gene']
    assert isinstance(gene_loader, DataFrame)
    assert not gene_loader.empty
    assert all([col in gene_loader.columns.tolist() for col in reqd_cols])