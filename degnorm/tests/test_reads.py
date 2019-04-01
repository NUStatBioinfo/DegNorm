import pytest
import os
import shutil
from pandas import DataFrame, read_csv
from random import choice
from pysam.libcalignmentfile import AlignmentFile
from degnorm.reads import *
from degnorm.gene_processing import GeneAnnotationProcessor, get_gene_overlap_structure

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------- #
# define fixtures for BamReadsCoverageProcessor
# ----------------------------------------------------- #
@pytest.fixture
def bam_setup(request):
    salt = ''.join([choice(['a', 'b', 'c', 'd', 'e', 'f', 'g', '1', '2', '3', '4', '5', '6', '7']) for _ in range(20)])
    outdir = os.path.join(THIS_DIR, 'data', 'reads_test_' + salt)
    os.makedirs(outdir)

    def teardown():
        shutil.rmtree(outdir)

    request.addfinalizer(teardown)

    bam_file_paired = os.path.join(THIS_DIR, 'data', 'hg_small_1.bam')  # sample paired reads data
    bai_file_paired = os.path.join(THIS_DIR, 'data', 'hg_small_1.bai')  # corresponding .bam index file
    bam_file_single = os.path.join(THIS_DIR, 'data', 'ff_small.bam')  # sample single-end reads data
    bai_file_single = os.path.join(THIS_DIR, 'data', 'ff_small.bai')  # corresponding .bai index file
    bam_processor_paired = BamReadsProcessor(bam_file_paired
                                             , index_file=bai_file_paired
                                             , n_jobs=1
                                             , output_dir=outdir)
    bam_processor_single = BamReadsProcessor(bam_file_single
                                             , index_file=bai_file_single
                                             , n_jobs=1
                                             , output_dir=outdir)
    return bam_processor_paired, bam_processor_single


# need .gtf file in order to test coverage parsing.
# GeneAnnotationLoader should have already been tested.
@pytest.fixture
def gtf_setup():
    gtf_file = os.path.join(THIS_DIR, 'data', 'chr1_small.gtf')
    gtf_processor = GeneAnnotationProcessor(gtf_file)
    exon_df = gtf_processor.run()
    return exon_df

# ----------------------------------------------------- #
# BamReadsCoverageProcessor tests
# ----------------------------------------------------- #

# test that .bam file header gets read correctly
def test_bam_header_paired(bam_setup):
    bam_setup = bam_setup[0]
    bamfile = bam_setup.loader.get_data()
    assert isinstance(bam_setup.header, DataFrame)
    assert isinstance(bamfile, AlignmentFile)
    assert not bam_setup.header.empty
    assert bam_setup.paired
    bamfile.close()


def test_bam_header_unpaired(bam_setup):
    bam_setup = bam_setup[1]
    bamfile = bam_setup.loader.get_data()
    assert isinstance(bam_setup.header, DataFrame)
    assert isinstance(bamfile, AlignmentFile)
    assert not bam_setup.header.empty
    assert not bam_setup.paired
    bamfile.close()


# test that paired read .bam files are loaded correctly.
def test_bam_load_paired(bam_setup):
    reqd_cols = ['qname', 'pos', 'cigar', 'qname_unpaired']
    bam_setup = bam_setup[0]
    reads_df = bam_setup.load_chromosome_reads('chr1')
    assert isinstance(reads_df, DataFrame)
    assert not reads_df.empty
    assert reads_df.shape[1] == 4
    assert all([col in reads_df.columns.tolist() for col in reqd_cols])


# test that single-end .bam files are loaded correctly.
def test_bam_load_single(bam_setup):
    reqd_cols = ['qname', 'pos', 'cigar']
    bam_setup = bam_setup[1]
    read_count_df = bam_setup.load_chromosome_reads('chr1')
    assert isinstance(read_count_df, DataFrame)
    assert not read_count_df.empty
    assert read_count_df.shape[1] == 3
    assert all([col in read_count_df.columns.tolist() for col in reqd_cols])


# test coverage / read count calculations on paired alignment file.
def test_bam_coverage_counts_paired(bam_setup, gtf_setup):
    bam_setup = bam_setup[0]
    exon_df = gtf_setup
    gene_df = exon_df[['chr', 'gene', 'gene_start', 'gene_end']].drop_duplicates().reset_index(drop=True)
    gene_overlap_dat = {'chr1': get_gene_overlap_structure(gene_df)}

    out = bam_setup.coverage_read_counts(gene_overlap_dat
                                         , gene_df=gene_df
                                         , exon_df=exon_df)

    output_files = os.listdir(bam_setup.save_dir)

    # check that chromosome coverage file and read counts file exist.
    assert 'chrom_coverage_hg_small_1_chr1.npz' in output_files
    assert 'read_counts_hg_small_1_chr1.csv' in output_files

    # check read counts file.
    reads_df = read_csv(os.path.join(bam_setup.save_dir, 'read_counts_hg_small_1_chr1.csv'))
    assert not reads_df.empty
    assert len(list(set(reads_df.columns.tolist()) - {'gene', 'hg_small_1'})) == 0


# test coverage / read count calculations on single-end reads alignment file.
def test_bam_coverage_counts_single(bam_setup, gtf_setup):
    bam_setup = bam_setup[1]
    exon_df = gtf_setup
    gene_df = exon_df[['chr', 'gene', 'gene_start', 'gene_end']].drop_duplicates().reset_index(drop=True)
    gene_overlap_dat = {'chr1': get_gene_overlap_structure(gene_df)}

    out = bam_setup.coverage_read_counts(gene_overlap_dat
                                         , gene_df=gene_df
                                         , exon_df=exon_df)

    output_files = os.listdir(bam_setup.save_dir)
    print('OUTPUT FILES:')
    print(output_files)

    # check that chromosome coverage file and read counts file exist.
    assert 'chrom_coverage_ff_small_chr1.npz' in output_files
    assert 'read_counts_ff_small_chr1.csv' in output_files

    # check read counts file.
    reads_df = read_csv(os.path.join(bam_setup.save_dir, 'read_counts_ff_small_chr1.csv'))
    assert not reads_df.empty
    assert len(list(set(reads_df.columns.tolist()) - {'gene', 'ff_small'})) == 0


# ----------------------------------------------------- #
# Other degnorm.reads module tests
# ----------------------------------------------------- #
def test_cigar_parser():

    # one match, all 100 base pairs covering positions 0 through 99 (inclusive)
    cigar_1 = '100M'

    # two matches, 0 -> 12 (inclusive), skip 20, 33 -> 132 (inclusive)
    cigar_2 = '13M10X10D100M'

    cigar_1_parse = cigar_segment_bounds(cigar_1, 0)
    assert cigar_1_parse == [0, 99]

    cigar_2_parse = cigar_segment_bounds(cigar_2, 0)
    assert cigar_2_parse == [0, 12, 33, 132]


def test_fill_in_bounds():

    bounds_pass = np.array([10, 15, 40, 50, 60, 65])
    expected_vec = np.array([10, 11, 12, 13, 14, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 60, 61, 62, 63, 64])
    fill_vec = fill_in_bounds(bounds_pass)
    assert np.array_equal(fill_vec, expected_vec)

    bounds_fail = bounds_pass[0:5]
    ValueError('bounds_vec must have even number of values.')
    with pytest.raises(ValueError, message='bounds_vec must have even number of values.'):
        fill_in_bounds(bounds_fail)

