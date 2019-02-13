import pytest
import shutil
from random import choice
from degnorm.reads import *
from degnorm.gene_processing import GeneAnnotationProcessor, get_gene_overlap_structure
from degnorm.reads_coverage_merge import *

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------- #
# define fixtures
# ----------------------------------------------------- #
@pytest.fixture
def bam_setup(request):
    salt = ''.join([choice(['a', 'b', 'c', 'd', 'e', 'f', 'g', '1', '2', '3', '4', '5', '6', '7']) for _ in range(20)])
    outdir = os.path.join(THIS_DIR, 'data', 'reads_test_' + salt)
    os.makedirs(outdir)

    def teardown():
        shutil.rmtree(outdir)

    request.addfinalizer(teardown)

    bam_file_1 = os.path.join(THIS_DIR, 'data', 'hg_small_1.bam')  # sample paired reads file 1
    bai_file_1 = os.path.join(THIS_DIR, 'data', 'hg_small_1.bai')  # corresponding .bam index file
    bam_file_2 = os.path.join(THIS_DIR, 'data', 'hg_small_2.bam')  # sample paired reads file 2
    bai_file_2 = os.path.join(THIS_DIR, 'data', 'hg_small_2.bai')  # corresponding .bam index file
    bam_processor_1 = BamReadsProcessor(bam_file_1
                                        , index_file=bai_file_1
                                        , n_jobs=1
                                        , output_dir=outdir)
    bam_processor_2 = BamReadsProcessor(bam_file_2
                                        , index_file=bai_file_2
                                        , n_jobs=1
                                        , output_dir=outdir)
    return [outdir, bam_processor_1, bam_processor_2]


# need .gtf file in order to test coverage parsing.
# GeneAnnotationLoader should have already run.
@pytest.fixture
def gtf_setup():
    gtf_file = os.path.join(THIS_DIR, 'data', 'chr1_small.gtf')
    gtf_processor = GeneAnnotationProcessor(gtf_file)
    exon_df = gtf_processor.run()
    return exon_df


# ----------------------------------------------------- #
# reads_merge module tests
# BamReadsCoverageProcessor tests should have already run.
# ----------------------------------------------------- #

def test_bam_coverage_counts_paired(bam_setup, gtf_setup):
    exon_df = gtf_setup
    gene_df = exon_df[['chr', 'gene', 'gene_start', 'gene_end']].drop_duplicates().reset_index(drop=True)
    gene_overlap_dat = {'chr1': get_gene_overlap_structure(gene_df)}

    cov_files = OrderedDict()
    read_count_files = OrderedDict()
    sample_ids = list()

    for i in [1, 2]:
        sample_id = bam_setup[i].sample_id
        sample_ids.append(sample_id)
        cov_files[sample_id], read_count_files[sample_id] = bam_setup[i].coverage_read_counts(gene_overlap_dat
                                                                                              , gene_df=gene_df
                                                                                              , exon_df=exon_df)

    read_counts_df = merge_read_count_files(read_count_files
                                            , chroms=['chr1'])

    gene_cov_dict = merge_gene_coverage_files(cov_files
                                              , exon_df=exon_df
                                              , n_jobs=2
                                              , output_dir=bam_setup[0])

    assert isinstance(read_counts_df, DataFrame)
    assert isinstance(gene_cov_dict, OrderedDict)
    assert gene_cov_dict.get(list(gene_cov_dict.keys())[0]).ndim == 2
    assert all(read_counts_df.columns == ['chr', 'gene'] + sample_ids)
    assert not read_counts_df.empty
    assert all([gene_cov_dict[x].ndim == 2 for x in gene_cov_dict])



