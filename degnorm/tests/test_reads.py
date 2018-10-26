import pytest
import os
from pandas import DataFrame
from random import choice
from degnorm.reads import ReadsProcessor

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------- #
# ReadsCoverageProcessor tests
# ----------------------------------------------------- #
@pytest.fixture
def rp_setup(request):
    salt = ''.join([choice(['a', 'b', 'c', 'd', 'e', 'f', 'g', '1', '2', '3', '4', '5', '6', '7']) for _ in range(20)])
    outdir = os.path.join(THIS_DIR, 'data', 'reads_test_' + salt)
    os.makedirs(outdir)

    def teardown():
        os.rmdir(outdir)

    request.addfinalizer(teardown)

    sam_file = os.path.join(THIS_DIR, 'data', 'hg_small_1.sam')  # sample paired reads data
    rp = ReadsProcessor(sam_file, n_jobs=1, output_dir=outdir)
    return rp


def test_rp_load(rp_setup):
    print('TESTING ReadsProcessor.load()')
    rp_setup.load()
    assert isinstance(rp_setup.header, DataFrame)
    assert isinstance(rp_setup.data, DataFrame)
    assert not rp_setup.header.empty
    assert not rp_setup.data.empty
    assert rp_setup.paired


def test_cigar_parser(rp_setup):
    print('TESTING ReadsProcessor._cigar_segment_bounds')

    # one match, all 100 base pairs covering positions 0 through 99 (inclusive)
    cigar_1 = '100M'

    # two matches, 0 -> 12 (inclusive), skip 20, 33 -> 132 (inclusive)
    cigar_2 = '13M10X10D100M'

    cigar_1_parse = rp_setup._cigar_segment_bounds(cigar_1, 0)
    assert cigar_1_parse == [0, 99]

    cigar_2_parse = rp_setup._cigar_segment_bounds(cigar_2, 0)
    assert cigar_2_parse == [0, 12, 33, 132]