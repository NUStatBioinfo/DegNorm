import pytest
import os
from pandas import DataFrame
from random import choice
from degnorm.reads import ReadsCoverageProcessor

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------- #
# ReadsCoverageProcessor tests
# ----------------------------------------------------- #
@pytest.fixture
def rcp_setup(request):
    salt = ''.join([choice(['a', 'b', 'c', 'd', 'e', 'f', 'g', '1', '2', '3', '4', '5', '6', '7']) for _ in range(20)])
    outdir = os.path.join(THIS_DIR, 'data', 'reads_test_' + salt)
    os.makedirs(outdir)

    def teardown():
        os.rmdir(outdir)

    request.addfinalizer(teardown)

    sam_file = os.path.join(THIS_DIR, 'data', 'hg_small_1.sam')
    rcp = ReadsCoverageProcessor(sam_file, n_jobs=1, tmp_dir=outdir)
    return rcp


def test_rcp_load(rcp_setup):
    rcp_setup.load()
    assert isinstance(rcp_setup.header, DataFrame)
    assert isinstance(rcp_setup.data, DataFrame)
    assert not rcp_setup.header.empty
    assert not rcp_setup.data.empty


def test_cigar_parser(rcp_setup):

    # one match, all 100 base pairs covering positions 0 through 99 (inclusive)
    cigar_1 = '100M'

    # two matches, 0 -> 12 (inclusive), skip 20, 33 -> 132 (inclusive)
    cigar_2 = '13M10X10D100M'

    cigar_1_parse = rcp_setup._cigar_segment_bounds(cigar_1, 0)
    assert cigar_1_parse[1] == 1
    assert cigar_1_parse[0] == [0, 99]

    cigar_2_parse = rcp_setup._cigar_segment_bounds(cigar_2, 0)
    assert cigar_2_parse[1] == 2
    assert cigar_2_parse[0] == [0, 12, 33, 132]
