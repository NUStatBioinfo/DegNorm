import pytest
import os
import subprocess
import glob
from random import choice
from shutil import rmtree

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------- #
# Whole-pipeline test
# ----------------------------------------------------- #
@pytest.fixture
def run_setup(request):

    # create temporary directory for test.
    salt = ''.join([choice(['a', 'b', 'c', 'd', 'e', 'f', 'g', '1', '2', '3', '4', '5', '6', '7']) for _ in range(20)])
    outdir = os.path.join('degnorm_test_' + salt)
    os.makedirs(outdir)

    def teardown():
        rmtree(outdir)

    request.addfinalizer(teardown)

    return outdir


def test_pipeline(run_setup):
    print('TESTING DegNorm PIPELINE')

    # run degnorm command with test data.
    cmd = 'degnorm -i {0} {1} -g {2} -o {3} --genes {4}'.format(os.path.join(THIS_DIR, 'data', 'hg_small_1.sam')
                                                                , os.path.join(THIS_DIR, 'data', 'hg_small_2.sam')
                                                                , os.path.join(THIS_DIR, 'data', 'chr1.gtf')
                                                                , run_setup
                                                                , os.path.join(THIS_DIR, 'data', 'genes.txt'))

    out = subprocess.run([cmd]
                         , shell=True)
    assert out.returncode == 0