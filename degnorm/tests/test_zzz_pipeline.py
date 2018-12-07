import pytest
import os
import subprocess
import shutil
from random import choice

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

        # determine whether or not to remove pipeline test output.
        if os.environ['DEGNORM_TEST_CLEANUP'] == 'True':
            shutil.rmtree(outdir)

    request.addfinalizer(teardown)

    return outdir


def test_pipeline(run_setup):

    # run degnorm command with test data, with downsampling.
    cmd = 'degnorm --bam-files {0} {1} --bai-files {2} {3} -g {4} -o {5} --plot-genes {6} --nmf-iter 50 -d 50'\
        .format(os.path.join(THIS_DIR, 'data', 'hg_small_1.bam')
                , os.path.join(THIS_DIR, 'data', 'hg_small_2.bam')
                , os.path.join(THIS_DIR, 'data', 'hg_small_1.bai')
                , os.path.join(THIS_DIR, 'data', 'hg_small_2.bai')
                , os.path.join(THIS_DIR, 'data', 'chr1_small.gtf')
                , run_setup
                , os.path.join(THIS_DIR, 'data', 'genes.txt'))

    out = subprocess.run([cmd]
                         , shell=True
                         , stderr=subprocess.PIPE)
    assert out.returncode == 0

    # run degnorm command with test data, without downsampling, without specifying .bai files directly either.
    cmd = 'degnorm --bam-files {0} {1} -g {2} -o {3} --plot-genes {4} --nmf-iter 50'\
        .format(os.path.join(THIS_DIR, 'data', 'hg_small_1.bam')
                , os.path.join(THIS_DIR, 'data', 'hg_small_2.bam')
                , os.path.join(THIS_DIR, 'data', 'chr1_small.gtf')
                , run_setup
                , os.path.join(THIS_DIR, 'data', 'genes.txt'))

    out = subprocess.run([cmd]
                         , shell=True
                         , stderr=subprocess.PIPE)
    assert out.returncode == 0


def test_mpi_pipeline(run_setup):

    # run degnorm_mpi command with test data, with downsampling.
    cmd = 'mpiexec -n 2 degnorm_mpi --bam-files {0} {1} -g {2} -o {3} -p 2 --plot-genes {4} --nmf-iter 50'\
        .format(os.path.join(THIS_DIR, 'data', 'hg_small_1.bam')
                , os.path.join(THIS_DIR, 'data', 'hg_small_2.bam')
                , os.path.join(THIS_DIR, 'data', 'chr1_small.gtf')
                , run_setup
                , os.path.join(THIS_DIR, 'data', 'genes.txt'))

    # if MPI available and > 1 node available, run pipeline test for degnorm_mpi.
    try:
        from mpi4py import MPI

        # try running degnorm_mpi test if there are > 1 nodes available. O/w skip it.
        try:
            out = subprocess.run([cmd]
                                 , shell=True
                                 , stderr=subprocess.PIPE)

            assert out.returncode == 0

        except RuntimeError:
            pass

    except ImportError:
        pass