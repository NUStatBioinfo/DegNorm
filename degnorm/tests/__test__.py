import argparse
import os
import pytest
import pkg_resources

def parse_args():
    """
    Obtain degnorm CLI parameters.

    :return: argparse.ArgumentParser object with runtime parameters required to run DegNorm pipeline.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--keep-output'
                        , action='store_true'
                        , required=False
                        , help='If specified keep degnorm_test output directory')

    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    # determine whether or not to remove degnorm_test output directory after testing pipeline.
    os.environ['DEGNORM_TEST_CLEANUP'] = 'True'
    if args.keep_output:
        os.environ['DEGNORM_TEST_CLEANUP'] = 'False'

    # run existing unit tests.
    print('RUNNING DegNorm TESTS...')
    tests_dir = pkg_resources.resource_filename('degnorm', 'tests')
    pytest.main(['-x', tests_dir])
    print('ALL TESTS PASSING.')


if __name__ == "__main__":
    main()