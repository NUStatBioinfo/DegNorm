import sys
import pytest
import pkg_resources


def main():

    # run existing unit tests.
    print('RUNNING TESTS...')
    tests_dir = pkg_resources.resource_filename('degnorm', 'tests')

    pytest.main(['-x', tests_dir])


if __name__ == "__main__":
    main()
