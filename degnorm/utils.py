import multiprocessing as mp
import logging
import sys
import numpy as np

logging.basicConfig(stream=sys.stdout
                    , level=logging.DEBUG
                    , format='%(asctime)s ---- %(message)s'
                    , datefmt='%m/%d/%Y %I:%M:%S')


def subset_to_chrom(df, chrom):
    return df[df['chr'] == chrom]


def max_cpu():
    return mp.cpu_count() - 1


def flatten_2d(lst2d):
    """
    Flatten a 2-dimensional list of lists or list of numpy arrays into a single numpy array.

    :param lst2d: 2-dimensional list of lists or list of numpy arrays
    :return: 1-dimensional numpy array
    """
    arr1d = np.array([elt for lst1d in lst2d for elt in lst1d])

    return arr1d