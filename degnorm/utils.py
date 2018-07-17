import multiprocessing as mp
import logging
import sys
import numpy as np

logging.basicConfig(stream=sys.stdout
                    , level=logging.DEBUG
                    , format='%(asctime)s ---- %(message)s'
                    , datefmt='%m/%d/%Y %I:%M:%S')


def subset_to_chrom(df, chrom, reindex=False):
    """
    Subset a pandas.DataFrame with a 'chr' (chromosome) column to a particular
    chromosome. Reset the index if desired.

    :param df: pandas.DataFrame with a 'chr' column
    :param chrom: str chromosome name
    :param reindex: bool indicator reset index? Default False: do not reset index.
    :return: pandas.DataFrame
    """
    if not reindex:
        return df[df['chr'] == chrom]
    else:
        return df[df['chr'] == chrom].reset_index(drop=True)


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