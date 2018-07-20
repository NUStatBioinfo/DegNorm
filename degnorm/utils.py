import multiprocessing as mp
import logging
import sys
import os
import numpy as np
import argparse

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


def parse_args():
    """
    Obtain degnorm CLI parameters.

    :return: argparse.ArgumentParser object with runtime parameters required to run DegNorm pipeline.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-f'
                        , '--input-files'
                        , nargs='+'
                        , default=None
                        , required=True
                        , help='Input data files (an be multiple; separate with space).'
                               'Can be .bam or .sam files containing paired reads from RNA-seq experiments.')
    parser.add_argument('-d'
                        , '--input-dir'
                        , default=None
                        , required=False
                        , help='Input data directory. Use if not specifying individual .sam or .bam files.'
                               'Must be used in coordination with --input-type flag to specify .sam or .bam file type.')
    parser.add_argument('-t'
                        , '--input-type'
                        , default=None
                        , required=False
                        , choices=['bam', 'sam']
                        , help='Input RNA-seq experiment file type. Can be one of "bam" or "sam".'
                               'Only in coordination with --input-dir, i.e. when specifying an input data directory.')
    parser.add_argument('-g'
                        , '--genome-annotation'
                        , type=str
                        , default=None
                        , required=True
                        , help='Genome annotation file.'
                               'Must have extension .gtf or .gff.'
                               'All non-exon regions will be removed, along with exons that appear in'
                               'multiple chromosomes and exons that overlap with multiple genes.')
    parser.add_argument('-o'
                        , '--output-dir'
                        , type=str
                        , default=os.getcwd()
                        , required=False
                        , help='Output directory.'
                               'A directory for storing DegNorm analyses, visualizations, '
                               'and data will be created. Default to the current working directory.')
    parser.add_argument('-c'
                        , '--save-coverage'
                        , type=bool
                        , default=True
                        , help='Save coverage matrices flag.'
                               'If True, save coverage matrices per chromosome to .pkl files.'
                               'Default to True.')
    parser.add_argument('-c'
                        , '--cpu'
                        , type=int
                        , default=max_cpu()
                        , help='Number of cores for running DegNorm pipeline in parallel.'
                               'Defaults to the number of available cores - 1.')
    parser.add_argument('-v'
                        , '--version'
                        , action='version'
                        , version='DegNorm version 0.0.1'
                        , help='Display DegNorm package version and exit.')
    parser.add_argument('-h'
                        , '--help'
                        , action='help'
                        , default=argparse.SUPPRESS
                        , help='RNA-seq Degradation Normalization (DegNorm) pipeline package.'
                               'Accompanies our 2018 research paper "Normalization of generalized '
                               'transcript degradation improves accuracy in RNA-seq analysis."')

    args = parser.parse_args()

    # check validity of cores selection.
    avail_cores = max_cpu() + 1
    if args.cpu > avail_cores:
        logging.warning('{0} is greater than the number of available cores. Reducing to {1}'
                        .format(args.cpu, avail_cores))
        args.cpu = avail_cores

    # check validity of file i/o selection.
    if args.input_dir:
        if not args.input_type:
            raise ValueError('If input-dir is specified, you must also specify input-type (.sam or .bam files).')
        if not os.path.isdir(args.input_dir):
            raise IOError('Cannot find input-dir {0}'.format(args.input_dir))

        input_files = list()
        for f in os.listdir(args.input_dir):
            if f.endswith(args.input_type):
                input_files.append(os.path.join(args.input_dir, f))

        if not input_files:
            raise ValueError('No {0} files found in input-dir {1}'.format(args.input_type, args.input_dir))

        args.input_files = input_files

    # ensure that input files are uniquely named.
    args.input_files = list(set(args.input_files))

    # ensure that all files can be found.
    for f in args.input_files + [args.genome_annotation]:
        if not os.path.isfile(f):
            raise IOError('File {0} not found.'.format(f))

    # ensure only .sam or .bam files were supplied.
    extensions = list(set(map(lambda x: x.split('.')[-1], args.input_files)))
    if len(extensions) > 1:
        raise ValueError('Selected files contain multiple file extension types.'
                         'Must supply exclusively .sam or .bam files.')
    args.input_type = extensions[0]

    # check validity of output directory.
    if not os.path.isdir(args.output_dir):
        raise IOError('Cannot find output-dir {0}'.format(args.output_dir))

    return args