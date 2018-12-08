import multiprocessing as mp
import logging
import warnings
import sys
import re
import numpy as np
import subprocess
import os
from datetime import datetime
import time
import argparse
import pkg_resources
import gc


def configure_logger(output_dir=None, mpi=False):
    """
    Configure DegNorm logger. Save log to file in output directory and route to stdout.

    :param output_dir: str path to DegNorm run output dir where degnorm.log file to be written.
    :param mpi: Bool is user running degnorm_mpi?
    """
    handlers = [logging.StreamHandler()]
    if output_dir:
        handlers += [logging.FileHandler(os.path.join(output_dir, 'degnorm.log'))]

    fmt = 'DegNorm (%(asctime)s) ---- %(message)s'
    if mpi:
        fmt = 'DegNorm MPI (%(asctime)s) ---- %(message)s'

    logging.basicConfig(level=logging.DEBUG
                        , format=fmt
                        , handlers=handlers
                        , datefmt='%m/%d/%Y %I:%M:%S')


def welcome():
    """
    Welcome our user with DegNorm ascii art.
    """
    resources_dir = pkg_resources.resource_filename('degnorm', 'resources')
    with open(os.path.join(resources_dir, 'welcome.txt'), 'r') as f:
        welcome = f.readlines()
        welcome += '\n' + 'version {0}'.format(pkg_resources.get_distribution('degnorm').version)

    logging.info('\n' + ''.join(welcome) + '\n'*4)


def create_output_dir(output_dir):
    """
    Create a DegNorm output directory.

    :param output_dir: str desired path to DegNorm output directory
    :return: str path to newly created output directory
    """
    dir_to_make = os.path.join(output_dir, 'DegNorm_' + datetime.now().strftime('%m%d%Y_%H%M%S'))
    while os.path.isdir(dir_to_make):
        time.sleep(2)
        dir_to_make = os.path.join(output_dir, 'DegNorm_' + datetime.now().strftime('%m%d%Y_%H%M%S'))

    os.makedirs(dir_to_make)
    return dir_to_make


def subset_to_chrom(df, chrom, reindex=False):
    """
    Subset a pandas.DataFrame with a 'chr' (chromosome) column to a particular
    chromosome. Reset the index if desired.

    :param df: pandas.DataFrame with a 'chr' column
    :param chrom: str or list; chromosome(s) with which to subset df
    :param reindex: bool indicator reset index? Default False: do not reset index.
    :return: pandas.DataFrame
    """
    if not isinstance(chrom, list):
        chrom = [chrom]

    if not reindex:
        sub_df = df[df['chr'].isin(chrom)]
    else:
        sub_df = df[df['chr'].isin(chrom)].reset_index(drop=True)

    if sub_df.empty:
        raise ValueError('Chromosome subsetting resulted in an empty DataFrame!')

    return sub_df


def max_cpu():
    """
    :return: int number of CPUs available on a node minus 1
    """
    return mp.cpu_count() - 1


def flatten_2d(lst2d, arr=True):
    """
    Flatten a 2-dimensional list of lists or list of numpy arrays into a single list or numpy array.

    :param lst2d: 2-dimensional list of lists or list of 1-d numpy arrays
    :param arr: Bool return numpy array or list?
    :return: 1-dimensional list or numpy array
    """
    lst1d = [elt for lst1d in lst2d for elt in lst1d]
    return np.array(lst1d) if arr else lst1d


def find_software(software='samtools'):
    """
    Determine if a software is in a $PATH.

    :return: True if software is in path, otherwise False
    """
    out = subprocess.run(['which {0}'.format(software)]
                         , shell=True)
    if out.returncode != 0:
        return False

    return True


def bai_from_bam_file(bam_file):
    """
    Simple helper function to change the file extension of a .bam file to .bai.
    """
    if not bam_file.endswith('.bam'):
        raise ValueError('{0} must have a .bam extension.'.format(bam_file))

    return bam_file[:-3] + 'bai'


def create_index_file(bam_file):
    """
    Create a BAM index file with samtools.

    :param bam_file: str realpath to .bam file for which we desire a .bai index file
    :return: str realpath to the created .bai file
    """
    bai_file = bai_from_bam_file(bam_file)

    # check if samtools is available in $PATH.
    samtools_avail = find_software('samtools')

    # Note that samtools is only available for Linux and Mac OS:
    # https://github.com/samtools/samtools/blob/develop/INSTALL
    if not samtools_avail:
        raise EnvironmentError('samtools not in found in PATH. samtools is required to convert {0} -> {1}'
                               .format(bam_file, bai_file))

    # run samtools index
    cmd = 'samtools index {0} {1}'.format(bam_file, bai_file)
    out = subprocess.run([cmd]
                         , shell=True)

    if out.returncode != 0:
        raise ValueError('{0} was not successfully converted into a .bai file'.format(bam_file))


def split_into_chunks(x, n):
    """
    Split a list into a set of n approximately evenly-sized sublists.

    :param x: iterable - a list, numpy array, tuple, etc. Not a generator.
    :param n: int number of chunks. If n >= len(x), split x into sublists of size 1.
    :return: list of lists
    """
    csize = int(np.ceil(len(x) / n))
    out = list()
    
    i = 0
    while i*csize < len(x):
        out.append(x[(i * csize):(i * csize + csize)])
        i += 1

    return out


def argparser():
    """
    Obtain degnorm CLI parameters.

    :return: argparse.ArgumentParser object with runtime parameters required to run DegNorm pipeline.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--bam-files'
                        , nargs='+'
                        , default=None
                        , required=False
                        , help='Aligned read files (can be multiple; separate with space) in .bam format. '
                               '.bam files may be for paired or single-end read experiments.')
    parser.add_argument('--bai-files'
                        , nargs='+'
                        , default=None
                        , required=False
                        , help='Input .bam index files (Can be multiple; separate with space). '
                               'Only use with --bam-files.'
                               '.bai files for files passed to --bam-files. Assumed .bai file order corresponds to '
                               'supplied .bam files. If --bam-files is supplied and not --bai-files, it is assumed '
                               'that you have .bai files in the same location as --bam-files, and with the same '
                               'basename, e.g. (A01.bam, A01.bai), (sample_27.bam, sample_27.bai)')
    parser.add_argument('--bam-dir'
                        , default=None
                        , required=False
                        , help='Input .bam/.bai data directory. Use instead of, or in addition to, specifying individual '
                               '.bam files. All .bam files (and .bai files with the same basename) in this directory '
                               'will be considered input to DegNorm.')
    parser.add_argument('-w'
                        , '--warm-start-dir'
                        , default=None
                        , required=False
                        , help='Previous DegNorm run output directory. Use to source '
                               'gene coverage matrices, read counts, and parsed genome annotation data. '
                               'Use this to avoid duplicating costly preprocessing over multiple runs.')
    parser.add_argument('-g'
                        , '--genome-annotation'
                        , type=str
                        , default=None
                        , required=False
                        , help='Genome annotation file.'
                               'Must have extension .gtf or .gff.'
                               'All non-exon regions will be removed, along with exons that appear in '
                               'multiple chromosomes and exons that overlap with multiple genes.')
    parser.add_argument('-o'
                        , '--output-dir'
                        , type=str
                        , default=os.getcwd()
                        , required=False
                        , help='Output directory.'
                               'A directory for storing DegNorm analyses, visualizations, '
                               'and data will be created. Default to the current working directory.')
    parser.add_argument('--plot-genes'
                        , nargs='+'
                        , type=str
                        , default=None
                        , required=False
                        , help='List of gene names or a text file (with extension .txt) specifying a set'
                               'of genes for which you would like pre- and post-DegNorm coverage curve plots rendered.')
    parser.add_argument('-d'
                        , '--downsample-rate'
                        , type=int
                        , default=1
                        , required=False
                        , help='Gene nucleotide downsample rate for systematic sampling of gene coverage curves. '
                               'Integer-valued. Specifies a \'take every\' interval. '
                               'Larger value -> fewer bases sampled. '
                               'Resulting coverage matrix estimates (and plots) will be re-interpolated back '
                               'to original size. '
                               'Use to speed up computation. Default is NO downsampling (i.e. downsample rate == 1.')
    parser.add_argument('--nmf-iter'
                        , type=int
                        , default=100
                        , required=False
                        , help='Number of iterations to perform per NMF-OA computation per gene. '
                               'Different than number of DegNorm iterations (--iter flag). Default = 100.')
    parser.add_argument('--iter'
                        , type=int
                        , default=5
                        , required=False
                        , help='Number of DegNorm iterations to perform. Default = 5. '
                               'Different than number of NMF-OA iterations (--nmf-iter flag).')
    parser.add_argument('--minimax-coverage'
                        , type=int
                        , default=1
                        , required=False
                        , help='Minimum maximum read coverage for a gene to be included in DegNorm Pipeline. ')
    parser.add_argument('-s'
                        , '--skip-baseline-selection'
                        , action='store_true'
                        , help='Skip baseline selection while computing coverage matrix estimates. '
                               'This will speed up degradation index score computation but may make '
                               'scores less accurate.')
    parser.add_argument('-u'
                        , '--unique-alignments'
                        , action='store_true'
                        , help='Only retain reads that were uniquely aligned. All reads with '
                               'the flag "NH:i:<x>" with x > 1 will be dropped.')
    parser.add_argument('-p'
                        , '--proc-per-node'
                        , type=int
                        , required=False
                        , default=max_cpu()
                        , help='Number of processes to spawn per node, for within-node parallelization.'
                               'Defaults to the number of available cores (on the worker node) - 1. '
                               'If too high for a given node, reduce to the node\'s default.')
    parser.add_argument('-v'
                        , '--version'
                        , action='version'
                        , version='DegNorm version {0}'.format(pkg_resources.get_distribution('degnorm').version)
                        , help='Display DegNorm package version and exit.')
    parser.add_argument('-h'
                        , '--help'
                        , action='help'
                        , default=argparse.SUPPRESS
                        , help='RNA-seq Degradation Normalization (DegNorm) pipeline package. '
                               'Accompanies our 2018 research paper "Normalization of generalized '
                               'transcript degradation improves accuracy in RNA-seq analysis."')

    return parser


def parse_args(mpi=False):
    """
    Parse command line arguments.

    :param mpi: Boolean is DegNorm being run in MPI mode?
    :return: parsed argparse.ArgumentParser
    """
    parser = argparser()
    args = parser.parse_args()

    # check validity of cores selection.
    if not mpi:
        max_ppn = max_cpu() + 1
        if args.proc_per_node > max_ppn:
            warnings.warn('{0} is greater than the number of available cores ({1}). Reducing to {2}'
                          .format(args.proc_per_node, max_ppn, max_ppn - 1))
            args.proc_per_node = max_ppn - 1

    # check validity of output directory.
    if not os.path.isdir(args.output_dir):
        raise NotADirectoryError('Cannot find output directory {0} for saving output'.format(args.output_dir))

    # ensure that user has supplied fresh .bam/.bai files or a warm start directory.
    if (not args.bam_files and not args.bam_dir) and (not args.warm_start_dir):
        raise ValueError('Must specify either --bam-files, --bam-dir, or --warm-start-dir as a data input option.')

    # quality control on DegNorm parameters.
    if (args.nmf_iter < 1) or (args.iter < 1) or (args.downsample_rate < 1):
        raise ValueError('--nmf-iter, --iter, and --downsample-rate must all be >= 1.')

    # if --plot-genes is specified, parse input for any .txt file(s) in addition to possible cli-specified genes.
    if args.plot_genes:
        genes = list()

        # check for text files that will hold names of genes.
        gene_files = list(filter(lambda x: x.endswith('.txt'), args.plot_genes))
        for gene_file in gene_files:
            with open(gene_file, 'r') as f:
                genes_in_file = f.readlines()
                genes += [x.strip() for x in genes_in_file]

        # include any individually-specified genes.
        genes += list(set(args.plot_genes) - set(gene_files))

        # replace input with parsed list of unique genes.
        args.plot_genes = list(set(genes))

    # basic checks on warm start directory if supplied.
    if args.warm_start_dir:

        if not os.path.isdir(args.warm_start_dir):
            raise NotADirectoryError('Cannot find --warm-start-dir {0}'.format(args.warm_start_dir))

        # warn user if they have also supplied read alignments and/or a .gtf file, that they will
        # be ignored for the content within the warm start directory.
        if args.bam_files or args.bam_dir or args.genome_annotation:
            logging.warning('Using warm-start directory. Supplied .bam files, .bam directory, '
                            'and genome annotation file will be ignored.')

            args.bam_files = None
            args.bai_files = None
            args.create_bai_files = None
            args.bam_dir = None
            args.genome_annotation = None

    # if not using a warm-start, parse input RNA-Seq + genome annotation files.
    else:

        # check validity of gene annotation file selection.
        if not args.genome_annotation:
            raise ValueError('If warm-start directory not specified, gene annotation file must be specified!')

        else:
            if not os.path.isfile(args.genome_annotation):
                raise FileNotFoundError('Gene annotation file {0} not found.'.format(args.genome_annotation))

        # check validity of file i/o selection.
        bam_files = list()
        bai_files = list()
        create_bai_files = list()

        # INPUT OPTION 1: a --bam-dir was specified.
        if args.bam_dir:

            # if user used both --bam-dir and --bam-files and/or --bai-files, yell at them. (only use one method).
            if args.bam_files is not None or args.bai_files is not None:
                raise ValueError('Do not specify both a --bam-dir and either --bam-files and/or --bai-files.'
                                 'Use one input selection method or the other.')

            # check that the dir actually exists.
            if not os.path.isdir(args.bam_dir):
                raise NotADirectoryError('Cannot find --bam-dir {0}'.format(args.bam_dir))

            # scan directory for .bam files.
            for f in os.listdir(args.bam_dir):
                if f.endswith('.bam'):
                    bam_files.append(os.path.join(args.bam_dir, f))

            # check that we found enough .bam files.
            if len(bam_files) < 2:
                raise ValueError('Only found {0} .bam {1} within directory {2}'
                                 .format(len(bam_files), 'files' if len(bam_files) > 0 else 'files', args.bam_dir))

            # search for .bai files in the --bam-dir. If they don't exist, try to make them.
            for bam_file in bam_files:
                bai_file = re.sub('.bam$', '.bai', bam_file)

                # if .bai file under same basename as .bam file doesn't exist,
                # add it to list of .bai files that need to be created.
                if not os.path.isfile(bai_file):
                    bai_files.append(bai_from_bam_file(bam_file))
                    create_bai_files.append(bam_file)
                else:
                    bai_files.append(bai_file)

        # INPUT OPTION 2: --bam-files and possibly --bai-files were specified.
        else:
            # ensure .bam files are actually .bam files.
            for bam_file in args.bam_files:
                if not bam_file.endswith('.bam'):
                    raise ValueError('{0} is not a .bam file.'.format(bam_file))
                elif not os.path.isfile(bam_file):
                    raise FileNotFoundError('Count not find .bam file {0}'.format(bam_file))
                else:
                    bam_files.append(bam_file)

            # case where user has specified .bai files to accompany .bam files.
            if args.bai_files is not None:
                # if user has supplied an incorrect number of bai files, fail out.
                if len(args.bai_files) != len(bam_files):
                    raise ValueError('Number of supplied .bai files does not match number of supplied .bam files.')

                # ensure .bai files are actually .bai files.
                for bai_file in args.bai_files:
                    if not bai_file.endswith('.bai'):
                        raise ValueError('{0} is not a .bai file.'.format(bai_file))
                    elif not os.path.isfile(bai_file):
                        raise FileNotFoundError('Count not find .bai file {0}'.format(bai_file))
                    else:
                        bai_files.append(bai_file)

            # if user has not supplied any bai files: look for them under the same name
            # as each of the .bam files, or create new .bai files with samtools (if possible).
            else:
                for bam_file in bam_files:
                    bai_file = re.sub('.bam$', '.bai', bam_file)

                    # if .bai file under same name as .bam file doesn't exist,
                    # add it to list of .bam files for which we need to create a .bai file.
                    if not os.path.isfile(bai_file):
                        bai_files.append(bai_from_bam_file(bam_file))
                        create_bai_files.append(bam_file)
                    else:
                        bai_files.append(bai_file)

        # ensure we have at least two .bam files.
        if len(bam_files) < 2:
            raise ValueError('Fewer than 2 .bam files were found. Not sufficiently many to run DegNorm.')

        # ensure that input files are uniquely named.
        if len(bam_files) != len(set(bam_files)):
            raise ValueError('Supplied .bam files are not uniquely named!')

        # create parser attributes for bam/index files.
        args.bam_files = bam_files
        args.bai_files = bai_files
        args.create_bai_files = create_bai_files

    return args