import multiprocessing as mp
import logging
import sys
import numpy as np
import subprocess
import os
from datetime import datetime
import time
import argparse
import pkg_resources
import gc


logging.basicConfig(stream=sys.stdout
                    , level=logging.DEBUG
                    , format='DegNorm (%(asctime)s) ---- %(message)s'
                    , datefmt='%m/%d/%Y %I:%M:%S')

# welcome our user.
def welcome():
    resources_dir = pkg_resources.resource_filename('degnorm', 'resources')
    with open(os.path.join(resources_dir, 'welcome.txt'), 'r') as f:
        welcome = f.readlines()
        welcome += '\nversion {0}'.format(pkg_resources.get_distribution('degnorm').version)

    sys.stdout.write(''.join(welcome) + '\n'*4)


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
    return mp.cpu_count() - 1


def flatten_2d(lst2d):
    """
    Flatten a 2-dimensional list of lists or list of numpy arrays into a single numpy array.

    :param lst2d: 2-dimensional list of lists or list of numpy arrays
    :return: 1-dimensional numpy array
    """
    arr1d = np.array([elt for lst1d in lst2d for elt in lst1d])

    return arr1d


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


def bam_to_sam(bam_file):
    """
    Convert a .bam file to a .sam file with samtools.

    :param bam_file: str realpath to .bam file to be converted to .sam file format.
    :return: str realpath to the created .sam file
    """
    if not bam_file.endswith('.bam'):
        raise ValueError('{0} is not a .bam file'.format(bam_file))

    output_dir = os.path.dirname(bam_file)
    sam = os.path.basename(bam_file).split('.')[:-1][0] + '.sam'
    sam_file = os.path.join(output_dir, sam)

    # check if a .sam file already exists; if it does, add salt by time.
    while os.path.isfile(sam_file):
        time.sleep(2)
        sam = os.path.basename(bam_file).split('.')[:-1][0] + '_' + datetime.now().strftime('%m%d%Y_%H%M%S') + '.sam'
        sam_file = os.path.join(output_dir, sam)

    cmd = 'samtools view -h -o {0} {1}'.format(sam_file, bam_file)
    out = subprocess.run([cmd], shell=True)

    if out.returncode != 0:
        raise ValueError('{0} was not successfully converted into a .sam file'.format(bam_file))

    return sam_file


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


def parse_args():
    """
    Obtain degnorm CLI parameters.

    :return: argparse.ArgumentParser object with runtime parameters required to run DegNorm pipeline.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-i'
                        , '--input-files'
                        , nargs='+'
                        , default=None
                        , required=False
                        , help='Input data files (an be multiple; separate with space).'
                               'Can be .bam or .sam files containing paired reads from RNA-seq experiments.')
    parser.add_argument('--input-dir'
                        , default=None
                        , required=False
                        , help='Input data directory. Use instead of, or in addition to, specifying individual '
                               '.sam files. All .sam files in this directory will be considered input to DegNorm.')
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
                        , help='Gene nucleotide downsample rate for systematic sampling of gene coverage curves.'
                               'Specifies a \'take every\' interval. Larger value -> fewer bases sampled. '
                               'Resulting coverage matrix estimates (and plots) will be re-interpolated back '
                               'to original size.'
                               'Use to speed up computation. Default is NO downsampling.')
    parser.add_argument( '--nmf-iter'
                        , type=int
                        , default=100
                        , required=False
                        , help='Number of iterations to perform per NMF-OA computation per gene.'
                               'Different than number of DegNorm iterations (--iter flag). Default = 100.')
    parser.add_argument( '--iter'
                        , type=int
                        , default=5
                        , required=False
                        , help='Number of DegNorm iterations to perform. Default = 5.'
                               'Different than number of NMF-OA iterations (--nmf-iter flag).')
    parser.add_argument('--minimax-coverage'
                        , type=int
                        , default=20
                        , required=False
                        , help='Minimum maximum read coverage for a gene to be included in DegNorm Pipeline. '
                               'Default is 20.')
    parser.add_argument('-c'
                        , '--cpu'
                        , type=int
                        , default=max_cpu()
                        , help='Number of cores for running DegNorm pipeline in parallel.'
                               'Defaults to the number of available cores - 1.')
    parser.add_argument('-t'
                        , '--input-type'
                        , default='sam'
                        , required=False
                        , choices=['bam', 'sam']
                        , help='Input RNA-seq experiment file type. Can be one of "bam" or "sam".'
                               'Only in coordination with --input-dir, i.e. when specifying an input data directory.')
    parser.add_argument('-v'
                        , '--version'
                        , action='version'
                        , version='DegNorm version {0}'.format(pkg_resources.get_distribution('degnorm').version)
                        , help='Display DegNorm package version and exit.')
    parser.add_argument('-h'
                        , '--help'
                        , action='help'
                        , default=argparse.SUPPRESS
                        , help='RNA-seq Degradation Normalization (DegNorm) pipeline package.'
                               'Accompanies our 2018 research paper "Normalization of generalized '
                               'transcript degradation improves accuracy in RNA-seq analysis."')

    args = parser.parse_args()

    if not args.input_files and not args.input_dir:
        raise ValueError('Must specify one of input-files or input-dir')

    # check validity of cores selection.
    avail_cores = max_cpu() + 1
    if args.cpu > avail_cores:
        logging.warning('{0} is greater than the number of available cores. Reducing to {1}'
                        .format(args.cpu, avail_cores))
        args.cpu = avail_cores

    # check validity of file i/o selection.
    if args.input_dir:
        if not os.path.isdir(args.input_dir):
            raise IOError('Cannot find input-dir {0}'.format(args.input_dir))

        # scan directory for .sam files.
        input_files = list()
        for f in os.listdir(args.input_dir):
            if f.endswith('.sam'):
                input_files.append(os.path.join(args.input_dir, f))

        if not input_files:
            raise ValueError('No {0} files found in input-dir {1}'.format(args.input_type, args.input_dir))

        # if user used -i/--input-files, append contents of directory to individually specified files.
        if args.input_files:
            args.input_files += input_files

        else:
            args.input_files = input_files

    # ensure that input files are uniquely named.
    args.input_files = list(set(args.input_files))

    # ensure that all files can be found.
    for f in args.input_files:
        if not os.path.isfile(f):
            raise IOError('Input file {0} not found.'.format(f))

    if not os.path.isfile(args.genome_annotation):
        raise IOError('Gene annotation file {0} not found.'.format(args.genome_annotation))

    # ensure there are at least 2 experiment files.
    if len(args.input_files) == 1:
        raise ValueError('Must input >= 2 unique RNA-Seq experiment files! Cannot estimate coverage curve matrix '
                         'approximations from a single experiment.')

    # if .bam files supplied, make sure samtools is installed.
    if any([x.split('.')[-1] == 'bam' for x in args.input_files]):
        samtools_avail = find_software('samtools')

        # Note that samtools is only available for Linux and Mac OS:
        # https://github.com/samtools/samtools/blob/develop/INSTALL
        if not samtools_avail:
            raise EnvironmentError('samtools is not installed or is not in your PATH.'
                                   'samtools is required to convert .bam -> .sam files'
                                   'Either use .sam files or install samtools.')

    # check validity of output directory.
    if not os.path.isdir(args.output_dir):
        raise IOError('Cannot find output-dir {0}'.format(args.output_dir))

    # if --genes is specified, parse input in case file(s) are given.
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

    # quality control on DegNorm parameters.
    if (args.nmf_iter < 1) or (args.iter < 1) or (args.downsample_rate < 1):
        raise ValueError('--nmf-iter, --iter, and --downsample-rate must all be >= 1.')

    welcome()

    return args


if __name__ == '__main__':
    print(parse_args())