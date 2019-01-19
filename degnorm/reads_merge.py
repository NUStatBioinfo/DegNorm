from pandas import read_csv, concat
from collections import OrderedDict
import pickle as pkl
import numpy as np
import os
import re
import gc


def merge_read_count_files(file_dict, chroms):
    """
    Merge set of RNA-Seq samples' chromosome gene coverage count files into one pandas.DataFrame with
    one row per gene, columns are `chr`, `gene`, <sample IDs>.

    See reads.BamReadsProcessor.coverage_read_counts method.

    Example:
    file_dict = {('sample123': ['data/read_counts_sample123_chr1.csv', 'data/read_counts_sample123_chr2.csv']),
                 ('sample124': ['data/read_counts_sample124_chr1.csv', 'data/read_counts_sample124_chr2.csv'])} ->

    +-----------+---------+-----------------+-----------------+
    |    chr    |  gene   |    sample123    |    sample124    |
    +===========+=========+=================+=================+
    |   chr1    |  ATX2   |      275        |       101       |
    +-----------+---------+-----------------+-----------------+
    |    ...    |   ...   |      ...        |       ...       |
    +-----------+---------+-----------------+-----------------+
    |   chr2    |  GET4   |      301        |       255       |
    +-----------+---------+-----------------+-----------------+

    :param files: OrderedDict, keys are RNA-Seq sample IDs, values are lists of str filenames
    of read count .csv files, each file containing one chromosome's gene read counts for sample ID.
    :param chroms: list of str names of chromosomes for which to load read counts
    :return: pandas.DataFrame containing gene read counts across samples. Columns are `chr` (chromosome), `gene`,
    <sample IDs>
    """
    sample_ids = list(file_dict.keys())
    chrom_df_list = list()

    for chrom in chroms:
        for i in range(len(sample_ids)):

            # identify one (chromosome, sample ID) combination.
            sample_id = sample_ids[i]
            r = re.compile('read_counts_{0}_{1}.csv'.format(sample_id, chrom))
            counts_file = list(filter(r.search, file_dict[sample_id]))[0]

            # load sample's chromosome's read counts.
            sample_chrom_counts_df = read_csv(counts_file)

            # join together other samples' read counts for this chromosome.
            if i == 0:
                chrom_counts_df = sample_chrom_counts_df
            else:
                chrom_counts_df = chrom_counts_df.merge(sample_chrom_counts_df
                                                        , on='gene')

        # append chromosome column.
        chrom_counts_df['chr'] = chrom

        # save this chromosome's read count DataFrame with consistent ordering of column names.
        chrom_df_list.append(chrom_counts_df[['chr', 'gene'] + sample_ids])

    # vertically stack chromosomes read count DataFrames.
    chrom_counts_df = concat(chrom_df_list)

    # remove separated chromosome read count DataFrames.
    del chrom_df_list
    gc.collect()

    return chrom_counts_df


def merge_gene_coverage_files(file_dict, chroms, save_dir=None):
    """
    Merge set of RNA-Seq samples' chromosome gene coverage arrays into full coverage matrices across
    genes, across samples.

    See reads.BamReadsProcessor.coverage_read_counts method.

    Example:
    file_dict = {('sample123': ['data/coverage_sample123_chr1.csv', 'data/coverage_sample123_chr2.csv']),
                 ('sample124': ['data/coverage_sample124_chr1.csv', 'data/coverage_sample124_chr2.csv'])} ->

    {('gene A'): <L1 x 2 coverage array>,
     ('gene B'): <L2 x 2 coverage array>,
     ...
     ('gene N'): <LN x 2 coverage array>}

    :param files: OrderedDict, keys are RNA-Seq sample IDs, values are lists of str filenames
    of coverage .npz files, each file containing one chromosomes' genes's 1-d coverage array for sample ID.
    :param chroms: list of str names of chromosomes for which to load coverage matrices.
    :param save_dir: (optional) str, if specified save each chromosome's coverage arrays in a dictionary to a
    serialized pickle file, <save_dir>/<chromosome>/coverage_matrices_<chromosome>.pkl
    :return: OrderedDict of (gene, 2-d numpy coverage array) key-value pairs, one pair for each gene
    in genome.
    """
    sample_ids = list(file_dict.keys())
    gene_cov_dict = OrderedDict()

    for chrom in chroms:
        for i in range(len(sample_ids)):

            # identify one (chromosome, sample ID) combination.
            sample_id = sample_ids[i]
            r = re.compile('coverage_{0}_{1}.pkl'.format(sample_id, chrom))
            cov_file = list(filter(r.search, file_dict[sample_id]))[0]

            # load sample's chromosome's gene coverage arrays.
            with open(cov_file, 'rb') as f:
                cov_dat = pkl.load(f)

            # join together other samples' read counts for this chromosome.
            for gene in cov_dat:
                cov_array = cov_dat[gene]

                # initialize wide-format coverage matrix if on first sample.
                if i == 0:
                    gene_cov_dict[gene] = np.zeros([len(sample_ids), len(cov_array)])

                # load coverage array into appropriate column.
                gene_cov_dict[gene][i, :] = cov_array

        # if saving coverage arrays, save chromosome's data to separate directories.
        if save_dir:
            chrom_dir = os.path.join(save_dir, chrom)

            if not os.path.isdir(chrom_dir):
                os.makedirs(chrom_dir)

            chrom_file = os.path.join(chrom_dir, 'coverage_matrices_{0}.pkl'.format(chrom))
            with open(chrom_file, 'wb') as f:
                pkl.dump({gene: gene_cov_dict[gene] for gene in cov_dat}, f)

    return gene_cov_dict
