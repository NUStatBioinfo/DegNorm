from pandas import read_csv, concat
from collections import OrderedDict
from degnorm.utils import *
from scipy import sparse
from joblib import Parallel, delayed
import pickle as pkl
import numpy as np
import os
import re
import gc
import tqdm


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


def dice_chrom_coverage(file_dict, chrom_exon_df, output_dir=None, verbose=True):
    """
    Slice up a coverage matrix for a chromosome into a dictionary of per-gene
    coverage matrices based on exon positioning for that chromosome.

    :param file_dict: OrderedDict, keys are RNA-Seq sample IDs, values are lists of str filenames
    of coverage .npz files, each file containing one chromosomes' genes's 1-d coverage array for sample ID.
    :param chrom_exon_df: pandas.DataFrame outlining exon positions within a single chromosome; has columns 'chr',
    'start' (exon start), 'end' (exon end), 'gene' (gene name), 'gene_end', and 'gene_start'
    :param output_dir: str (optional) if specified, save gene-level coverage matrices to a .pkl file
    under a chromosome sub-directory of output_dir.
    :param verbose: bool indicator should progress be written with logger?
    :return: dictionary of the form {gene_name: coverage numpy array} for specified chromosome
    """
    # output storage: (gene name, coverage matrix) key-value pairs.
    gene_cov_dict = dict()

    # identify the chromosome in question.
    unique_chrom = chrom_exon_df.chr.unique()
    if len(unique_chrom) > 1:
        raise ValueError('chrom_exon_df contains exon data for more than one chromosome!')

    chrom = unique_chrom[0]

    # Keep memory manageable:
    # attempt to break genes into groups so that each group's total coverage matrix
    # is ~ 500Mb. mem_splits dictates size of gene groups for load procedure.
    # use a randomly sampled chromosome coverage file to determine size of gene groups.
    random_sample_id = np.random.choice(list(file_dict.keys())
                                        , size=1)[0]
    r = re.compile('coverage_{0}_{1}.npz'.format(random_sample_id, chrom))
    npz_file = list(filter(r.search, file_dict[random_sample_id]))[0]
    cov_vec_sp = sparse.load_npz(npz_file)
    mem_splits = int(np.ceil(len(file_dict) * cov_vec_sp.asfptype().todense().nbytes / 500e6))
    del cov_vec_sp

    # sort genes by end position so we won't need to have entire chromosome coverage vectors loaded at once,
    # then break genes up into mem_splits subsets.
    chrom_exon_df = chrom_exon_df.sort_values('gene_end'
                                              , axis=0)
    genes = chrom_exon_df['gene'].unique()
    n_genes = len(genes)

    # if no breaks (e.g. if set of genes very small), the breaks are [0, number of genes]
    if mem_splits == 1:
        gene_splits = [0, n_genes]
    else:
        gene_splits = np.linspace(0
                                  , stop=n_genes
                                  , num=mem_splits
                                  , endpoint=True
                                  , dtype=int).tolist()
        gene_splits = list(set(gene_splits))
        gene_splits.sort()

    mem_splits = len(gene_splits) - 1

    use_pbar = False
    if verbose:
        logging.info('CHROMOSOME {0}: begin processing coverage matrices. \n'
                     'Using {1} gene splits for memory efficiency.'.format(chrom, mem_splits))

        # Instantiate progress bar if parsing non-negligible number of genes. Update in intervals of 5%.
        use_pbar = n_genes > 100
        if use_pbar:
            gene_idx = 0
            pbar_step_size = int(np.ceil(n_genes / 10))
            pbar = tqdm.tqdm(total=100
                             , leave=False
                             , desc='CHROMOSOME {0}: gene coverage matrix progress'.format(chrom)
                             , unit='%')

    # create the coverage matrix for each subset of genes.
    for i in range(mem_splits):

        # subset exon data to current gene subset.
        sub_genes = genes[gene_splits[i]:gene_splits[i + 1]].tolist()
        sub_chrom_exon_df = chrom_exon_df[chrom_exon_df.gene.isin(sub_genes)]

        # determine gene span: we only need a subset of the chromosome's coverage for gene subset.
        start_pos = int(sub_chrom_exon_df.gene_start.min() - 1)
        end_pos = int(sub_chrom_exon_df.gene_end.max() + 1)

        # load up gene span's coverage matrix.
        idx = 0
        for sample_id in file_dict:
            r = re.compile('coverage_{0}_{1}.npz'.format(sample_id, chrom))
            npz_file = list(filter(r.search, file_dict[sample_id]))[0]

            # load the gene span of a sample's compressed sparse row chromosome coverage array
            cov_vec_sp = sparse.load_npz(npz_file).transpose()[start_pos:end_pos, :]

            # initialize coverage matrix (len(chrom) x p) with copy of first experiment's coverage array.
            if idx == 0:
                cov_mat = cov_vec_sp.copy()

            # column-append sparse coverage vector to existing coverage matrix.
            else:
                cov_mat = sparse.hstack([cov_mat, cov_vec_sp], dtype=int)

            idx += 1

        # clean up coverage array data.
        del cov_vec_sp
        gc.collect()

        # convert coverage matrix to dense matrix for speed in splicing,
        # should blow up coverage matrix size to about 500Mb, on average.
        cov_mat = cov_mat.asfptype().todense()

        # tear out each gene's coverage matrix from loaded chromosome coverage sub-matrix.
        for ii in range(len(sub_genes)):

            gene = sub_genes[ii]

            # subset chromosome's gene/exon data to a single gene and
            # identify gene's start and end positions.
            single_gene_df = sub_chrom_exon_df[sub_chrom_exon_df.gene == gene]

            # Slice up cov_mat based on relative exon positions within a gene while remembering to
            # shift starts and ends based on the start position of the current gene span.
            e_starts, e_ends = single_gene_df.start.values - start_pos - 1, single_gene_df.end.values - start_pos - 1
            slices = [np.arange(e_starts[j], e_ends[j]) for j in range(len(e_starts))]

            # in case exons are overlapping, take union of their covered regions.
            slicing = np.unique(flatten_2d(slices))

            # Save transposed coverage matrix so that shape is p x Li.
            gene_cov_dict[gene] = np.array(cov_mat[slicing, :].T).astype(np.float_)

            if use_pbar:
                if (gene_idx % pbar_step_size == 0) and (gene_idx > 0):
                    pbar.update(10)
                    gene_idx += 1

    # close progress bar if using one.
    if use_pbar:
        pbar.close()

    # free up memory allocation for dense coverage matrix, exon data subset.
    del cov_mat, sub_chrom_exon_df
    gc.collect()

    if verbose:
        logging.info('CHROMOSOME {0}: successfully parsed {1} gene coverage matrices.'
                     .format(chrom, len(gene_cov_dict)))

    # if a save location is specified, save {gene: coverage matrix} data per chromosome in a
    # new directory named after the chromosome.
    if output_dir:
        output_dir = os.path.join(output_dir, chrom)

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        # save per-gene coverage matrices to .pkl files
        gene_cov_file = os.path.join(output_dir, 'coverage_matrices_{0}.pkl'.format(chrom))
        if verbose:
            logging.info('CHROMOSOME {0}: Saving gene coverage matrices to {1}'
                         .format(chrom, gene_cov_file))

        with open(gene_cov_file, 'wb') as f:
            pkl.dump(gene_cov_dict, f)

    return gene_cov_dict


def merge_gene_coverage_files(file_dict, exon_df, n_jobs=max_cpu(),
                              output_dir=None, verbose=True):
    """
    For each chromosome, load the coverage arrays resulting from each alignment file, join them,
    and then slice the joined coverage array into per-gene coverage matrices. Run process in parallel over
    chromosomes with the n_jobs argument.

    Example:
    file_dict = {('sample123': ['data/coverage_sample123_chr1.npz', 'data/coverage_sample123_chr2.npz']),
                 ('sample124': ['data/coverage_sample124_chr1.npz', 'data/coverage_sample124_chr2.npz'])} ->

    {('gene A'): <L1 x 2 coverage array>,
     ('gene B'): <L2 x 2 coverage array>,
     ...
     ('gene N'): <LN x 2 coverage array>}

    :param file_dict: OrderedDict, keys are RNA-Seq sample IDs, values are lists of str filenames
    of read count .csv files, each file containing one chromosome's gene read counts for sample ID.
    :param exon_df: pandas.DataFrame outlining exon positions for an entire genome; has columns 'chr',
    'start' (exon start), 'end' (exon end), 'gene' (gene name), 'gene_end', and 'gene_start'
    :param n_jobs: int number of cores used for distributing gene coverage merge process over different chromosomes.
    :param output_dir: str (optional) if specified, save chromosome gene coverage matrix dictionaries
     to serialized .pkl files of the form `<output_dir>/<chromosome>/coverage_matrices_<chromosome>.pkl`
    :param verbose: bool indicator should progress be written with logger?
    :return: OrderedDict of the form {gene: 2-d numpy coverage array} for all genes present in exon_df.
    """
    chroms = exon_df.chr.unique()
    gene_cov_dict = OrderedDict()

    # get list of chromosome gene coverage dictionaries.
    chrom_gene_cov_dicts = Parallel(n_jobs=min(n_jobs, len(chroms))
                              , verbose=0
                              , backend='threading')(delayed(dice_chrom_coverage)(
        file_dict=file_dict,
        chrom_exon_df=subset_to_chrom(exon_df, chrom=chrom),
        output_dir=output_dir,
        verbose=verbose) for chrom in chroms)

    # concatenate chromosome coverage dictionaries into one ordered dictionary.
    for chrom_gene_cov_dict in chrom_gene_cov_dicts:
        for gene in chrom_gene_cov_dict:
            gene_cov_dict[gene] = chrom_gene_cov_dict[gene]

    del chrom_gene_cov_dicts, chrom_gene_cov_dict
    gc.collect()

    return gene_cov_dict


# if __name__ == '__main__':
    # from degnorm.gene_processing import GeneAnnotationProcessor
    #
    # data_path = '/Users/fineiskid/nu/jiping_research/degnorm_test_files'
    #
    # cov_files = [os.path.join(data_path, 'hg_small_1', 'coverage_hg_small_1_chr1.npz')
    #              , os.path.join(data_path, 'hg_small_2', 'coverage_hg_small_2_chr1.npz')]
    #
    # file_dict = {'hg_small_1': [cov_files[0]]
    #              , 'hg_small_2': [cov_files[1]]}
    #
    # gtf_file = '/Users/fineiskid/nu/jiping_research/DegNorm/degnorm/tests/data/chr1_small.gtf'
    #
    # gtf_processor = GeneAnnotationProcessor(gtf_file)
    # exon_df = gtf_processor.run()
    #
    # # dice_chrom_coverage(file_dict, chrom_exon_df=subset_to_chrom(exon_df, 'chr1')
    # #                     , output_dir=data_path, verbose=True)
    #
    # merge_gene_coverage_files(file_dict
    #                           , exon_df=exon_df
    #                           , n_jobs=1
    #                           , output_dir=data_path,
    #                           verbose=True)
