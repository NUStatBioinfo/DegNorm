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


def merge_read_counts(data_dir, sample_ids, chroms):
    """
    Merge set of RNA-Seq samples' chromosome gene coverage count files into one pandas.DataFrame with
    one row per gene, columns are `chr`, `gene`, <sample IDs> by scanning data_dir for sample ID subdirectories
    and extracting chromosome read count .csv files.

    See reads.BamReadsProcessor.coverage_read_counts method.

    Example:

    Suppose data_dir is comprised of a file tree structure like this:
    |-- data_dir
    |   |-- sample123
    |   |   |-- read_counts_sample123_chr1.csv
    |   |   |-- read_counts_sample123_chr2.csv
    |   |-- sample124
    |   |   |-- read_counts_sample124_chr1.csv
    |   |   |-- read_counts_sample124_chr2.csv

    ->> merge_read_counts(data_dir, ['sample123', 'sample124'], ['chr1', 'chr2']) ->>

    +-----------+---------+-----------------+-----------------+
    |    chr    |  gene   |    sample123    |    sample124    |
    +===========+=========+=================+=================+
    |   chr1    |  ATX2   |      275        |       101       |
    +-----------+---------+-----------------+-----------------+
    |    ...    |   ...   |      ...        |       ...       |
    +-----------+---------+-----------------+-----------------+
    |   chr2    |  GET4   |      301        |       255       |
    +-----------+---------+-----------------+-----------------+

    :param data_dir: str path of directory containing RNA-Seq sample ID subdirectories, one per sample ID contained
    in sample_ids, each subdirectory containing one read count .csv file per chromosome, named in the fashion
    "read_counts_<sample ID>_<chromosome>.csv"
    :param sample_ids: list of str names RNA Seq samples, i.e. basenames of various alignment files.
    :param chroms: list of str names of chromosomes for which to load read counts
    :return: pandas.DataFrame containing gene read counts across samples. Columns are `chr` (chromosome), `gene`,
    <sample IDs>
    """
    chrom_df_list = list()

    for chrom in chroms:
        for i in range(len(sample_ids)):

            # identify one (chromosome, sample ID) combination, and therefore, path to
            # of read counts .csv containing this chromosome's gene read counts for this sample.
            sample_id = sample_ids[i]
            counts_file = os.path.join(data_dir
                                       , sample_id
                                       , 'read_counts_{0}_{1}.csv'.format(sample_id, chrom))

            if not os.path.isfile(counts_file):
                raise IOError('read counts file {0} not available!'.format(counts_file))

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


def merge_overlap_gene_coverage(data_dir, sample_ids, chrom):
    """
    For single chromosome, join multiple RNA Seq alignment files' gene coverage vectors for genes in genome
    that overlap others on the chromosome of interest. Similar in spirit to merge_chrom_coverage,
    but join is over individual genes' coverage vectors, not entire chromosomes' coverage vectors.

    See reads.BamReadsProcessor.coverage_read_counts method.

    Example:

    Suppose data_dir is comprised of a file tree structure like this:
    |-- data_dir
    |   |-- sample123
    |   |   |-- overlap_coverage_sample123_chr1.pkl
    |   |   |-- overlap_coverage_sample123_chr2.pkl
    |   |-- sample124
    |   |   |-- overlap_coverage_sample124_chr1.pkl
    |   |   |-- overlap_coverage_sample124_chr2.pkl

    ->> merge_overlap_gene_coverage(data_dir, ['sample123', 'sample124'], 'chrj') ->>

    {('gene Aj'): <L1 x 2 coverage array>,
     ('gene Bj'): <L2 x 2 coverage array>,
     ...
     ('gene Nj'): <LNj x 2 coverage array>}

    :param data_dir: str path of directory containing RNA-Seq sample ID subdirectories, one per sample ID contained
    in sample_ids, each subdirectory containing one read count .csv file per chromosome, named in the fashion
    "overlap_coverage_<sample ID>_<chromosome>.csv"
    :param sample_ids: list of str names RNA Seq samples, i.e. basenames of various alignment files.
    :param chrom: str name of chromosome
    :return: dictionary of the form {gene_name: coverage numpy array} for genes in genome that overlap others
    on the chromosome of interest.
    """
    # output storage: (gene name, coverage matrix) key-value pairs.
    sample_cov_dict = dict()
    gene_cov_dict = dict()
    n_samples = len(sample_ids)

    # sample by sample, build gene coverage matrices.
    for i in range(n_samples):

        # identify one (chromosome, sample ID) combination, and therefore, path to
        # of read counts .csv containing this chromosome's gene read counts for this sample.
        sample_id = sample_ids[i]
        cov_file = os.path.join(data_dir
                                , sample_id
                                , 'overlap_coverage_{0}_{1}.pkl'.format(sample_id, chrom))

        # if there are no overlapping genes for this chromosome, return empty iterable.
        if not os.path.isfile(cov_file):
            return dict()

        # load (sample, chromosome) gene coverage vector dictionary.
        with open(cov_file, 'rb') as f:
            sample_cov_dict = pkl.load(f)

        for gene in sample_cov_dict:
            cov_vec = sample_cov_dict[gene]

            # if loading the first sample's coverage vectors, initialize gene coverage matrices.
            if i == 0:
                gene_cov_dict[gene] = np.zeros(shape=[n_samples, len(cov_vec)]
                                               , dtype=np.float_)

            # update sample's coverage within coverage matrix.
            gene_cov_dict[gene][i, :] = cov_vec

    del sample_cov_dict
    gc.collect()

    return gene_cov_dict


def merge_chrom_coverage(data_dir, sample_ids,
                         chrom_exon_df, verbose=True):
    """
    Join multiple RNA Seq alignment files' chromosome coverage vectors into a dictionary of per-gene
    overage matrices based on exon positioning for that chromosome.

    Example:

    Suppose data_dir is comprised of a file tree structure like this:

    |-- data_dir
    |   |-- sample123
    |   |   |-- chrom_coverage_sample123_chr1.pkl
    |   |-- sample124
    |   |   |-- chrom_coverage_sample124_chr1.pkl

    and suppose chrom_exon_df is a pandas.DataFrame looking like this:

    +-----------+----------+-----------------+-----------------+
    |    chr    |   gene   |   gene_start    |     gene_end    |
    +===========+==========+=================+=================+
    |   chr1    |  ZNF326  |      12275      |      14533      |
    +-----------+----------+-----------------+-----------------+
    |    ...    |    ...   |      ...        |       ...       |
    +-----------+----------+-----------------+-----------------+
    |   chr1    |   GORAB  |     3098838     |     3121055     |
    +-----------+----------+-----------------+-----------------+

    ->> merge_overlap_gene_coverage(data_dir, ['sample123', 'sample124'], chrom_exon_df) ->>

    {('gene Aj'): <L1 x 2 coverage array>,
     ('gene Bj'): <L2 x 2 coverage array>,
     ...
     ('gene Nj'): <LN x 2 coverage array>}

    :param data_dir: str path of directory containing RNA-Seq sample ID subdirectories,
     each subdirectory containing the chromosome of interest's coverage array in an .npz file, named in the fashion
    "chrom_coverage_<sample ID>_<chromosome>.csv"
    :param sample_ids: list of str names RNA Seq samples, i.e. basenames of various alignment files.
    :param chrom_exon_df: pandas.DataFrame outlining exon positions within a single chromosome; has columns 'chr',
    'start' (exon start), 'end' (exon end), 'gene' (gene name), 'gene_end', and 'gene_start'
    :param verbose: bool indicator should progress be written with logger?
    :return: dictionary of the form {gene_name: coverage numpy array} for isolated genes within specified chromosome
    """
    # output storage: (gene name, coverage matrix) key-value pairs.
    gene_cov_dict = dict()

    # identify the chromosome in question.
    unique_chrom = chrom_exon_df.chr.unique()
    if len(unique_chrom) > 1:
        raise ValueError('chrom_exon_df contains exon data for more than one chromosome!')

    chrom = unique_chrom[0]

    # identify all sample coverage arrays for this chromosome.
    npz_files = [os.path.join(data_dir, x, 'chrom_coverage_{0}_{1}.npz'.format(x, chrom)) for x in sample_ids]

    # determine if this chromosome even has chromosome coverage from isolated genes to parse.
    # This should never not be the case, but check for testing purposes.
    random_cov_file = np.random.choice(npz_files
                                       , size=1)[0]

    # if chromosome coverage file is missing for one sample, it will be missing for all samples,
    # so skip this chromosome because there will be no chromosome coverage vectors to dice.
    if not os.path.exists(random_cov_file):
        if verbose:
            logging.info('CHR {0}: no chromosome coverage files available.'.format(chrom))

        # return an *empty* gene coverage matrix dictionary.
        return dict()

    # Keep memory manageable:
    # break genes into groups so that each group's total coverage matrix
    # is ~ 500Mb. mem_splits dictates size of gene groups for load procedure.
    # use a randomly sampled chromosome coverage file to determine size of gene groups.
    cov_vec_sp = sparse.load_npz(random_cov_file)
    mem_splits = int(np.ceil(len(sample_ids) * cov_vec_sp.asfptype().todense().nbytes / 500e6))
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
        logging.info('CHR {0}: begin coverage matrix processing. \n'
                     'Using {1} gene splits for memory efficiency.'.format(chrom, mem_splits))

        # Instantiate progress bar if parsing non-negligible number of genes. Update in intervals of 5%.
        use_pbar = n_genes > 100
        if use_pbar:
            gene_idx = 0
            pbar_step_size = int(np.ceil(n_genes / 10))
            pbar = tqdm.tqdm(total=100
                             , leave=False
                             , desc='CHR {0}: coverage matrix progress'.format(chrom)
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
        for npz_file in npz_files:

            # load the gene span of a sample's compressed sparse row chromosome coverage array
            cov_vec_sp = sparse.load_npz(npz_file).transpose()[start_pos:end_pos, :]

            # initialize coverage matrix (len(chrom) x p) with copy of first experiment's coverage array.
            if idx == 0:
                cov_mat = cov_vec_sp.copy()

            # column-append sparse coverage vector to existing coverage matrix.
            else:
                cov_mat = sparse.hstack([cov_mat, cov_vec_sp]
                                        , dtype=int)

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
            # Coverage vectors are 0-indexed so take off 1 from 1-indexed gene positions,
            # but remembering to include exon end positions.
            e_starts, e_ends = single_gene_df.start.values - start_pos - 1, single_gene_df.end.values - start_pos
            slicing = [np.arange(e_starts[j], e_ends[j]) for j in range(len(e_starts))]

            # in case exons are overlapping, take union of their covered regions.
            slicing = np.unique(flatten_2d(slicing))

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
        logging.info('CHR {0} -- obtained {1} coverage matrices.'
                     .format(chrom, len(gene_cov_dict)))

    return gene_cov_dict


def merge_coverage(data_dir, sample_ids, exon_df, n_jobs=max_cpu(),
                   output_dir=None, verbose=True):
    """
    For each chromosome, load the coverage arrays resulting from each alignment file, join them,
    and then slice the joined coverage array into per-gene coverage matrices. Run process in parallel over
    chromosomes with the n_jobs argument.

    Example output:

    {('gene A'): <L1 x 2 coverage array>,
     ('gene B'): <L2 x 2 coverage array>,
     ...
     ('gene N'): <LN x 2 coverage array>}

    :param data_dir: str directory containing subdirectories named after the alignment sample IDs in sample_ids
    list, each containing files named `overlap_coverage_<sample ID>_<chromosome>.pkl` and
    `chrom_coverage_<sample ID>_<chromosome>.npz`.
    :param sample_ids: list of str names RNA Seq samples, i.e. basenames of various alignment files.
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

    # get list of gene coverage matrix dictionaries from joining chromosome-wide coverage arrays.
    chrom_gene_cov_dicts = Parallel(n_jobs=min(n_jobs, len(chroms))
                              , verbose=0
                              , backend='threading')(delayed(merge_chrom_coverage)(
        data_dir=data_dir,
        sample_ids=sample_ids,
        chrom_exon_df=subset_to_chrom(exon_df, chrom=chrom),
        verbose=verbose) for chrom in chroms)

    # get list of gene coverage matrix dictionaries upon joining gene coverage arrays for genes in overlapping groups.
    if verbose:
        logging.info('Joining overlapping genes\' coverage vectors into coverage matrices.')

    overlap_gene_cov_dicts = Parallel(n_jobs=min(n_jobs, len(chroms))
                              , verbose=0
                              , backend='threading')(delayed(merge_overlap_gene_coverage)(
        data_dir=data_dir,
        sample_ids=sample_ids,
        chrom=chrom) for chrom in chroms)

    # (1) concatenate chromosome coverage dictionaries into one ordered dictionary.
    # (2) for each chromosome, save gene coverage matrices to .pkl file in per-chromosome subdirectories.
    for i in range(len(chroms)):

        chrom = chroms[i]

        # concatenate coverage matrices from different sources (isolated or gene overlap genes).
        chrom_cov_dict = {**chrom_gene_cov_dicts[i], **overlap_gene_cov_dicts[i]}

        # move chromosome's gene coverage matrices into overall gene coverage matrix storage.
        for gene in chrom_cov_dict:
            gene_cov_dict[gene] = chrom_cov_dict[gene]

        # save {gene: coverage matrix} dictionary data per chromosome,
        # in a new directory named after the chromosome.
        if output_dir:
            save_dir = os.path.join(output_dir, chrom)

            if not os.path.isdir(save_dir):
                os.makedirs(save_dir)

            # save per-gene coverage matrices to .pkl file, one .pkl file per chromosome.
            chrom_cov_file = os.path.join(save_dir, 'coverage_matrices_{0}.pkl'.format(chrom))
            if verbose:
                logging.info('CHR {0} -- saving coverage matrices to {1}'
                             .format(chrom, chrom_cov_file))

            with open(chrom_cov_file, 'wb') as f:
                pkl.dump(chrom_cov_dict, f)

    del chrom_gene_cov_dicts, overlap_gene_cov_dicts, chrom_cov_dict
    gc.collect()

    return gene_cov_dict

