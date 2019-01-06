from degnorm.utils import *
from scipy import sparse
import pickle as pkl
import re
import tqdm


def gene_coverage(exon_df, chrom, coverage_files, output_dir=None, verbose=True):
    """
    Slice up a coverage matrix for a chromosome into a dictionary of per-gene
    coverage matrices based on exon positioning for that chromosome.

    :param exon_df: pandas.DataFrame outlining exon positions in a gene; has columns 'chr',
    'start' (exon start), 'end' (exon end), 'gene' (gene name), 'gene_end', and 'gene_start'
    :param chrom: str name of chromosome to break up for gene-level coverage matrices
    :param coverage_files: OrderedDict of {sample ID: list of .npz files} specifying, per RNA-seq experiment,
     the paths to compressed numpy chromosome coverage array files.
    :param output_dir: str if specified save gene-level coverage matrices to binary .pkl files.
    :param verbose: bool indicator should progress be written with logger?
    Default is None (do not save)
    :return: dictionary with {gene_name: coverage numpy array} data
    """
    exon_chrom_df = subset_to_chrom(exon_df, chrom=chrom)

    if exon_chrom_df.empty:
        raise ValueError('Chromosome {0} not present in exon_df!'.format(chrom))

    # Keep memory manageable:
    # attempt to break genes into groups so that each group's total coverage matrix
    # is ~ 500Mb. mem_splits dictates size of gene groups for load procedure.
    random_sample_id = np.random.choice(list(coverage_files.keys())
                                        , size=1)[0]
    r = re.compile('sample_{0}_{1}.npz'.format(random_sample_id, chrom))
    npz_file = list(filter(r.search, coverage_files[random_sample_id]))[0]
    cov_vec_sp = sparse.load_npz(npz_file)
    mem_splits = int(np.ceil(len(coverage_files) * cov_vec_sp.asfptype().todense().nbytes / 500e6))
    del cov_vec_sp

    # sort genes by end position so we won't need to have entire chromosome coverage vectors loaded at once,
    # then break genes up into mem_splits subsets.
    exon_chrom_df = exon_chrom_df.sort_values('gene_end'
                                              , axis=0)
    genes = exon_chrom_df['gene'].unique()
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

    # output storage: (gene name, coverage matrix) key-value pairs.
    gene_cov_dict = dict()

    use_pbar = False
    if verbose:
        logging.info('CHROMOSOME {0}: begin processing coverage matrices. \n'
                     'Using {1} gene splits for memory efficiency.'.format(chrom, mem_splits))

        # Instantiate progress bar if parsing non-negligible number of genes. Update in intervals of 5%.
        use_pbar = n_genes > 100
        if use_pbar:
            gene_idx = 0
            pbar_step_size = int(np.ceil(n_genes / 20))
            pbar = tqdm.tqdm(total=100
                             , leave=False
                             , desc='CHROMOSOME {0}: gene coverage matrix progress'.format(chrom)
                             , unit='%')

    # create the coverage matrix for each subset of genes.
    for i in range(mem_splits):

        # subset exon data to current gene subset.
        sub_genes = genes[gene_splits[i]:gene_splits[i + 1]].tolist()
        sub_exon_chrom_df = exon_chrom_df[exon_chrom_df.gene.isin(sub_genes)]

        # determine gene span: we only need a subset of the chromosome's coverage for gene subset.
        start_pos = int(sub_exon_chrom_df.gene_start.min())
        end_pos = int(sub_exon_chrom_df.gene_end.max() + 1)

        # load up gene span's coverage matrix.
        idx = 0
        for sample_id in coverage_files:
            r = re.compile('sample_{0}_{1}.npz'.format(sample_id, chrom))
            npz_file = list(filter(r.search, coverage_files[sample_id]))[0]

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
            single_gene_df = sub_exon_chrom_df[sub_exon_chrom_df.gene == gene]

            # Slice up cov_mat based on relative exon positions within a gene while remembering to
            # shift starts and ends based on the start position of the current gene span.
            e_starts, e_ends = single_gene_df.start.values - start_pos, single_gene_df.end.values - start_pos
            slices = [np.arange(e_starts[j], e_ends[j]) for j in range(len(e_starts))]

            # in case exons are overlapping, take union of their covered regions.
            slicing = np.unique(flatten_2d(slices))

            # Save transposed coverage matrix so that shape is p x Li.
            gene_cov_dict[gene] = np.array(cov_mat[slicing, :].T).astype(np.float_)

            if use_pbar:
                if (gene_idx % pbar_step_size == 0) and (gene_idx > 0):
                    pbar.update(5)
                    gene_idx += 1

    # close progress bar if using one.
    if use_pbar:
        pbar.close()

    # free up memory allocation for dense coverage matrix, exon data subset.
    del cov_mat, sub_exon_chrom_df
    gc.collect()

    # if a target location is specified, save {gene: coverage matrix} data per chromosome in a
    # new directory named after the chromosome.
    if output_dir:
        output_dir = os.path.join(output_dir, chrom)

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        # save per-gene coverage matrices to .pkl files
        gene_cov_output_file = os.path.join(output_dir, 'coverage_matrices_{0}.pkl'.format(chrom))
        if verbose:
            logging.info('CHROMOSOME {0}: successfully parsed gene coverage matrices. Saving to {1}'
                         .format(chrom, gene_cov_output_file))

        with open(gene_cov_output_file, 'wb') as f:
            pkl.dump(gene_cov_dict, f)

    return gene_cov_dict