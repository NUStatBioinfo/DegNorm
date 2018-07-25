from degnorm.utils import *
import re
import pickle as pkl


def gene_coverage(exon_df, chrom, coverage_files, output_dir=None):
    """
    Slice up a coverage matrix for a chromosome into a dictionary of per-gene
    coverage matrices based on exon positioning for that chromosome.

    :param exon_df: pandas.DataFrame outlining exon positions in a gene; has columns 'chr',
    'start' (exon start), 'end' (exon end), 'gene' (gene name), 'gene_end', and 'gene_start'
    :param chrom: str name of chromosome to break up for gene-level coverage matrices
    :param coverage_files: dict of {sample ID: list of .npz files} specifying, per RNA-seq experiment, the paths to
    compressed numpy chromosome coverage array files.
    :param output_dir: str if specified save gene-level coverage matrices to binary .pkl files.
    Default is None (do not save)
    :return: dictionary with {gene_name: coverage numpy array} data
    """
    exon_chrom_df = subset_to_chrom(exon_df, chrom=chrom)

    if exon_chrom_df.empty:
        raise ValueError('Chromosome {0} not present in exon_df!'.format(chrom))

    genes = exon_chrom_df['gene'].unique()
    n_genes = len(genes)

    idx = 0
    for sample_id in coverage_files:
        r = re.compile('sample_{0}_{1}.npz'.format(sample_id, chrom))
        npz_file = list(filter(r.search, coverage_files[sample_id]))[0]
        cov_vec = np.load(npz_file)['cov']

        if idx == 0:
            cov_mat = np.zeros([len(cov_vec), len(coverage_files)])

        cov_mat[:, idx] = cov_vec
        idx += 1

    logging.info('CHROMOSOME {0} -- coverage matrix shape: {1}'.format(chrom, cov_mat.shape))

    # store coverage matrices in a dictionary with gene name keys
    gene_cov_dict = dict()

    for gene_idx in range(n_genes):
        gene = genes[gene_idx]

        if gene_idx % 500 == 0 and gene_idx > 0:
            logging.info('GENE {0} -- {1} / {2}'.format(gene, gene_idx, n_genes))

        # subset chromosome's gene/exon data to a single gene and
        # identify gene's start and end positions.
        single_gene_df = exon_chrom_df[exon_chrom_df.gene == gene]

        # Slice up cov_mat based on relative exon positions within a gene.
        e_starts, e_ends = single_gene_df.start.values, single_gene_df.end.values
        slices = [np.arange(e_starts[i], e_ends[i]) for i in range(len(e_starts))]
        slicing = np.unique(flatten_2d(slices))

        gene_cov_dict[gene] = cov_mat[slicing, :]

        if gene_idx % 500 == 0 and gene_idx > 0:
            logging.info('GENE {0} -- coverage matrix shape: {1}'.format(gene, gene_cov_dict[gene].shape))
            logging.info('GENE {0} -- mean coverage by sample: {1}'.format(gene, gene_cov_dict[gene].mean(axis=0)))

    if output_dir:
        output_file = os.path.join(output_dir, 'coverage_matrices_{0}.pkl'.format(chrom))

        with open(output_file, 'wb') as f:
            logging.info('Saving gene coverage matrices to binary pickle file -- {0}'.format(output_file))
            pkl.dump(gene_cov_dict, f)

    return gene_cov_dict, chrom