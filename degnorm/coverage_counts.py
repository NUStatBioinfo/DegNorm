from degnorm.utils import *
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
    :param coverage_files: dict of {sample ID: list of .npz files} specifying, per RNA-seq experiment, the paths to
    compressed numpy chromosome coverage array files.
    :param output_dir: str if specified save gene-level coverage matrices to binary .pkl files.
    :param verbose: bool indicator should progress be written with logger?
    Default is None (do not save)
    :return: dictionary with {gene_name: coverage numpy array} data
    """
    exon_chrom_df = subset_to_chrom(exon_df, chrom=chrom)

    if exon_chrom_df.empty:
        raise ValueError('Chromosome {0} not present in exon_df!'.format(chrom))

    genes = exon_chrom_df['gene'].unique()
    n_genes = len(genes)

    if verbose:
        logging.info('CHROMOSOME {0} -- begin loading coverage matrix.'.format(chrom))

    idx = 0
    for sample_id in coverage_files:
        r = re.compile('sample_{0}_{1}.npz'.format(sample_id, chrom))
        npz_file = list(filter(r.search, coverage_files[sample_id]))[0]
        cov_vec = np.load(npz_file)['cov']

        if idx == 0:
            cov_mat = np.zeros([len(cov_vec), len(coverage_files)])

        cov_mat[:, idx] = cov_vec
        idx += 1

    if verbose:
        logging.info('CHROMOSOME {0} -- coverage matrix shape: {1}'.format(chrom, cov_mat.shape))

    # store coverage matrices in a dictionary with gene name keys
    gene_cov_dict = dict()

    # Instantiate progress bar.
    pbar = tqdm.tqdm(total=n_genes
                     , leave=False
                     , desc='CHROMOSOME {0} -- coverage matrix progress'.format(chrom))

    for gene_idx in range(n_genes):
        gene = genes[gene_idx]

        # subset chromosome's gene/exon data to a single gene and
        # identify gene's start and end positions.
        single_gene_df = exon_chrom_df[exon_chrom_df.gene == gene]

        # Slice up cov_mat based on relative exon positions within a gene.
        e_starts, e_ends = single_gene_df.start.values, single_gene_df.end.values
        slices = [np.arange(e_starts[i], e_ends[i]) for i in range(len(e_starts))]
        slicing = np.unique(flatten_2d(slices))

        # Save transposed coverage matrix so that shape is p x Li.
        gene_cov_dict[gene] = cov_mat[slicing, :].T

        pbar.update()

    pbar.close()

    # if a target location is specified, save {gene: coverage matrix} data per chromosome in a
    # new directory named after the chromosome.
    if output_dir:
        output_dir = os.path.join(output_dir, chrom)

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        # save per-gene coverage matrices to .pkl files
        gene_cov_output_file = os.path.join(output_dir, 'coverage_matrices_{0}.pkl'.format(chrom))
        if verbose:
            logging.info('Saving {0} coverage matrices to {1}'.format(chrom, gene_cov_output_file))

        with open(gene_cov_output_file, 'wb') as f:
            pkl.dump(gene_cov_dict, f)

    return gene_cov_dict, chrom


def read_counts(reads_df, genes_df):
    """
    Determine the read count for gene i in sample j, defined as the number of paired reads
    that land within the exonic region of the ith gene (in sample j).

    :param reads_df: pandas.DataFrame with columns "chr," "pos", and "qname_unpaired" designating the
    query name of a paired read, but without the pairing index
    :param genes_df: pandas.DataFrame with columns "chr", "gene", "gene_start", and "gene_end" designating
    the start and end position of a gene on a chromosome.
    :return: 1-d numpy array of read counts, length is equal to number of genes
    """
    counts = np.zeros(genes_df.shape[0])
    chroms = genes_df.chr.unique().tolist()
    chroms.sort()

    count_idx = 0
    for chrom in chroms:
        reads_sub_df = subset_to_chrom(reads_df, chrom=chrom)
        genes_sub_df = subset_to_chrom(genes_df, chrom=chrom)
        dat = genes_sub_df[['gene_start', 'gene_end']].values

        for i in range(genes_sub_df.shape[0]):
            counts_df = reads_sub_df[reads_sub_df.pos.between(dat[i, 0], dat[i, 1])]
            counts[count_idx + i] = counts_df.shape[0] / 2

        count_idx += i

    return counts