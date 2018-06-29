import os
import sys
import logging

# TODO write setup.py for install
sys.path.append('/Users/fineiskid/nu/jiping_research/degnorm')

from degnorm.coverage import *
from degnorm.gene_groupby import *
from degnorm.genome_preprocessing import *
from degnorm.utils import *


if __name__ == '__main__':

    logging.basicConfig(stream=sys.stdout
                        , level=logging.DEBUG
                        , format='%(asctime)s ---- %(message)s'
                        , datefmt='%m/%d/%Y %I:%M:%S')

    logging.info('Loading data')
    # LOAD
    # TODO: write pybedtools data converters, .gtf file manipulators
    # TODO: write command line argument parsers
    # TODO: write loaders
    data_dir = '~/nu/jiping_research/data/rna_seq'

    r1_file = os.path.join(data_dir, 'A01_R01.coverage')
    r2_file = os.path.join(data_dir, 'A01_R02.coverage')
    gene_file = os.path.join(data_dir, 'genes_exon.bed')

    r1_df = pd.read_table(r1_file, header=None)
    r1_df.columns = ['chr', 'start', 'end', 'count']

    r2_df = pd.read_table(r2_file, header=None)
    r2_df.columns = ['chr', 'start', 'end', 'count']

    gene_df = pd.read_table(gene_file, sep='\s+', header=None)
    gene_df.columns = ['chr', 'start', 'end', 'gene']

    # remove genes present in multiple chromosomes
    logging.info('Removing genes that appear in multiple chromosomes')
    gene_df = remove_multichrom_genes(gene_df)

    # create map of genes within genome
    logging.info('Establishing gene, exon positioning within genome')
    gene_df = get_gene_outline(gene_df)
    exon_df = get_exon_outline(gene_df)


    # cut out exons that overlap with multiple genes.
    logging.info('Removing exons that overlap with multiple genes')
    exon_df = remove_multigene_exons(exon_df
                                     , g_df = gene_df[['chr', 'gene_start', 'gene_end']]
                                     , n_jobs=max_cpu())

    # make list of pd.DataFrames from RNA-seq experiments
    sample_dfs = [r1_df, r2_df]
    n_samples = len(sample_dfs)

    # identify sample-union of chromosomes
    logging.info('Determining union of chromosomes over samples...')
    chroms = set(sample_dfs[0]['chr'].unique())
    for i in range(1, n_samples):
        chroms = set(sample_dfs[i]['chr'].unique()).intersection(chroms)

    chroms = list(chroms)
    n_chroms = len(chroms)
    logging.info('{0} chromosomes total'.format(n_chroms))

    # TODO: make gene coverage matrix storage smarter/what Bin wants
    gene_cov_dict = dict()

    # outer iteration loop: chromosomes
    for chrom_idx in range(n_chroms):
        chrom = chroms[chrom_idx]

        logging.info('Determining coverage matrices for genes in chromosome {0} -- {1} / {2}'
                     .format(chrom, chrom_idx + 1, n_chroms))

        chrom_sample_df_dict = dict()

        # subset RNA-seq sample data to chromosome
        for sample_idx in range(n_samples):
            chrom_sample_df_dict[sample_idx] = subset_to_chrom(sample_dfs[sample_idx], chrom=chrom)

        # subset genes to chromosome, merge processed exon regions with gene-level metadata.
        exon_sub_df = subset_to_chrom(exon_df, chrom=chrom)
        gene_sub_df = exon_sub_df.merge(gene_df[['chr', 'gene', 'gene_end', 'gene_start']]
                                        , on = ['chr', 'gene'])

        genes = gene_sub_df['gene'].unique()
        n_genes = len(genes)

        # TODO: parallelize gene loop (specific to a particular chromosome)
        for gene_idx in range(len(genes)):
            gene = genes[gene_idx]

            if gene_idx % 100 == 0 and gene_idx > 0:
                logging.info('On gene {0} -- {1} / {2}'.format(gene, gene_idx, n_genes))

            # subset chromosome's gene/exon data to a single gene.
            single_gene_df = gene_sub_df[gene_sub_df.gene == gene]

            rng = (single_gene_df['gene_start'].iloc[0], single_gene_df['gene_end'].iloc[0])
            cov_mat = np.zeros([rng[1] - rng[0], n_samples])

            for sample_idx in range(n_samples):
                cov_mat[:, sample_idx] = relative_sample_coverage(chrom_sample_df_dict[sample_idx]
                                                                  , rng=rng)

            # Slice up cov_mat based on relative exon positions within a gene.
            e_starts, e_ends = single_gene_df['exon_start'].values, single_gene_df['exon_end'].values
            slices = [np.arange(e_starts[i], e_ends[i]) - rng[0] for i in range(len(e_starts))]
            slicing = np.unique(np.array([sl for subslice in slices for sl in subslice]))

            gene_cov_dict[gene] = cov_mat[slicing, :]

            if gene_idx % 300 == 0 and gene_idx > 0:
                logging.info('Coverage matrix shape: {0}'.format(cov_mat.shape))
                logging.info('Coverage matrix mean coverage by sample: {0}'.format(cov_mat.mean(axis=0)))