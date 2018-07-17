import os
import sys
from degnorm.reads_coverage import *
from degnorm.gene_coverage import *
from degnorm.gene_groupby import *
from degnorm.genome_preprocessing import *
from degnorm.utils import *


def main():

    # TODO: write .gtf file manipulators
    # TODO: write command line argument parsers
    data_dir = '/Users/fineiskid/nu/jiping_research/data/rna_seq'

    n_jobs = 2
    sam_files = [os.path.join(data_dir, 'A01_R01_chr1_2.sam')]
    n_samples = len(sam_files)
    samples = list(range(n_samples))

    # initialize storage for chromosome names and coverage array file paths.
    chroms = list()
    sample_cov_files = dict()

    # iterate over .sam files; compute each sample's chromosomes' coverage arrays
    # and save them to .npz files.
    logging.info('Begin loading RNA-seq data...')
    for idx in samples:
        sam_reader = ReadsCoverageParser(sam_file=sam_files[idx]
                                         , verbose=True
                                         , n_jobs=n_jobs)
        sample_cov_files[idx] = sam_reader.coverage()
        chroms += [os.path.basename(f).split('.npz')[0] for f in sample_cov_files[idx]]

    chroms = list(set(chroms))
    n_chroms= len(chroms)
    logging.info('Total chromosomes: {0}'.format(n_chroms))

    logging.info('Loading gene annotation exon file')
    gene_file = os.path.join(data_dir, 'genes_exon.bed')
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
                                     , g_df=gene_df[['chr', 'gene_start', 'gene_end']]
                                     , n_jobs=n_jobs)

    # outer iteration loop: chromosomes
    for chrom_idx in range(n_chroms):
        chrom = chroms[chrom_idx]

        logging.info('Loading total coverage matrix for chromosome {0} from .npz files. {1} / {2}'
                     .format(chrom, chrom_idx + 1, n_chroms))
        for sample_idx in samples:
            r = re.compile('{0}.npz'.format(chrom))
            npz_file = list(filter(r.search, sample_cov_files[sample_idx]))[0]
            cov_vec = np.load(npz_file)['cov']

            if sample_idx == 0:
                cov_mat = np.zeros([len(cov_vec), n_samples])

            cov_mat[:, sample_idx] = cov_vec

        logging.info('CHROMOSOME {0} -- coverage matrix shape: {1}'.format(chrom, cov_mat.shape))

        # subset genes to chromosome, merge processed exon regions with gene-level metadata.
        exon_sub_df = subset_to_chrom(exon_df, chrom=chrom)
        gene_sub_df = exon_sub_df.merge(gene_df[['chr', 'gene', 'gene_end', 'gene_start']]
                                        , on=['chr', 'gene'])

        genes = gene_sub_df['gene'].unique()
        n_genes = len(genes)

        # TODO: parallelize gene loop (specific to a particular chromosome)
        for gene_idx in range(len(genes)):
            gene = genes[gene_idx]

            if gene_idx % 100 == 0 and gene_idx > 0:
                logging.info('GENE {0} -- {1} / {2}'.format(gene, gene_idx, n_genes))

            # subset chromosome's gene/exon data to a single gene.
            # identify gene's start and end positions.
            single_gene_df = gene_sub_df[gene_sub_df.gene == gene]
            rng = (single_gene_df['gene_start'].iloc[0], single_gene_df['gene_end'].iloc[0])

            # initialize gene's coverage matrix.
            gene_cov_mat = np.zeros([rng[1] - rng[0], n_samples])

            # # iterate over samples, building gene's coverage matrix.
            # for sample_idx in samples:
            #     cov_array = rna_dat_dict[sample_idx][chrom][np.arange(rng[0], rng[1])]
            #     cov_mat[:, sample_idx] = cov_array

            # Slice up cov_mat based on relative exon positions within a gene.
            e_starts, e_ends = single_gene_df['exon_start'].values, single_gene_df['exon_end'].values
            # slices = [np.arange(e_starts[i], e_ends[i]) - rng[0] for i in range(len(e_starts))]
            slices = [np.arange(e_starts[i], e_ends[i]) for i in range(len(e_starts))]
            slicing = np.unique(flatten_2d(slices))

            gene_cov_dict[gene] = gene_cov_mat[slicing, :]

            if gene_idx % 300 == 0 and gene_idx > 0:
                logging.info('Coverage matrix shape: {0}'.format(gene_cov_dict[gene].shape))
                logging.info('Coverage matrix mean coverage by sample: {0}'.format(gene_cov_dict[gene].mean(axis=0)))