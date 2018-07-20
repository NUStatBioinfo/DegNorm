from degnorm.reads_coverage import *
from degnorm.gene_processing import *
from degnorm.gene_groupby import *
from degnorm.genome_preprocessing import *
from degnorm.utils import *
from datetime import datetime
import time


def main():

    # ---------------------------------------------------------------------------- #
    # Load CLI arguments and initialize storage, output directory
    # ---------------------------------------------------------------------------- #

    args = parse_args()
    n_jobs = args.cpu
    n_samples = len(args.input_files)
    samples = list(range(n_samples))
    chroms = list()
    sample_cov_files = dict()

    output_dir = os.path.join(args.output_dir, 'DegNorm_' + datetime.now().strftime('%m%d%Y_%H%M%S'))
    if os.path.isdir(output_dir):
        time.sleep(5)
        output_dir = os.path.join(args.output_dir, 'DegNorm_' + datetime.now().strftime('%m%d%Y_%H%M%S'))

    os.makedirs(output_dir)

    # ---------------------------------------------------------------------------- #
    # Load .sam or .bam files; parse into coverage arrays and store them.
    # ---------------------------------------------------------------------------- #

    # iterate over .sam files; compute each sample's chromosomes' coverage arrays
    # and save them to .npz files.
    for idx in samples:
        logging.info('Loading RNA-seq data file {0} / {1}'.format(idx + 1, n_samples))
        if args.input_type == 'sam':
            reader = ReadsCoverageProcessor(sam_file=args.input_files[idx]
                                            , n_jobs=args.cpu
                                            , tmp_dir=output_dir
                                            , verbose=True)

        sample_cov_files[idx] = reader.coverage()
        chroms += [os.path.basename(f).split('.npz')[0] for f in sample_cov_files[idx]]

    chroms = list(set(chroms))
    n_chroms= len(chroms)
    logging.info('Total chromosomes: {0}'.format(n_chroms))

    # ---------------------------------------------------------------------------- #
    # Load .gtf or .gff files and run processing pipeline.
    # ---------------------------------------------------------------------------- #

    # load and process gene annotation file
    logging.info('Begin genome annotation file processing...')
    gap = GeneAnnotationProcessor(args.genome_annotation
                                  , n_jobs=n_jobs
                                  , verbose=True)
    exon_df = gap.run()

    def gene_coverage(exon_df, chrom, coverage_files, samples, output_dir):

        for sample_idx in samples:
            r = re.compile('{0}.npz'.format(chrom))
            npz_file = list(filter(r.search, sample_cov_files[sample_idx]))[0]
            cov_vec = np.load(npz_file)['cov']

            if sample_idx == 0:
                cov_mat = np.zeros([len(cov_vec), n_samples])

            cov_mat[:, sample_idx] = cov_vec

        logging.info('CHROMOSOME {0} -- coverage matrix shape: {1}'.format(chrom, cov_mat.shape))
        exon_chrom_df = subset_to_chrom(exon_df, chrom=chrom)
        genes = exon_chrom_df['gene'].unique()
        n_genes = len(genes)

        # store coverage matrices in a dictionary with gene name keys
        gene_cov_dict = dict()

        for gene_idx in range(len(genes)):
            gene = genes[gene_idx]

            if gene_idx % 100 == 0 and gene_idx > 0:
                logging.info('GENE {0} -- {1} / {2}'.format(gene, gene_idx, n_genes))

            # subset chromosome's gene/exon data to a single gene and
            # identify gene's start and end positions.
            single_gene_df = exon_chrom_df[exon_chrom_df.gene == gene]
            # rng = (single_gene_df['gene_start'].iloc[0], single_gene_df['gene_end'].iloc[0])

            # Slice up cov_mat based on relative exon positions within a gene.
            e_starts, e_ends = single_gene_df.start.values, single_gene_df.end.values
            # slices = [np.arange(e_starts[i], e_ends[i]) - rng[0] for i in range(len(e_starts))]
            slices = [np.arange(e_starts[i], e_ends[i]) for i in range(len(e_starts))]
            slicing = np.unique(flatten_2d(slices))

            gene_cov_dict[gene] = cov_mat[slicing, :]

            if gene_idx % 300 == 0 and gene_idx > 0:
                logging.info('GENE {0} -- coverage matrix shape: {1}'.format(gene, gene_cov_dict[gene].shape))
                logging.info('GENE {0} -- mean coverage by sample: {1}'.format(gene, gene_cov_dict[gene].mean(axis=0)))


    # subset genes to chromosome, merge processed exon regions with gene-level metadata.
    # outer iteration loop: chromosomes
    # for chrom_idx in range(n_chroms):
    #     chrom = chroms[chrom_idx]
    #
    #     logging.info('Loading total coverage matrix for chromosome {0} from .npz files. {1} / {2}'
    #                  .format(chrom, chrom_idx + 1, n_chroms))
    #     for sample_idx in samples:
    #         r = re.compile('{0}.npz'.format(chrom))
    #         npz_file = list(filter(r.search, sample_cov_files[sample_idx]))[0]
    #         cov_vec = np.load(npz_file)['cov']
    #
    #         if sample_idx == 0:
    #             cov_mat = np.zeros([len(cov_vec), n_samples])
    #
    #         cov_mat[:, sample_idx] = cov_vec

        # logging.info('CHROMOSOME {0} -- coverage matrix shape: {1}'.format(chrom, cov_mat.shape))

        # exon_chrom_df = subset_to_chrom(exon_df, chrom=chrom)
        # genes = exon_chrom_df['gene'].unique()
        # n_genes = len(genes)

        # TODO: parallelize gene loop (specific to a particular chromosome)
        # for gene_idx in range(len(genes)):
        #     gene = genes[gene_idx]

            # if gene_idx % 100 == 0 and gene_idx > 0:
            #     logging.info('GENE {0} -- {1} / {2}'.format(gene, gene_idx, n_genes))
            #
            # # subset chromosome's gene/exon data to a single gene.
            # # identify gene's start and end positions.
            # single_gene_df = exon_chrom_df[exon_chrom_df.gene == gene]
            # rng = (single_gene_df['gene_start'].iloc[0], single_gene_df['gene_end'].iloc[0])
            #
            # # initialize gene's coverage matrix.
            # gene_cov_mat = np.zeros([rng[1] - rng[0], n_samples])
            #
            # # # iterate over samples, building gene's coverage matrix.
            # # for sample_idx in samples:
            # #     cov_array = rna_dat_dict[sample_idx][chrom][np.arange(rng[0], rng[1])]
            # #     cov_mat[:, sample_idx] = cov_array
            #
            # # Slice up cov_mat based on relative exon positions within a gene.
            # e_starts, e_ends = single_gene_df['exon_start'].values, single_gene_df['exon_end'].values
            # # slices = [np.arange(e_starts[i], e_ends[i]) - rng[0] for i in range(len(e_starts))]
            # slices = [np.arange(e_starts[i], e_ends[i]) for i in range(len(e_starts))]
            # slicing = np.unique(flatten_2d(slices))
            #
            # gene_cov_dict[gene] = gene_cov_mat[slicing, :]
            #
            # if gene_idx % 300 == 0 and gene_idx > 0:
            #     logging.info('Coverage matrix shape: {0}'.format(gene_cov_dict[gene].shape))
            #     logging.info('Coverage matrix mean coverage by sample: {0}'.format(gene_cov_dict[gene].mean(axis=0)))