from degnorm.reads import *
from degnorm.coverage_counts import *
from degnorm.gene_processing import *
from degnorm.utils import *
from degnorm.nmf import *
from datetime import datetime
from collections import OrderedDict
import time


def main():

    # ---------------------------------------------------------------------------- #
    # Load CLI arguments and initialize storage, output directory
    # ---------------------------------------------------------------------------- #
    args = parse_args()
    n_jobs = args.cpu
    n_samples = len(args.input_files)
    sample_ids = list()
    chroms = list()
    cov_files = dict()
    reads_dict = dict()

    output_dir = os.path.join(args.output_dir, 'DegNorm_' + datetime.now().strftime('%m%d%Y_%H%M%S'))
    if os.path.isdir(output_dir):
        time.sleep(2)
        output_dir = os.path.join(args.output_dir, 'DegNorm_' + datetime.now().strftime('%m%d%Y_%H%M%S'))

    os.makedirs(output_dir)

    # ---------------------------------------------------------------------------- #
    # Load .sam or .bam files; parse into coverage arrays and store them.
    # ---------------------------------------------------------------------------- #

    # iterate over .sam files; compute each sample's chromosomes' coverage arrays
    # and save them to .npz files.
    for idx in range(n_samples):
        logging.info('Loading RNA-seq data file {0} / {1}'.format(idx + 1, n_samples))
        sam_file = args.input_files[idx]

        # if actually working with .bam files, convert them to .sam.
        if sam_file.endswith('.bam'):
            logging.info('Converting {0} into .sam file format...'
                         .format(sam_file))
            sam_file = bam_to_sam(sam_file)

        reader = ReadsCoverageProcessor(sam_file=sam_file
                                        , n_jobs=args.cpu
                                        , tmp_dir=output_dir
                                        , verbose=True)
        sample_id = reader.sample_id
        sample_ids.append(sample_id)
        cov_files[sample_id] = reader.coverage()
        reads_dict[sample_id] = reader.data

        chroms += [os.path.basename(f).split('.npz')[0].split('_')[-1] for f in cov_files[sample_id]]

    chroms = list(set(chroms))
    logging.info('Obtained reads coverage for {0} chromosomes:\n'
                 '{1}'.format(len(chroms), ', '.join(chroms)))

    # ---------------------------------------------------------------------------- #
    # Load .gtf or .gff files and run processing pipeline.
    # Run in parallel over chromosomes.
    # ---------------------------------------------------------------------------- #
    logging.info('Begin genome annotation file processing...')
    gap = GeneAnnotationProcessor(args.genome_annotation
                                  , n_jobs=n_jobs
                                  , verbose=True
                                  , chroms=chroms
                                  , genes=args.genes)
    exon_df = gap.run()
    genes_df = exon_df[['chr', 'gene', 'gene_start', 'gene_end']].drop_duplicates().reset_index(drop=True)

    # ---------------------------------------------------------------------------- #
    # Obtain read count matrix: X, an n (genes) x p (samples) matrix
    # Run in parallel over samples.
    # ---------------------------------------------------------------------------- #
    logging.info('Obtaining read count matrix...')

    p = mp.Pool(processes=n_jobs)
    read_count_vecs = [p.apply_async(read_counts
                                     , args=(reads_dict[sample_id], genes_df)) for sample_id in sample_ids]
    p.close()
    read_count_vecs = [x.get() for x in read_count_vecs]

    if n_samples > 1:
        X = np.vstack(read_count_vecs).T
    else:
        X = np.reshape(read_count_vecs[0], (-1, 1))

    logging.info('Successfully obtained read count matrix -- shape: {0}'.format(X.shape))

    # free up some memory: delete reads data
    del reads_dict, read_count_vecs, reader
    
    # ---------------------------------------------------------------------------- #
    # Slice up genome coverage matrix for each gene according to exon positioning.
    # Run in parallel over chromosomes.
    # ---------------------------------------------------------------------------- #
    cov_output = None if args.disregard_coverage else output_dir
    p = mp.Pool(processes=n_jobs)
    gene_cov_mats = [p.apply_async(gene_coverage
                                   , args=(exon_df, chrom, cov_files, cov_output, True)) for chrom in chroms]
    p.close()
    gene_cov_mats = [x.get() for x in gene_cov_mats]

    # ---------------------------------------------------------------------------- #
    # Process gene coverage matrix output prior to running NMF.
    # ---------------------------------------------------------------------------- #

    # convert list of tuples into 2-d dictionary: {chrom: {gene: coverage matrix}}, and
    # initialized an OrderedDict to store
    chrom_gene_cov_dict = {gene_cov_mats[i][1]: gene_cov_mats[i][0] for i in range(len(gene_cov_mats))}
    gene_cov_dict = OrderedDict()

    # Determine for which genes to run DegNorm, and for genes where we will run DegNorm,
    # which transcript regions to filter out prior to running DegNorm.
    logging.info('Determining genes include in DegNorm coverage curve approximation.')
    delete_idx = list()
    for i in range(genes_df.shape[0]):
        chrom = genes_df.chr.iloc[i]
        gene = genes_df.gene.iloc[i]
        cov_mat = chrom_gene_cov_dict[chrom][gene]

        if any(cov_mat.sum(axis=0) == 0):
            delete_idx.append(i)

        else:
            gene_cov_dict[gene] = cov_mat.T

    if delete_idx:
        X = np.delete(X
                      , obj=delete_idx
                      , axis=0)
        genes_df = genes_df.drop(delete_idx
                                 , axis=0).reset_index(drop=True)

    logging.info('DegNorm will run on {0} genes.'.format(len(gene_cov_dict)))

    # check that read counts and coverage matrices contain data for same number of genes.
    if len(gene_cov_dict.keys()) != X.shape[0]:
        raise ValueError('Number of coverage matrices not equal to number of genes in read count matrix!')

    # save gene annotation metadata.
    gene_output_file = os.path.join(output_dir, 'gene_metadata.csv')
    logging.info('Saving gene metadata.')
    genes_df.to_csv(gene_output_file
                    , index=False)

    # save read counts.
    logging.info('Saving read counts.')
    np.savetxt(os.path.join(output_dir, 'read_counts.csv')
               , X=X
               , delimiter=',')

    # free up more memory: delete exon / genome annotation data
    del chrom_gene_cov_dict, gene_cov_mats

    # ---------------------------------------------------------------------------- #
    # Run NMF.
    # ---------------------------------------------------------------------------- #

    logging.info('Executing NMF over-approximation algorithm...')
    nmfoa = GeneNMFOA(nmf_iter=20, grid_points=2000, n_jobs=n_jobs)
    nmfoa.fit_transform(gene_cov_dict, reads_dat=X)

    # extract approximated coverage curves.
    estimates = [nmfoa.cov_dat[gene]['estimated_coverage'] for gene in nmfoa.baseline_genes]

    # ---------------------------------------------------------------------------- #
    # Save
    # ---------------------------------------------------------------------------- #

    # adjust read counts based on DI scores computed from final coverage estimates.
    logging.info('Computing degradation index scores, adjusting read counts.')
    nmfoa.compute_scale_factors(estimates)
    X_adj = nmfoa.x / (1 - nmfoa.rho)

    logging.info('Saving degradation index scores, adjusted read counts.')
    np.savetxt(os.path.join(output_dir, 'degradation_index.csv')
               , X=nmfoa.rho
               , delimiter=',')

    np.savetxt(os.path.join(output_dir, 'adjusted_read_counts.csv')
               , X=X_adj
               , delimiter=',')

    # ---------------------------------------------------------------------------- #
    # Generate plots.
    # ---------------------------------------------------------------------------- #





if __name__ == "__main__":
    main()