from degnorm.reads import *
from degnorm.coverage import *
from degnorm.gene_processing import *
from degnorm.visualizations import *
from degnorm.nmf import *
from degnorm.report import render_report
from datetime import datetime
from collections import OrderedDict
import time
import sys


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
    read_count_dict = dict()

    output_dir = os.path.join(args.output_dir, 'DegNorm_' + datetime.now().strftime('%m%d%Y_%H%M%S'))
    if os.path.isdir(output_dir):
        time.sleep(2)
        output_dir = os.path.join(args.output_dir, 'DegNorm_' + datetime.now().strftime('%m%d%Y_%H%M%S'))

    os.makedirs(output_dir)

    # ---------------------------------------------------------------------------- #
    # Overhead: if supplied with .bam files, convert to .sam;
    # further, determine intersection of chromosomes across samples.
    # ---------------------------------------------------------------------------- #
    for idx in range(n_samples):
        sam_file = args.input_files[idx]

        # if file actually .bam, convert to .sam.
        if sam_file.endswith('.bam'):
            logging.info('Converting input file {0} into .sam file format...'
                         .format(sam_file))
            args.input_files[idx] = bam_to_sam(sam_file)

        # Take intersection of chromosomes over samples.
        if idx == 0:
            chroms = SamLoader(args.input_files[idx]).find_chromosomes()

        else:
            chroms = np.intersect1d(chroms, SamLoader(args.input_files[idx]).find_chromosomes()).tolist()

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

    # take intersection of chromosomes available in genome annotation file and those in the reads data,
    # if for some reason annotation file only contains subset.
    chroms = np.intersect1d(chroms, genes_df.chr.unique()).tolist()
    logging.info('Found {0} chromosomes in intersection of all experiments and gene annotation data:\n'
                 '\t{1}'.format(len(chroms), ', '.join(chroms)))

    # ---------------------------------------------------------------------------- #
    # Load .sam or .bam files, parse into coverage arrays, store them, and
    # compute gene read counts.
    # ---------------------------------------------------------------------------- #

    # iterate over .sam files; compute each sample's chromosomes' coverage arrays
    # and save them to .npz files.
    for idx in range(n_samples):
        logging.info('Loading RNA-seq data file {0} / {1}'.format(idx + 1, n_samples))
        sam_file = args.input_files[idx]

        reader = ReadsProcessor(sam_file=sam_file
                                , chroms=chroms
                                , n_jobs=n_jobs
                                , tmp_dir=output_dir
                                , verbose=True)

        sample_id = reader.sample_id
        sample_ids.append(sample_id)
        cov_files[sample_id], read_count_dict[sample_id] = reader.coverage_read_counts(genes_df)

    logging.info('Successfully processed chromosome read coverage and gene read counts for all {0} experiments'
                 .format(len(sample_ids)))

    del reader

    # ---------------------------------------------------------------------------- #
    # Merge per-sample gene read count matrices:
    # obtain read count DataFrame containing X, an n (genes) x p (samples) matrix.
    # ---------------------------------------------------------------------------- #

    read_count_df = read_count_dict[sample_ids[0]]
    read_count_df.rename(columns={'read_count': sample_ids[0]}, inplace=True)

    for sample_id in sample_ids[1:]:
        read_count_df = pd.merge(read_count_df
                                 , read_count_dict[sample_id].rename(columns={'read_count': sample_id})
                                 , how='inner'
                                 , on=['chr', 'gene'])

    logging.info('Successfully merged sample read counts -- shape: {0}'.format(read_count_df.shape))

    del read_count_dict

    # subset genes, exons to genes in intersection of experiments.
    genes_df = genes_df[genes_df.gene.isin(read_count_df.gene.unique())]
    exon_df = exon_df[exon_df.gene.isin(read_count_df.gene.unique())]

    # re-order genes, read counts so that they're parsimonious in chromosomes + gene ordering.
    genes_df.sort_values(['chr', 'gene'], inplace=True)
    genes_df.reset_index(inplace=True, drop=True)
    read_count_df.sort_values(['chr', 'gene'], inplace=True)
    read_count_df.reset_index(inplace=True, drop=True)

    # quality control.
    if genes_df.shape[0] != read_count_df.shape[0]:
        raise ValueError('Genes DataFrame and read counts DataFrame do not have same number of rows!')

    # ---------------------------------------------------------------------------- #
    # Slice up genome coverage matrix for each gene according to exon positioning.
    # Run in parallel over chromosomes.
    # ---------------------------------------------------------------------------- #
    # cov_output_dir = None if args.disregard_coverage else output_dir
    p = mp.Pool(processes=n_jobs)
    gene_cov_mats = [p.apply_async(gene_coverage
                                   , args=(exon_df, chrom, cov_files, output_dir, True)) for chrom in chroms]
    p.close()
    gene_cov_mats = [x.get() for x in gene_cov_mats]

    # ---------------------------------------------------------------------------- #
    # Process gene coverage matrix output prior to running NMF.
    # ---------------------------------------------------------------------------- #
    # convert list of tuples into 2-d dictionary: {chrom: {gene: coverage matrix}}, and
    # initialize an OrderedDict to store coverage matrices that pass DegNorm inclusion
    # criteria in order.
    chrom_gene_cov_dict = {gene_cov_mats[i][1]: gene_cov_mats[i][0] for i in range(len(gene_cov_mats))}
    gene_cov_dict = OrderedDict()

    # Determine for which genes to run DegNorm, and for genes where we will run DegNorm,
    # which transcript regions to filter out prior to running DegNorm.
    logging.info('Determining genes to include in DegNorm coverage curve approximation.')
    delete_idx = list()

    for i in range(genes_df.shape[0]):
        chrom = genes_df.chr.iloc[i]
        gene = genes_df.gene.iloc[i]

        # if user-defined gene subset is specified, only add gene if in subset.
        if args.genes:
            if gene not in args.genes:
                delete_idx.append(i)
                continue

        # extract gene's p x Li coverage matrix.
        cov_mat = chrom_gene_cov_dict[chrom][gene]

        # do not add gene if there are any 100%-zero coverage samples.
        if any(cov_mat.sum(axis=1) == 0):
            delete_idx.append(i)

        else:
            gene_cov_dict[gene] = cov_mat

    if delete_idx:
        read_count_df = read_count_df.drop(delete_idx
                                           , axis=0).reset_index(drop=True)
        genes_df = genes_df.drop(delete_idx
                                 , axis=0).reset_index(drop=True)

    # quality control.
    if (read_count_df.shape[0] == 0) or (genes_df.empty) or (len(gene_cov_dict) == 0):
        raise ValueError('No genes available to run through DegNorm!\n'
                         'Check that your requested genes are in genome annotation file.')

    # check that read counts and coverage matrices contain data for same number of genes.
    if len(gene_cov_dict.keys()) != read_count_df.shape[0]:
        raise ValueError('Number of coverage matrices not equal to number of genes in read count DataFrame!')

    logging.info('DegNorm will run on {0} genes.'.format(len(gene_cov_dict)))

    # save gene annotation metadata.
    exon_output_file = os.path.join(output_dir, 'gene_exon_metadata.csv')
    logging.info('Saving gene-exon metadata to {0}'.format(exon_output_file))
    exon_df.to_csv(exon_output_file
                   , index=False)

    # save read counts.
    read_count_file = os.path.join(output_dir, 'read_counts.csv')
    logging.info('Saving original read counts to {0}'.format(read_count_file))
    read_count_df.to_csv(read_count_file
                         , index=False)

    # free up more memory: delete exon / genome annotation data
    del chrom_gene_cov_dict, gene_cov_mats

    # ---------------------------------------------------------------------------- #
    # Run NMF.
    # ---------------------------------------------------------------------------- #
    logging.info('Executing NMF-OA over-approximation algorithm...')
    nmfoa = GeneNMFOA(iter=args.iter
                      , nmf_iter=args.nmf_iter
                      , grid_points=args.downsample_rate
                      , n_jobs=n_jobs)
    nmfoa.fit_transform(gene_cov_dict
                        , reads_dat=read_count_df[sample_ids].values)

    # ---------------------------------------------------------------------------- #
    # Save results.
    # ---------------------------------------------------------------------------- #
    logging.info('Saving NMF-OA output:\n'
                 '-- degradation index scores --\n'
                 '-- adjusted read counts --\n'
                 '-- coverage curve estimates --')
    nmfoa.save_results(genes_df
                       , output_dir=output_dir
                       , sample_ids=sample_ids)

    # ---------------------------------------------------------------------------- #
    # Generate coverage curve plots.
    # ---------------------------------------------------------------------------- #
    logging.info('Generating coverage curve plots.')

    # subset chromosomes to just those corresponding to genes run through degnorm
    chroms = genes_df[genes_df.gene.isin(nmfoa.genes)].chr.unique().tolist()

    p = mp.Pool(processes=n_jobs)
    out = [p.apply_async(save_chrom_coverage
                         , args=(os.path.join(output_dir, chrom, 'coverage_matrices_{0}.pkl'.format(chrom))
                                 , os.path.join(output_dir, chrom, 'estimated_coverage_matrices_{0}.pkl'.format(chrom))
                                 , exon_df[exon_df.chr == chrom]
                                 , sample_ids
                                 , [10, 6]
                                 , os.path.join(output_dir, chrom))) for chrom in chroms]
    p.close()

    # Execute parallel work.
    out = [x.get() for x in out]

    # ---------------------------------------------------------------------------- #
    # Run summary report and exit.
    # ---------------------------------------------------------------------------- #
    logging.info('Rendering summary report.')
    render_report(data_dir=output_dir
                  , genenmfoa=nmfoa
                  , gene_manifest_df=genes_df
                  , input_files=args.input_files
                  , sample_ids=sample_ids
                  , top_n_genes=5
                  , output_dir=output_dir)

    logging.info('DegNorm pipeline complete! Exiting...')
    sys.exit(0)


if __name__ == "__main__":
    main()