from degnorm.reads import *
from degnorm.coverage import *
from degnorm.gene_processing import *
from degnorm.visualizations import *
from degnorm.nmf import *
from degnorm.warm_start import *
from degnorm.report import render_report
from collections import OrderedDict
import sys


def main():

    # ---------------------------------------------------------------------------- #
    # Load CLI arguments and create output directory
    # ---------------------------------------------------------------------------- #
    args = parse_args()
    n_jobs = args.cpu
    output_dir = create_output_dir(args.output_dir)

    # ---------------------------------------------------------------------------- #
    # Path 1: warm-start path.
    # If supplied a warm-start directory, load and copy previously parsed coverage matrices,
    # genome annotation data, read counts into new output dir.
    # ---------------------------------------------------------------------------- #
    if args.warm_start_dir:
        logging.info('WARM-START: loading data from previous DegNorm run contained in {0}'
                     .format(args.warm_start_dir))

        load_dat = load_from_previous(args.warm_start_dir
                                      , new_dir=output_dir)
        chrom_gene_cov_dict = load_dat['chrom_gene_cov_dict']
        read_count_df = load_dat['read_count_df']
        exon_df = load_dat['exon_df']
        genes_df = load_dat['genes_df']
        sample_ids = load_dat['sample_ids']

    # ---------------------------------------------------------------------------- #
    # Path 2: preprocessing path.
    # If supplied with .bam files (NOT using a warm-start),
    # convert them to .sam; further, determine intersection of chromosomes across samples.
    # ---------------------------------------------------------------------------- #
    else:

        sample_ids = list()
        chroms = list()
        cov_files = dict()
        read_count_dict = dict()
        n_samples = len(args.input_files)

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
                                      , chroms=chroms)
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
        gc.collect()

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
        gc.collect()

        # ---------------------------------------------------------------------------- #
        # Save gene annotation metadata and original read counts.
        # ---------------------------------------------------------------------------- #

        # subset genes, exons to genes in intersection of experiments.
        genes_df = genes_df[genes_df.gene.isin(read_count_df.gene.unique())]
        exon_df = exon_df[exon_df.gene.isin(read_count_df.gene.unique())]

        # re-order genes, read counts so that they're parsimonious in chromosome + gene ordering.
        genes_df.sort_values(['chr', 'gene'], inplace=True)
        genes_df.reset_index(inplace=True, drop=True)
        read_count_df.sort_values(['chr', 'gene'], inplace=True)
        read_count_df.reset_index(inplace=True, drop=True)

        # quality control.
        if genes_df.shape[0] != read_count_df.shape[0]:
            raise ValueError('Genes DataFrame and read counts DataFrame do not have same number of rows!')

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

        # ---------------------------------------------------------------------------- #
        # Slice up genome coverage matrix for each gene according to exon positioning.
        # Run in parallel over chromosomes.
        # ---------------------------------------------------------------------------- #
        gene_cov_dicts = Parallel(n_jobs=n_jobs
                       , verbose=0
                       , backend='threading')(delayed(gene_coverage)(
            exon_df=exon_df,
            chrom=chrom,
            coverage_files=cov_files,
            output_dir=output_dir,
            verbose=True) for chrom in chroms)

        # convert list of tuples into 2-d dictionary: {chrom: {gene: coverage matrix}}
        chrom_gene_cov_dict = {chroms[idx]: gene_cov_dicts[idx] for idx in range(len(chroms))}

        del gene_cov_dicts

    # ---------------------------------------------------------------------------- #
    # PATHS MERGE: warm-start and input-file paths meet.
    # Process gene coverage matrices prior to running NMF.
    # ---------------------------------------------------------------------------- #

    # initialize an OrderedDict to store coverage matrices that pass DegNorm inclusion
    # criteria in order.
    gene_cov_dict = OrderedDict()

    # Determine for which genes to run DegNorm, and for genes where we will run DegNorm,
    # which transcript regions to filter out prior to running DegNorm.
    logging.info('Determining genes to include in DegNorm coverage curve approximation.')
    delete_idx = list()

    for i in range(genes_df.shape[0]):
        chrom = genes_df.chr.iloc[i]
        gene = genes_df.gene.iloc[i]

        # extract gene's p x Li coverage matrix.
        cov_mat = chrom_gene_cov_dict[chrom][gene]

        # do not add gene if there are any 100%-zero coverage samples.
        # do not add gene if maximum coverage is < minimum coverage threshold.
        # do not add gene if downsample rate low enough s.t. take-every > length of gene.
        if any(cov_mat.sum(axis=1) == 0) or (cov_mat.max() < args.minimax_coverage) \
                or (cov_mat.shape[1] <= args.downsample_rate):
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

    # free up more memory: delete exon / genome annotation data
    del chrom_gene_cov_dict
    gc.collect()

    logging.info('DegNorm will run on {0} genes with downsampling rate = 1 / {1}'
                 .format(len(gene_cov_dict), args.downsample_rate))

    # ---------------------------------------------------------------------------- #
    # Run NMF.
    # ---------------------------------------------------------------------------- #

    # joblib overhead: specify temp folder if not in environment.
    joblib_folder = os.environ.get('JOBLIB_TEMP_FOLDER')
    if not joblib_folder:
        os.environ['JOBLIB_TEMP_FOLDER'] = output_dir

    logging.info('Executing NMF-OA over-approximation algorithm...')
    nmfoa = GeneNMFOA(iter=args.iter
                      , nmf_iter=args.nmf_iter
                      , downsample_rate=args.downsample_rate
                      , n_jobs=n_jobs)
    nmfoa.fit_transform(gene_cov_dict
                        , reads_dat=read_count_df[sample_ids].values)

    # restore original environment.
    if not joblib_folder:
        del os.environ['JOBLIB_TEMP_FOLDER']

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
    # Generate coverage curve plots if -g/--genes were specified.
    # ---------------------------------------------------------------------------- #
    if args.plot_genes:
        plot_genes = np.intersect1d(args.plot_genes, nmfoa.genes)

        if len(plot_genes) > 0:

            plot_exon_df = exon_df[exon_df.gene.isin(plot_genes)]
            logging.info('Generating coverage curve plots for specified genes.')

            out = Parallel(n_jobs=n_jobs
                           , verbose=0
                           , backend='threading')(delayed(save_chrom_coverage)(
                coverage_file=os.path.join(output_dir, chrom, 'coverage_matrices_{0}.pkl'.format(chrom)),
                estimates_file=os.path.join(output_dir, chrom, 'estimated_coverage_matrices_{0}.pkl'.format(chrom)),
                exon_df=plot_exon_df[plot_exon_df.chr == chrom],
                sample_ids=sample_ids,
                figsize=[10, 6],
                output_dir=os.path.join(output_dir, chrom))
                for chrom in plot_exon_df.chr.unique())

    # ---------------------------------------------------------------------------- #
    # Run summary report and exit.
    # ---------------------------------------------------------------------------- #
    logging.info('Rendering DegNorm summary report.')
    render_report(data_dir=output_dir
                  , genenmfoa=nmfoa
                  , input_files=args.input_files
                  , sample_ids=sample_ids
                  , top_n_genes=5
                  , output_dir=output_dir)

    logging.info('DegNorm pipeline complete! Exiting...')
    sys.exit(0)


if __name__ == "__main__":
    main()