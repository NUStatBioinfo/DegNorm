#!/usr/bin/env python3

# ---------------------------------------------------------------------------- #
# DegNorm CLI entrypoint for use on a single node with hyperthreading.
# ---------------------------------------------------------------------------- #

from degnorm.reads import *
from degnorm.reads_coverage_merge import *
from degnorm.gene_processing import *
from degnorm.data_access import *
from degnorm.nmf import *
from degnorm.warm_start import *
from degnorm.report import render_report


def main():

    # ---------------------------------------------------------------------------- #
    # Load CLI arguments, display welcome message, create output directory
    # ---------------------------------------------------------------------------- #
    args = parse_args()
    n_jobs = args.proc_per_node
    unique_alignments = not args.non_unique_alignments
    output_dir = create_output_dir(args.output_dir)
    configure_logger(output_dir)
    welcome()
    logging.info('DegNorm output directory -- {0}'.format(output_dir))

    # ---------------------------------------------------------------------------- #
    # If any Bam index (.bai) files need to be created, do that now.
    # ---------------------------------------------------------------------------- #
    if args.create_bai_files:
        for file_idx in range(len(args.create_bai_files)):
            bam_file = args.create_bai_files[file_idx]
            logging.info('creating index file for {0} -- {1} / {2}'
                         .format(bam_file, file_idx + 1, len(args.create_bai_files)))
            out = create_index_file(bam_file)

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
        gene_cov_dict = load_dat['gene_cov_dict']
        read_count_df = load_dat['read_count_df']
        genes_df = load_dat['genes_df']
        sample_ids = load_dat['sample_ids']

    # ---------------------------------------------------------------------------- #
    # Path 2: .bam file preprocessing path.
    # Determine intersection of chromosomes across samples from .bam files
    # ---------------------------------------------------------------------------- #
    else:
        sample_ids = list()
        chroms = list()
        n_samples = len(args.bam_files)

        # load each .bam file's header, find joint intersection of read chromosomes.
        for idx in range(n_samples):
            header_dat = BamReadsProcessor(args.bam_files[idx]
                                            , index_file=args.bai_files[idx]).header
            new_chroms = header_dat.chr.values.tolist()

            if not chroms:
                chroms = new_chroms

            else:
                chroms = np.intersect1d(chroms, new_chroms).tolist()

        # ---------------------------------------------------------------------------- #
        # Load .gtf or .gff files and run processing pipeline.
        # Run in parallel over chromosomes.
        # ---------------------------------------------------------------------------- #
        logging.info('Begin genome annotation file processing...')
        gap = GeneAnnotationProcessor(args.genome_annotation
                                      , verbose=True
                                      , chroms=chroms)
        exon_df = gap.run()

        # take intersection of chromosomes available in genome annotation file and those in the reads data,
        # if for some reason annotation file only contains subset.
        chroms = np.intersect1d(chroms, exon_df.chr.unique()).tolist()

        # subset exon (and therefore gene) data based on chromosome set.
        exon_df = exon_df[exon_df.chr.isin(chroms)]
        genes_df = exon_df[['chr', 'gene', 'gene_start', 'gene_end']].drop_duplicates().reset_index(drop=True)

        logging.info('Found {0} chromosomes in intersection of all experiments and gene annotation data:\n'
                     '\t{1}'.format(len(chroms), ', '.join(chroms)))

        # break down gene overlap structures by chromosome, will need to feed it to coverage_read_counts method.
        logging.info('Determining gene overlap structure across chromosomes.')
        gene_overlap_dict = dict()
        n_overlap = 0
        n_isolated = 0
        for chrom in chroms:
            gene_overlap_dict[chrom] = get_gene_overlap_structure(subset_to_chrom(genes_df
                                                                                  , chrom=chrom))
            # compute fraction of overlap genes to total genes.
            if gene_overlap_dict[chrom].get('overlap_genes') is not None:
                n_overlap += np.sum([len(x) for x in gene_overlap_dict[chrom]['overlap_genes']])

            if gene_overlap_dict[chrom].get('isolated_genes') is not None:
                n_isolated += len(gene_overlap_dict[chrom]['isolated_genes'])

        logging.info('Rate of gene overlap: {0} / {1}'.format(n_overlap, n_isolated + n_overlap))

        # ---------------------------------------------------------------------------- #
        # Load .bam files and parse them into coverage arrays, read counts.
        # ---------------------------------------------------------------------------- #

        # iterate over .bam files; compute each sample's chromosomes' coverage arrays
        # and save them to .npz files.
        for idx in range(n_samples):
            logging.info('Loading RNA-seq data file {0} / {1}'.format(idx + 1, n_samples))

            reader = BamReadsProcessor(bam_file=args.bam_files[idx]
                                       , index_file=args.bai_files[idx]
                                       , chroms=chroms
                                       , n_jobs=n_jobs
                                       , output_dir=output_dir
                                       , unique_alignment=unique_alignments
                                       , verbose=True)

            sample_id = reader.sample_id
            sample_ids.append(sample_id)

            # run simultaneous coverage, read counting procedure on alignment file.
            reader.coverage_read_counts(gene_overlap_dict
                                        , gene_df=genes_df
                                        , exon_df=exon_df)

        logging.info('Successfully processed chromosome read coverage and gene read counts for all {0} experiments'
                     .format(len(sample_ids)))

        del reader
        gc.collect()

        # ---------------------------------------------------------------------------- #
        # Merge, load per-sample files:
        # 1. obtain read count DataFrame containing X, an n (genes) x p (samples) matrix.
        # 2. per-sample gene coverage matrices,
        #    and save them back to disk in .pkl files on per-chromosome basis.
        # ---------------------------------------------------------------------------- #
        logging.info('Merging read counts across samples.')
        read_count_df = merge_read_counts(output_dir
                                          , sample_ids=sample_ids
                                          , chroms=chroms)
        logging.info('Read counts merge successful. Read count data shape: {0}'.format(read_count_df.shape))

        logging.info('Merging gene coverage arrays across samples and saving results to chromosome directories.')
        gene_cov_dict = merge_coverage(output_dir
                                       , sample_ids=sample_ids
                                       , exon_df=exon_df
                                       , n_jobs=n_jobs
                                       , output_dir=output_dir
                                       , verbose=True)

        logging.info('Complete coverage merge successful. Number of loaded coverage arrays: {0}'
                     .format(len(gene_cov_dict)))

        # remove per-sample raw sample coverage, read count files.
        for s_id in sample_ids:
            shutil.rmtree(os.path.join(output_dir, s_id))

        # ---------------------------------------------------------------------------- #
        # Save gene annotation metadata and original read counts.
        # ---------------------------------------------------------------------------- #
        genes = list(gene_cov_dict.keys())

        # order genes and read counts by order in coverage set by
        # setting index to gene names, loc[genes], making genes a column again by dropping index.
        genes_df.set_index('gene'
                           , inplace=True)
        read_count_df.set_index('gene'
                                , inplace=True)

        genes_df = genes_df.loc[genes]
        read_count_df = read_count_df.loc[genes]

        genes_df.reset_index(drop=False
                             , inplace=True)
        read_count_df.reset_index(drop=False
                                  , inplace=True)

        # only keep exons for genes in coverage set.
        exon_df = exon_df[exon_df.gene.isin(genes)]

        # quality control.
        if genes_df.shape[0] != read_count_df.shape[0]:
            raise ValueError('Genes DataFrame and read counts DataFrame do not have same number of rows!')

        # save gene annotation metadata.
        exon_output_file = os.path.join(output_dir, 'gene_exon_metadata.csv')
        logging.info('Saving gene-exon metadata to {0}'.format(exon_output_file))
        exon_df.to_csv(exon_output_file
                       , index=False)

        # save merged read counts.
        read_count_file = os.path.join(output_dir, 'read_counts.csv')
        logging.info('Saving original read counts to {0}'.format(read_count_file))
        read_count_df.to_csv(read_count_file
                             , index=False)

    # ---------------------------------------------------------------------------- #
    # PATHS MERGE: warm-start and input-file paths meet.
    # Process gene coverage matrices prior to running NMF.
    # ---------------------------------------------------------------------------- #

    # Determine for which genes to run DegNorm, and for genes where we will run DegNorm,
    # which transcript regions to filter out prior to running DegNorm.
    logging.info('Determining genes to include in DegNorm coverage curve approximation.')
    delete_idx = list()

    for i in range(genes_df.shape[0]):
        gene = genes_df.gene.iloc[i]

        # extract gene's p x Li coverage matrix.
        cov_mat = gene_cov_dict[gene]

        # do not run gene if maximum coverage is < minimum maximum coverage threshold.
        # do not run gene if downsample rate low enough s.t. take-every > length of gene.
        if (cov_mat.max() < args.minimax_coverage) or (cov_mat.shape[1] <= args.downsample_rate):
            delete_idx.append(i)
            del gene_cov_dict[gene]

    # drop genes from read counts, gene set if coverage was non conformant.
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

    # briefly summarize DegNorm input and settings.
    logging.info('RNA-seq sample identifiers: \n\t' + ', '.join(sample_ids))
    logging.info('DegNorm will run on {0} genes, downsampling rate = 1 / {1}, {2} baseline selection.'
                 .format(len(gene_cov_dict), args.downsample_rate, 'without' if args.skip_baseline_selection else 'with'))

    # ---------------------------------------------------------------------------- #
    # Run NMF-OA.
    # ---------------------------------------------------------------------------- #

    # joblib overhead: specify temp folder if not in environment.
    joblib_folder = os.environ.get('JOBLIB_TEMP_FOLDER')
    if not joblib_folder:
        os.environ['JOBLIB_TEMP_FOLDER'] = output_dir

    logging.info('Executing NMF-OA over-approximation algorithm...')
    nmfoa = GeneNMFOA(degnorm_iter=args.iter
                      , nmf_iter=args.nmf_iter
                      , downsample_rate=args.downsample_rate
                      , n_jobs=n_jobs
                      , skip_baseline_selection=args.skip_baseline_selection)
    estimates = nmfoa.run(gene_cov_dict
                          , reads_dat=read_count_df[sample_ids].values.astype(np.float_))

    # restore original environment.
    if not joblib_folder:
        del os.environ['JOBLIB_TEMP_FOLDER']

    # ---------------------------------------------------------------------------- #
    # Save results.
    # ---------------------------------------------------------------------------- #
    logging.info('Saving NMF-OA output:\n'
                 '\t-- degradation index scores --\n'
                 '\t-- adjusted read counts --\n'
                 '\t-- coverage curve estimates --')
    nmfoa.save_results(estimates
                       , gene_manifest_df=genes_df
                       , output_dir=output_dir
                       , sample_ids=sample_ids)

    # ---------------------------------------------------------------------------- #
    # Generate coverage curve plots if -g/--genes were specified.
    # ---------------------------------------------------------------------------- #
    if args.plot_genes:
        plot_genes = np.intersect1d(args.plot_genes, nmfoa.genes)

        if len(plot_genes) > 0:
            logging.info('Generating coverage curve plots for specified genes.')
            out = get_coverage_plots(plot_genes
                                     , degnorm_dir=output_dir
                                     , figsize=[10, 6]
                                     , save_dir=output_dir)

    # ---------------------------------------------------------------------------- #
    # Run summary report and exit.
    # ---------------------------------------------------------------------------- #
    logging.info('Rendering DegNorm summary report.')
    degnorm_dat = {'degnorm_iter': args.iter
                   , 'nmf_iter': args.nmf_iter
                   , 'downsample_rate': args.downsample_rate
                   , 'rho': nmfoa.rho
                   , 'genes': nmfoa.genes}

    render_report(data_dir=output_dir
                  , degnorm_data=degnorm_dat
                  , bam_files=args.bam_files if not args.warm_start_dir else [args.warm_start_dir]
                  , sample_ids=sample_ids
                  , top_n_genes=5
                  , output_dir=output_dir)

    logging.info('DegNorm pipeline complete! Exiting...')
    sys.exit(0)


if __name__ == '__main__':
    main()
