#!/usr/bin/env python3

# ---------------------------------------------------------------------------- #
# DegNorm CLI entrypoint for use on a single node with hyperthreading.
# ---------------------------------------------------------------------------- #

import sys
from degnorm.reads import *
from degnorm.coverage import *
from degnorm.gene_processing import *
from degnorm.data_access import *
from degnorm.nmf import *
from degnorm.warm_start import *
from degnorm.report import render_report
from collections import OrderedDict


def main():

    # ---------------------------------------------------------------------------- #
    # Load CLI arguments, display welcome message, create output directory
    # ---------------------------------------------------------------------------- #
    args = parse_args()
    n_jobs = args.proc_per_node
    output_dir = create_output_dir(args.output_dir)
    configure_logger(output_dir)
    welcome()
    logging.info('DegNorm output directory -- {0}'.format(output_dir))

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
        genes_df = load_dat['genes_df']
        sample_ids = load_dat['sample_ids']

    # ---------------------------------------------------------------------------- #
    # Path 2: .bam file preprocessing path.
    # Determine intersection of chromosomes across samples from .bam files
    # ---------------------------------------------------------------------------- #
    else:
        sample_ids = list()
        chroms = list()
        cov_files = OrderedDict()
        read_count_dict = dict()
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
        # Load .bam files while simultaneously parsing into coverage arrays. Store chromosome
        # coverage arrays and and compute gene read counts.
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

        # ensure chromosomes still reflective of genes to be processed.
        chroms = read_count_df.chr.unique().tolist()

        # ---------------------------------------------------------------------------- #
        # Slice up genome coverage matrix for each gene according to exon positioning.
        # Run in parallel over chromosomes.
        # ---------------------------------------------------------------------------- #
        gene_cov_dicts = Parallel(n_jobs=min(n_jobs, len(chroms))
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
        # do not add gene if maximum coverage is < minimum maximum coverage threshold.
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

    # briefly summarize DegNorm input and settings.
    logging.info('RNA-seq sample identifiers: \n\t' + ', '.join(sample_ids))
    logging.info('DegNorm will run on {0} genes, downsampling rate = 1 / {1}, {2} baseline selection.'
                 .format(len(gene_cov_dict), args.downsample_rate, 'without' if args.skip_baseline_selection else 'with'))

    # ---------------------------------------------------------------------------- #
    # Run NMF.
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


if __name__ == "__main__":
    main()