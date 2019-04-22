#!/usr/bin/env python3

# ---------------------------------------------------------------------------- #
# DegNorm CLI entrypoint for distributed pipeline run.
# ---------------------------------------------------------------------------- #

try:
    from mpi4py import MPI
except ImportError as e:
    raise e

from degnorm.reads import *
from degnorm.reads_coverage_merge import *
from degnorm.gene_processing import *
from degnorm.data_access import *
from degnorm.nmf_mpi import *
from degnorm.warm_start import *
from degnorm.report import render_report
import sys


COMM = MPI.COMM_WORLD
SIZE = COMM.size
RANK = COMM.rank
HOST = MPI.Get_processor_name()
WORKER_NODES = list(range(SIZE))

# ensure we have at least 2 nodes.
if SIZE <= 1:
    raise RuntimeError('Must run degnorm_mpi with at least 2 workers.')


def mpi_logging_info(msg):
    """
    Format a message from an MPI worker to a pre-configured stream logger, so that worker rank
    and comm size are logged as well.

    :param msg: message to send to logger
    """
    logging.info('({rank}/{size}) -- {msg}'.format(rank=RANK + 1, size=SIZE, msg=msg))


def main():

    # ---------------------------------------------------------------------------- #
    # Load CLI arguments (everyone)
    # display welcome message, create output directory (just master node)
    # ---------------------------------------------------------------------------- #
    args = parse_args()

    # number of processes to spawn within a compute node.
    n_jobs = args.proc_per_node

    # determine if we're only keeping uniquely-mapped reads
    unique_alignments = not args.non_unique_alignments

    # everyone sets up stdout logger.
    configure_logger(output_dir=None
                     , mpi=True)

    # master welcomes user, sets up output directory.
    if RANK == 0:
        welcome()
        output_dir = create_output_dir(args.output_dir)
        mpi_logging_info('DegNorm save directory: {0}'.format(output_dir))

    else:
        output_dir = None

    # Give everyone the output directory.
    output_dir = COMM.bcast(output_dir, root=0)

    # ---------------------------------------------------------------------------- #
    # If Bam index (.bai) files need to be created, do that now.
    # ---------------------------------------------------------------------------- #
    if args.create_bai_files:

        # distribute .bai files to workers.
        bai_file_chunks = split_into_chunks(args.create_bai_files
                                            , n=SIZE)

        # worker creates its .bai files.
        if RANK < len(bai_file_chunks):
            my_bai_files = bai_file_chunks[RANK]

            for file_idx in range(len(my_bai_files)):
                bam_file = my_bai_files[file_idx]
                mpi_logging_info('creating index file for {0}'.format(bam_file))
                out = create_index_file(bam_file)

        # have everyone wait up until .bai files are created.
        COMM.Barrier()

    # ---------------------------------------------------------------------------- #
    # Path 1: warm-start path.
    # If supplied a warm-start directory, load and copy previously parsed coverage matrices,
    # genome annotation data, read counts into new output dir.
    # ---------------------------------------------------------------------------- #
    if args.warm_start_dir:

        # only master needs genes, reads, and coverage matrix dictionary.
        if RANK == 0:
            mpi_logging_info('WARM-START: loading data from previous DegNorm run contained in {0}'
                             .format(args.warm_start_dir))
            load_dat = load_from_previous(args.warm_start_dir
                                          , new_dir=output_dir)
            gene_cov_dict = load_dat['gene_cov_dict']
            read_count_df = load_dat['read_count_df']
            genes_df = load_dat['genes_df']
            sample_ids = load_dat['sample_ids']

        else:
            sample_ids = None

        # Everyone needs sample_ids.
        sample_ids = COMM.bcast(sample_ids, root=0)

    # ---------------------------------------------------------------------------- #
    # Path 2: .bam file preprocessing path.
    # Determine intersection of chromosomes across samples from .bam files
    # ---------------------------------------------------------------------------- #
    else:
        n_samples = len(args.bam_files)

        # Have master node assess chromosomes to be included in DegNorm pipeline run.
        if RANK == 0:
            chroms = list()

            # load each .bam file's header, find joint intersection of read chromosomes.
            for idx in range(n_samples):
                header_dat = BamReadsProcessor(args.bam_files[idx]
                                                , index_file=args.bai_files[idx]).header
                new_chroms = header_dat.chr.values.tolist()

                if not chroms:
                    chroms = new_chroms

                else:
                    chroms = np.intersect1d(chroms, new_chroms).tolist()

        else:
            chroms = None

        # broadcast the chromosomes of interest.
        chroms = COMM.bcast(chroms, root=0)

        # ---------------------------------------------------------------------------- #
        # Load .gtf or .gff files and run processing pipeline.
        # Run in parallel on the master node over chromosomes.
        # (master does this)
        # ---------------------------------------------------------------------------- #
        if RANK == 0:
            mpi_logging_info('Begin genome annotation file processing...')
            gap = GeneAnnotationProcessor(args.genome_annotation
                                          , verbose=True
                                          , chroms=chroms)
            exon_df = gap.run()

            # take intersection of chromosomes available in genome annotation file and those in the reads data,
            # if for some reason annotation file only contains subset.
            chroms = np.intersect1d(chroms, exon_df.chr.unique()).tolist()
            exon_df = exon_df[exon_df.chr.isin(chroms)]
            genes_df = exon_df[['chr', 'gene', 'gene_start', 'gene_end']].drop_duplicates().reset_index(drop=True)

            mpi_logging_info('Found {0} chromosomes in intersection of all experiments and gene annotation data:\n'
                             '\t{1}'.format(len(chroms), ', '.join(chroms)))

        else:
            chroms = None
            exon_df = None
            genes_df = None

        # Broadcast the potentially updated set of chromosomes of interest, exons, and genes.
        chroms = COMM.bcast(chroms, root=0)
        exon_df = COMM.bcast(exon_df, root=0)
        genes_df = COMM.bcast(genes_df, root=0)

        # ---------------------------------------------------------------------------- #
        # For each chromosome, find groups of mutually overlapping genes and
        # groups of non-overlapping genes (isolates). Every worker will need
        # every chromosome's overlap structure to run BamReadsProcessor.coverage_read_counts method.
        # ---------------------------------------------------------------------------- #
        if RANK == 0:
            mpi_logging_info('Determining gene overlap structure for all {0} chromosomes'.format(len(chroms)))

        # divvy up chromosomes across workers for gene overlap structure work.
        my_chroms = split_into_chunks(chroms
                                      , n=SIZE)

        # storage for each worker's chromosome/gene overlap structure data
        gene_overlap_dict = dict()

        # distributed processing of alignment files
        if RANK < len(my_chroms):
            for chrom in my_chroms[RANK]:
                gene_overlap_dict[chrom] = get_gene_overlap_structure(subset_to_chrom(genes_df
                                                                                      , chrom=chrom))

        # everyone sends their gene overlap structure data to master.
        gene_overlap_dict = COMM.gather(gene_overlap_dict, root=0)

        # master concatenates chromosome gene overlap structures, sends back out to workers.
        if RANK == 0:

            # collapse list of gene overlap structure dictionaries into one dictionary.
            gene_overlap_dict = {k: v for d in gene_overlap_dict for k, v in d.items()}

            # display rate of gene overlap.
            n_overlap = 0
            n_isolated = 0
            for chrom in chroms:
                if gene_overlap_dict[chrom].get('overlap_genes') is not None:
                    n_overlap += np.sum([len(x) for x in gene_overlap_dict[chrom]['overlap_genes']])

                if gene_overlap_dict[chrom].get('isolated_genes') is not None:
                    n_isolated += len(gene_overlap_dict[chrom]['isolated_genes'])

            mpi_logging_info('Rate of gene overlap: {0} / {1}'.format(n_overlap, n_isolated + n_overlap))

        # master broadcasts out gene overlap structure dictionary to all workers.
        gene_overlap_dict = COMM.bcast(gene_overlap_dict, root=0)

        # make sure everyone is caught up and has gene_overlap_dict.
        COMM.Barrier()

        # ---------------------------------------------------------------------------- #
        # Distribute .bam files and parse them into coverage arrays, read counts.
        # (Each worker gets ~ same number of .bam files to process.)
        # ---------------------------------------------------------------------------- #

        # sample IDs (to be gathered up by master later)
        sample_ids = list()

        # iterate over node's .bam/.bai files; compute each sample's chromosomes' coverage arrays
        # and save them to .npz files.
        my_files_idx = split_into_chunks(range(n_samples)
                                         , n=SIZE)

        # iterate over .bam files; compute each sample's chromosomes' coverage arrays
        # and save them to .npz files.
        if RANK < len(my_files_idx):
            for idx in my_files_idx[RANK]:
                mpi_logging_info('Loading RNA-seq data file {0} -- {1}/{2}'
                                 .format(args.bam_files[idx], idx + 1, n_samples))

                reader = BamReadsProcessor(bam_file=args.bam_files[idx]
                                           , index_file=args.bai_files[idx]
                                           , chroms=chroms
                                           , n_jobs=n_jobs
                                           , output_dir=output_dir
                                           , unique_alignment=unique_alignments
                                           , verbose=True)

                sample_ids.append(reader.sample_id)

                # run simultaneous coverage, read counting procedure on alignment file.
                reader.coverage_read_counts(gene_overlap_dict
                                            , gene_df=genes_df
                                            , exon_df=exon_df)

                del reader
                gc.collect()

        # everyone gives their sample_ids to master.
        sample_ids = COMM.gather(sample_ids, root=0)

        # ---------------------------------------------------------------------------- #
        # Master to merge per-sample coverage, read counts across samples.
        # Obtain read count DataFrame, an n (genes) x p (samples) matrix,
        # gene_cov_dict, an OrderedDict containing coverage matrices for all n genes.
        # ---------------------------------------------------------------------------- #
        if RANK == 0:
            # collapse sample_ids from list of worker sample ID lists into a single 1-d list.
            # Important: this sets the order of sample IDs for the remainder of the DegNorm pipeline.
            sample_ids = flatten_2d(sample_ids).tolist()
            mpi_logging_info('Successfully computed gene coverage, read counts for all {0} experiments.'
                             .format(len(sample_ids)))
            mpi_logging_info('RNA-seq sample identifiers: \n\t' + ', '.join(sample_ids))

            mpi_logging_info('Merging read counts across samples.')
            read_count_df = merge_read_counts(output_dir
                                              , sample_ids=sample_ids
                                              , chroms=chroms)
            mpi_logging_info('Read counts merge successful. Read count data shape: {0}'.format(read_count_df.shape))

            mpi_logging_info('Merging gene coverage arrays across samples and saving results to chromosome directories.')
            gene_cov_dict = merge_coverage(output_dir
                                           , sample_ids=sample_ids
                                           , exon_df=exon_df
                                           , n_jobs=n_jobs
                                           , output_dir=output_dir
                                           , verbose=True)

            mpi_logging_info('Coverage merge successful. Number of loaded coverage matrices: {0}'
                             .format(len(gene_cov_dict)))

            # # remove raw sample coverage, read count files.
            for s_id in sample_ids:
                shutil.rmtree(os.path.join(output_dir, s_id))

            # ---------------------------------------------------------------------------- #
            # Reorder read counts, genes, exons according to the ordering of the genes
            # within gene_cov_dict.
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

            # ---------------------------------------------------------------------------- #
            # Save gene annotation metadata and original read counts.
            # ---------------------------------------------------------------------------- #

            # save gene annotation metadata.
            exon_output_file = os.path.join(output_dir, 'gene_exon_metadata.csv')
            mpi_logging_info('Saving gene-exon metadata to {0}'.format(exon_output_file))
            exon_df.to_csv(exon_output_file
                           , index=False)

            # save merged read counts.
            read_count_file = os.path.join(output_dir, 'read_counts.csv')
            mpi_logging_info('Saving original read counts to {0}'.format(read_count_file))
            read_count_df.to_csv(read_count_file
                                 , index=False)

        else:
            # only MASTER has correct genes_df, exon_df, read_count_df, gene_cov_dict.
            # Everyone else will get them later on.
            sample_ids = None

        # give everyone those sample_ids.
        sample_ids = COMM.bcast(sample_ids, root=0)

    # ---------------------------------------------------------------------------- #
    # PATHS MERGE: warm-start and input-file paths meet.
    # Process gene coverage matrices prior to running NMF.
    # ---------------------------------------------------------------------------- #
    if RANK == 0:

        # Determine for which genes to run DegNorm, and for genes where we will run DegNorm,
        # which transcript regions to filter out prior to running DegNorm.
        mpi_logging_info('Determining genes to include in DegNorm coverage curve approximation.')
        delete_idx = list()

        for i in range(genes_df.shape[0]):
            gene = genes_df.gene.iloc[i]

            # extract gene's p x Li coverage matrix. We got chrom_gene_cov_dict from prior gather.
            cov_mat = gene_cov_dict[gene]

            # do not add gene if maximum coverage is < minimum maximum coverage threshold.
            # do not add gene if it's unreasonably long, i.e. > 9 megabases.
            # do not add gene if downsample rate low enough s.t. take-every > length of gene.
            # do not add gene if max coverage is unreasonable, i.e. > 2^31.
            if (cov_mat.max() < args.minimax_coverage) or (cov_mat.shape[1] > 9e6) \
                    or (cov_mat.shape[1] <= args.downsample_rate) \
                    or (cov_mat.max() > 2147483647):

                delete_idx.append(i)
                del gene_cov_dict[gene]

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
        mpi_logging_info('DegNorm will run on {0} genes, downsampling rate = 1 / {1}, {2} baseline selection.'
                         .format(len(gene_cov_dict), args.downsample_rate, 'without' if args.skip_baseline_selection else 'with'))

        # write gene_cov_dict to disk, avoid MPI errors in trying to broadcast this large dictionary.
        with open(os.path.join(output_dir, 'TMP_gene_cov_dict.pkl'), 'wb') as f:
            pkl.dump(gene_cov_dict, f)

    else:
        read_count_df = None

    # master broadcasts read counts out to all workers.
    read_count_df = COMM.bcast(read_count_df, root=0)

    # everyone waits for master to finish dumping gene_cov_dict to disk.
    COMM.Barrier()

    # workers load up the pickled up gene coverage matrices saved by master.
    if RANK > 0:
        with open(os.path.join(output_dir, 'TMP_gene_cov_dict.pkl'), 'rb') as f:
            gene_cov_dict = pkl.load(f)

    # let's all wait up before going into DegNorm iterations.
    COMM.Barrier()

    # ---------------------------------------------------------------------------- #
    # Run NMF-OA.
    # ---------------------------------------------------------------------------- #
    if RANK == 0:
        # master deletes the saved coverage data that everyone has loaded by now.
        os.remove(os.path.join(output_dir, 'TMP_gene_cov_dict.pkl'))
        mpi_logging_info('Begin executing NMFOA algorithm...')

    nmfoa_output = run_gene_nmfoa_mpi(COMM
                                      , cov_dat=gene_cov_dict
                                      , reads_dat=read_count_df[sample_ids].values.astype(np.float_)
                                      , degnorm_iter=args.iter
                                      , nmf_iter=args.nmf_iter
                                      , downsample_rate=args.downsample_rate
                                      , n_jobs=n_jobs
                                      , skip_baseline_selection=args.skip_baseline_selection)

    # drop large data objects we don't need anymore.
    del gene_cov_dict, read_count_df

    # ---------------------------------------------------------------------------- #
    # Save results.
    # ---------------------------------------------------------------------------- #
    if RANK == 0:
        mpi_logging_info('Saving NMF-OA output:\n'
                         '\t-- degradation index scores --\n'
                         '\t-- adjusted read counts --\n'
                         '\t-- coverage curve estimates --\n')

        save_results(genes_df
                     , estimates=nmfoa_output['estimates']
                     , rho=nmfoa_output['rho']
                     , x_adj=nmfoa_output['x_adj']
                     , ran_baseline_selection=nmfoa_output['ran_baseline_selection']
                     , sample_ids=sample_ids
                     , output_dir=output_dir)

        genes = list(nmfoa_output['estimates'].keys())

        # split up any specified --plot-genes for parallel processing.
        if args.plot_genes:
            plot_genes = np.intersect1d(args.plot_genes, genes)
            plot_genes = split_into_chunks(plot_genes
                                           , n=SIZE)

            # if --plot-genes < number of workers...
            if len(plot_genes) < SIZE:
                plot_genes += [[None]] * (SIZE - len(plot_genes))

        else:
            plot_genes = [[None]] * SIZE

    else:
        plot_genes = None

    # master disseminates genes for plotting.
    plot_genes = COMM.scatter(plot_genes
                              , root=0)

    # ---------------------------------------------------------------------------- #
    # Generate coverage curve plots in parallel (only if --plot-genes specified)
    # ---------------------------------------------------------------------------- #
    if plot_genes[0] is not None:
        mpi_logging_info('Generating visualizations for --plot-genes')
        out = get_coverage_plots(plot_genes
                                 , degnorm_dir=output_dir
                                 , figsize=[10, 6]
                                 , save_dir=output_dir)

    # ---------------------------------------------------------------------------- #
    # Run summary report and exit.
    # ---------------------------------------------------------------------------- #
    if RANK == 0:
        mpi_logging_info('Rendering DegNorm summary report.')
        degnorm_dat = {'degnorm_iter': args.iter
                       , 'nmf_iter': args.nmf_iter
                       , 'downsample_rate': args.downsample_rate
                       , 'rho': nmfoa_output['rho']
                       , 'genes': list(nmfoa_output['estimates'].keys())}

        render_report(data_dir=output_dir
                      , degnorm_data=degnorm_dat
                      , bam_files=args.bam_files if not args.warm_start_dir else [args.warm_start_dir]
                      , sample_ids=sample_ids
                      , top_n_genes=5
                      , output_dir=output_dir)

        mpi_logging_info('DegNorm pipeline complete! Exiting...')

    COMM.Barrier()
    sys.exit(0)


if __name__ == '__main__':
    main()