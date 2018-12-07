#!/usr/bin/env python3

# ---------------------------------------------------------------------------- #
# DegNorm CLI entrypoint for use in a multi-node MPI-enabled setting.
# ---------------------------------------------------------------------------- #

try:
    from mpi4py import MPI
except ImportError as e:
    raise e

from degnorm.reads import *
from degnorm.coverage import *
from degnorm.gene_processing import *
from degnorm.data_access import *
from degnorm.nmf_mpi import *
from degnorm.warm_start import *
from degnorm.report import render_report
from collections import OrderedDict
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
    logging.info('({rank}/{size}) -- {msg}'.format(rank=RANK, size=SIZE, msg=msg))


def main():

    # ---------------------------------------------------------------------------- #
    # Load CLI arguments (everyone)
    # display welcome message, create output directory (just master node)
    # ---------------------------------------------------------------------------- #
    args = parse_args()
    n_samples = len(args.bam_files)

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
                logging.info('creating index file for {0}'
                             .format(bam_file, file_idx + 1, len(my_bai_files)))
                out = create_index_file(bam_file)

        # have everyone wait up until .bai files are created.
        COMM.Barrier()

    # assess number of processes to spawn within a compute node.
    n_jobs = min(args.proc_per_node, max_cpu() + 1) if args.proc_per_node else max_cpu()

    if RANK == 0:
        output_dir = create_output_dir(args.output_dir)
        configure_logger(output_dir=None
                         , mpi=True)
        welcome()
        mpi_logging_info('DegNorm save directory: {0}'.format(output_dir))

    else:
        output_dir = None

    # Broadcast the output directory.
    output_dir = COMM.bcast(output_dir, root=0)

    # ---------------------------------------------------------------------------- #
    # Path 1: warm-start path.
    # If supplied a warm-start directory, load and copy previously parsed coverage matrices,
    # genome annotation data, read counts into new output dir.
    # ---------------------------------------------------------------------------- #
    if args.warm_start_dir:

        if RANK == 0:
            mpi_logging_info('WARM-START: loading data from previous DegNorm run contained in {0}'
                             .format(args.warm_start_dir))
            load_dat = load_from_previous(args.warm_start_dir
                                          , new_dir=output_dir)
            chrom_gene_cov_dict = load_dat['chrom_gene_cov_dict']
            read_count_df = load_dat['read_count_df']
            genes_df = load_dat['genes_df']
            sample_ids = load_dat['sample_ids']

        else:
            chrom_gene_cov_dict = None
            read_count_df = None
            genes_df = None
            sample_ids = None

        # Broadcast the data everyone wants for DegNorm iterations.
        chrom_gene_cov_dict = COMM.bcast(chrom_gene_cov_dict, root=0)
        read_count_df = COMM.bcast(read_count_df, root=0)
        genes_df = COMM.bcast(genes_df, root=0)
        sample_ids = COMM.bcast(sample_ids, root=0)

    # ---------------------------------------------------------------------------- #
    # Path 2: .bam file preprocessing path.
    # Determine intersection of chromosomes across samples from .bam files
    # ---------------------------------------------------------------------------- #
    else:
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
                                          , n_jobs=n_jobs
                                          , verbose=True
                                          , chroms=chroms)
            exon_df = gap.run()
            genes_df = exon_df[['chr', 'gene', 'gene_start', 'gene_end']].drop_duplicates().reset_index(drop=True)

            # take intersection of chromosomes available in genome annotation file and those in the reads data,
            # if for some reason annotation file only contains subset.
            chroms = np.intersect1d(chroms, genes_df.chr.unique()).tolist()
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
        # Load .bam files while simultaneously parsing into coverage arrays. Store chromosome
        # coverage arrays and and compute gene read counts.
        # (Each worker gets ~ same number of .bam files to process.)
        # ---------------------------------------------------------------------------- #
        sample_ids = list()
        cov_files = dict()
        read_count_dict = dict()

        # iterate over node's .bam/.bai files; compute each sample's chromosomes' coverage arrays
        # and save them to .npz files.
        my_files_idx = split_into_chunks(range(n_samples)
                                         , n=SIZE)
        for idx in my_files_idx[RANK]:
            mpi_logging_info('Loading RNA-seq data file {0} -- {1}/{2}'
                             .format(args.bam_files[idx], idx + 1, n_samples))

            reader = BamReadsProcessor(bam_file=args.bam_files[idx]
                                       , index_file=args.bai_files[idx]
                                       , chroms=chroms
                                       , n_jobs=n_jobs
                                       , output_dir=output_dir
                                       , verbose=True)

            sample_id = reader.sample_id
            sample_ids.append(sample_id)
            cov_files[sample_id], read_count_dict[sample_id] = reader.coverage_read_counts(genes_df)

        del reader
        gc.collect()

        # everyone gives sample_ids, cov_file map, and read_counts to master.
        sample_ids = COMM.gather(sample_ids, root=0)
        cov_files = COMM.gather(cov_files, root=0)
        read_count_dict = COMM.gather(read_count_dict, root=0)

        # ---------------------------------------------------------------------------- #
        # Master to merge per-sample gene read count matrices:
        # obtain read count DataFrame containing X, an n (genes) x p (samples) matrix.
        # Master also to order cov_files according to the order of sample_ids.
        # ---------------------------------------------------------------------------- #
        if RANK == 0:

            # collapse sample_ids and read_count_dict from list of lists to something useful.
            sample_ids = flatten_2d(sample_ids).tolist()
            cov_files_unordered = {k: v for d in cov_files for k, v in d.items()}
            read_count_dict = {k: v for d in read_count_dict for k, v in d.items()}

            # order the cov_files by sample_id.
            cov_files = OrderedDict()

            for sample_id in sample_ids:
                cov_files[sample_id] = cov_files_unordered.get(sample_id)

            # Celebrate successful parsing of .bam files.
            mpi_logging_info('Successfully processed chromosome read coverage and gene read counts for all experiments.'
                             .format(len(sample_ids)))
            mpi_logging_info('RNA-seq sample identifiers: \n\t' + ', '.join(sample_ids))

            read_count_df = read_count_dict[sample_ids[0]]
            read_count_df.rename(columns={'read_count': sample_ids[0]}, inplace=True)

            for sample_id in sample_ids[1:]:
                read_count_df = pd.merge(read_count_df
                                         , read_count_dict[sample_id].rename(columns={'read_count': sample_id})
                                         , how='inner'
                                         , on=['chr', 'gene'])

            mpi_logging_info('Successfully merged sample read counts ({0} samples) -- shape: {1}'
                             .format(len(sample_ids), read_count_df.shape))

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
            mpi_logging_info('Saving gene-exon metadata to {0}'.format(exon_output_file))
            exon_df.to_csv(exon_output_file
                           , index=False)

            # save read counts.
            read_count_file = os.path.join(output_dir, 'read_counts.csv')
            mpi_logging_info('Saving original read counts to {0}'.format(read_count_file))
            read_count_df.to_csv(read_count_file
                                 , index=False)

            # ensure chromosomes still reflective of genes to be processed.
            chroms = read_count_df.chr.unique().tolist()

        else:
            # only MASTER has genes_df and read_count_df. Everyone else will get them later on.
            exon_df = None
            chroms = None

        # give everyone the exon data, updated chromosome set,
        # sample_ids, and coverage files for coverage matrix parsing
        exon_df = COMM.bcast(exon_df, root=0)
        chroms = COMM.bcast(chroms, root=0)
        cov_files = COMM.bcast(cov_files, root=0)
        sample_ids = COMM.bcast(sample_ids, root=0)

        # ---------------------------------------------------------------------------- #
        # Slice up genome coverage matrix for each gene according to exon positioning.
        # Run in parallel across workers.
        # ---------------------------------------------------------------------------- #
        chrom_chunks = split_into_chunks(chroms
                                         , n=SIZE)

        # if there are enough chromosome chunks to distribute to workers, disseminate chromosomes,
        # otherwise this worker gets no chromosomes to dice up.
        if RANK < len(chrom_chunks):
            my_chroms = chrom_chunks[RANK]
            gene_cov_dicts = Parallel(n_jobs=min(n_jobs, len(my_chroms))
                           , verbose=0
                           , backend='threading')(delayed(gene_coverage)(
                exon_df=exon_df,
                chrom=chrom,
                coverage_files=cov_files,
                output_dir=output_dir,
                verbose=True) for chrom in my_chroms)

            # convert list of tuples into 2-d dictionary: {chrom: {gene: coverage matrix}}
            chrom_gene_cov_dict = {my_chroms[idx]: gene_cov_dicts[idx] for idx in range(len(my_chroms))}

            del gene_cov_dicts, exon_df
            gc.collect()

            # workers give their coverage matrix data to master.
            if RANK > 0:
                COMM.send(chrom_gene_cov_dict
                          , dest=0
                          , tag=66)

        # host assembles chromosome gene coverage matrix dictionary.
        if RANK == 0:
            if len(chrom_chunks) > 1:
                for worker_id in range(1, len(chrom_chunks)):
                    chrom_gene_cov_dict.update(COMM.recv(source=worker_id
                                                         , tag=66))

    # ---------------------------------------------------------------------------- #
    # PATHS MERGE: warm-start and input-file paths meet.
    # Process gene coverage matrices prior to running NMF.
    # ---------------------------------------------------------------------------- #
    if RANK == 0:
        # initialize an OrderedDict to store coverage matrices that pass DegNorm inclusion
        # criteria in order.
        gene_cov_dict = OrderedDict()

        # Determine for which genes to run DegNorm, and for genes where we will run DegNorm,
        # which transcript regions to filter out prior to running DegNorm.
        mpi_logging_info('Determining genes to include in DegNorm coverage curve approximation.')
        delete_idx = list()

        for i in range(genes_df.shape[0]):
            chrom = genes_df.chr.iloc[i]
            gene = genes_df.gene.iloc[i]

            # extract gene's p x Li coverage matrix. We got chrom_gene_cov_dict from prior gather.
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
        mpi_logging_info('DegNorm will run on {0} genes, downsampling rate = 1 / {1}, {2} baseline selection.'
                         .format(len(gene_cov_dict), args.downsample_rate, 'without' if args.skip_baseline_selection else 'with'))

    else:
        gene_cov_dict = None
        read_count_df = None

    # broadcast per-gene coverage data and read counts.
    gene_cov_dict = COMM.bcast(gene_cov_dict, root=0)
    read_count_df = COMM.bcast(read_count_df, root=0)

    # let's all wait up before going into DegNorm iterations.
    COMM.Barrier()

    # ---------------------------------------------------------------------------- #
    # Run NMF.
    # ---------------------------------------------------------------------------- #
    if RANK == 0:
        mpi_logging_info('Begin executing DegNorm iterations...')

    nmfoa_output = run_gene_nmfoa_mpi(COMM
                                      , cov_dat=gene_cov_dict
                                      , reads_dat=read_count_df[sample_ids].values.astype(np.float_)
                                      , degnorm_iter=args.iter
                                      , nmf_iter=args.nmf_iter
                                      , downsample_rate=args.downsample_rate
                                      , n_jobs=n_jobs
                                      , skip_baseline_selection=args.skip_baseline_selection)

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

        if args.plot_genes:
            plot_genes = np.intersect1d(args.plot_genes, genes)
            plot_genes = split_into_chunks(plot_genes
                                           , n=SIZE)

            # if --plot-genes < number of workers...
            if len(plot_genes) < SIZE:
                plot_genes += [[None]] * (SIZE - len(plot_genes))

    else:
        plot_genes = [[None]] * SIZE

    # master disseminates genes for plotting.
    plot_genes = COMM.scatter(plot_genes
                              , root=0)

    # ---------------------------------------------------------------------------- #
    # Generate coverage curve plots (only if --plot-genes specified)
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
                       , 'genes': genes}

        render_report(data_dir=output_dir
                      , degnorm_data=degnorm_dat
                      , bam_files=args.bam_files if not args.warm_start_dir else [args.warm_start_dir]
                      , sample_ids=sample_ids
                      , top_n_genes=5
                      , output_dir=output_dir)

        mpi_logging_info('DegNorm pipeline complete! Exiting...')

    COMM.Barrier()
    sys.exit(0)


if __name__ == "__main__":
    main()