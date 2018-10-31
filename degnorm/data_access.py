from degnorm.visualizations import *
import pickle as pkl
import gc
import tqdm
from pandas import DataFrame
plt.rcParams.update({'figure.max_open_warning': 0})


class CoverageLoader(object):
    def __init__(self, data_dir):
        """
        Initialize a CoverageLoader object for loading coverage arrays from a DegNorm output directory.

        :param data_dir: str path to a DegNorm output directory
        """
        self.genes = None
        self.chroms = None
        self.sample_ids = None
        self.exon_df = None
        self.cov_dict = dict()

        # check that DegNorm dir exists.
        if not os.path.isdir(data_dir):
            raise NotADirectoryError('{0} is not a directory, cannot find DegNorm output.'.format(data_dir))

        self.data_dir = data_dir

        # check that data_dir contains the files associated with a DegNorm output dir.
        out = check_for_files(self.data_dir
                              , file_names=['gene_exon_metadata.csv', 'read_counts.csv', 'degradation_index_scores.csv'])

    def load(self, genes):
        """
        Load coverage arrays and exon positioning data for a set of genes from a DegNorm output directory,
        data to be used by other functions for saving arrays to .txt or plotting pre/post-Degnorm coverage.

        :param genes: str or list of str, names of genes for which to load up coverage arrays. Can also be
        the string "all" indicating we want coverage for all genes sent through DegNorm.
        """
        # discern between all/one/multiple input genes.
        all_genes = False
        if isinstance(genes, str):
            if genes.lower() == 'all':
                all_genes = True
            else:
                genes = [genes]

        self.genes = genes

        # read in required data: exon positioning data and sample ID's (saved in DI score data).
        self.exon_df = read_csv(os.path.join(self.data_dir, 'gene_exon_metadata.csv'))
        with open(os.path.join(self.data_dir, 'degradation_index_scores.csv'), 'r') as di:
            self.sample_ids = di.readline().strip().split(',')[2:]

        if all_genes:
            self.genes = self.exon_df.gene.unique().tolist()
            print('Loading coverage data for all {0} genes sent through DegNorm pipeline. This could take a minute.'
                  .format(len(self.genes)))

        # make genes case-insensitive: cast to uppercase.
        self.genes = [x.upper() for x in self.genes]
        self.exon_df.gene = self.exon_df.gene.apply(lambda x: x.upper())

        # check that every input gene is available in gene/exon metadata.
        if not all_genes:
            avail_genes = self.exon_df.gene.unique()
            gene_diff = list(set(self.genes) - set(avail_genes))

            # error out if some genes are not available.
            if gene_diff:
                raise ValueError('Genes {0} were not found in DegNorm output. Check that they were run through pipeline.'
                                 .format(', '.join(gene_diff)))

            # subset exon data to requested genes.
            self.exon_df = self.exon_df[self.exon_df.gene.isin(self.genes)]

        # iterate over unique chromosomes corresponding to genes requested.
        for chrom in self.exon_df.chr.unique():

            raw_file = os.path.join(self.data_dir, chrom, 'coverage_matrices_{0}.pkl'.format(chrom))
            ests_file = os.path.join(self.data_dir, chrom, 'estimated_coverage_matrices_{0}.pkl'.format(chrom))

            # load the payload: gene dictionaries with coverage curve numpy arrays.
            with open(raw_file, 'rb') as raw, open(ests_file, 'rb') as est:
                raw_dat = pkl.load(raw)
                est_dat = pkl.load(est)

            # cast gene keys of dictionaries to uppercase.
            raw_dat = {k.upper(): v for k, v in raw_dat.items()}
            est_dat = {k.upper(): v for k, v in est_dat.items()}

            # determine genes in this chromosome.
            chrom_genes = self.exon_df[self.exon_df.chr == chrom].gene.unique()

            # if gene is available in both raw + estimated coverage, append those matrices.
            for gene in chrom_genes:
                raw_cov = raw_dat.get(gene)
                est_cov = est_dat.get(gene)

                if (raw_cov is not None) and (est_cov is not None):
                    self.cov_dict[gene] = dict()
                    self.cov_dict[gene]['raw'] = raw_cov
                    self.cov_dict[gene]['estimate'] = est_cov

        # clean up (coverage dicts are large)
        del raw_dat, est_dat
        gc.collect()


def get_coverage_plots(genes, degnorm_dir, figsize=[10, 6], save_dir=None):
    """
    Generate pre- and post-DegNorm gene coverage plots on demand from a DegNorm output directory.

    By default, returns a list of matplotlib.figure.Figures, but save=True will
    cause gene coverage plots to be saved into corresponding chromosome directory and
    filepaths of each saved image will be returned.

    :param genes: str or list of str, gene names (case insensitive)
    :param degnorm_dir: str path to DegNorm pipeline run output directory
    :param figsize: [width (int), height (int)] dimensions of coverage curve plots.
    :param save_dir: If specified, save each plot to save_dir/<chromosome name>/<gene name>_coverage.png and return
    str filepaths of saved plots. If False (default) return a matplotlib.figure.Figures.
    :return: See save parameter.
    """
    cov_ldr = CoverageLoader(degnorm_dir)
    cov_ldr.load(genes)

    figs = list()

    # instantiate progress bar if making > 100 plots.
    use_pbar = len(cov_ldr.cov_dict) > 100
    if use_pbar:
        pbar_step_size = int(np.floor(len(cov_ldr.cov_dict) / 20))
        pbar = tqdm.tqdm(total=100
                         , leave=False
                         , desc='plotting progress'
                         , unit='%')

    # iterate over loaded genes.
    ctr = 0
    for gene in cov_ldr.cov_dict:

        # extract coverage.
        raw_cov = cov_ldr.cov_dict.get(gene)['raw']
        est_cov = cov_ldr.cov_dict.get(gene)['estimate']

        # extract exon positioning, chromosome name.
        gene_exon_df = cov_ldr.exon_df[cov_ldr.exon_df.gene == gene]
        x_exon = gene_exon_df[['start', 'end']].values
        chrom = gene_exon_df.chr.iloc[0]

        figs.append(plot_gene_coverage(est_cov
                                       , f=raw_cov
                                       , x_exon=x_exon
                                       , gene=gene
                                       , chrom=chrom
                                       , sample_ids=cov_ldr.sample_ids
                                       , save_dir=degnorm_dir if not save_dir else save_dir
                                       , figsize=figsize))
        ctr += 1

        # update progress bar in 5% intervals
        if use_pbar:
            if (ctr % pbar_step_size == 0) and (ctr > 0):
                pbar.update(5)

    # close progress bar if using one.
    if use_pbar:
        pbar.close()

    return figs


def get_coverage_data(genes, degnorm_dir, save_dir=None):
    """
    Access raw and DegNorm-estimated coverage matrices of a set of genes run through the DegNorm pipeline.

    By default, returns two lists:
        - raw coverage pandas.DataFrames
        - DegNorm-estimated coverage pandas.DataFrames

    Each coverage pandas.DataFrame is in *long* format: (base positions x samples)

    +-------------+-----------+---------------+
    |  <sample 1> |   <...>   |  <sample p>   |
    +=============+===========+===============+
    |     0.0     |    ...    |      11.0     |
    +-------------+-----------+---------------+

    :param genes: str or list of str, gene names (case insensitive)
    :param degnorm_dir: str path to DegNorm pipeline run output directory
    :param save_dir: If specified, save raw and estimated coverage arrays to
    save_dir/<chromosome name>/<gene name>_raw_coverage.txt, save_dir/<chromosome name>/<gene name>_estimated_coverage.txt.
    :return: nested dictionary:
    gene (key): dict('raw': raw coverage DataFrame, 'estimate': estimated coverage DataFrame) (value)
    """
    cov_ldr = CoverageLoader(degnorm_dir)
    cov_ldr.load(genes)

    cov_dat_output = dict()

    # instantiate progress bar if making > 100 plots.
    use_pbar = len(cov_ldr.cov_dict) > 100
    if use_pbar:
        pbar_step_size = int(np.floor(len(cov_ldr.cov_dict) / 20))
        pbar = tqdm.tqdm(total=100
                         , leave=False
                         , desc='save progress'
                         , unit='%')

    # iterate over loaded genes.
    ctr = 0
    for gene in cov_ldr.cov_dict:

        # extract coverage.
        raw_cov = cov_ldr.cov_dict[gene]['raw']
        est_cov = cov_ldr.cov_dict[gene]['estimate']

        # store coverage as Li x p DataFrames.
        cov_dat_output[gene] = dict()
        cov_dat_output[gene]['raw'] = DataFrame(raw_cov.T
                                                , columns=cov_ldr.sample_ids)
        cov_dat_output[gene]['estimate'] = DataFrame(est_cov.T
                                                     , columns=cov_ldr.sample_ids)

        # protocol to save raw and estimated coverage to disk...
        if save_dir:

            # extract gene's chromosome name.
            chrom = cov_ldr.exon_df[cov_ldr.exon_df.gene == gene].chr.iloc[0]

            # ensure writability of coverage data: create missing directories.
            if not os.path.isdir(os.path.join(save_dir, chrom)):
                os.makedirs(os.path.join(save_dir, chrom))

            raw_save_file = os.path.join(save_dir, chrom, '{0}_raw_coverage.txt'.format(gene))
            est_save_file = os.path.join(save_dir, chrom, '{0}_estimated_coverage.txt'.format(gene))

            # write coverage.
            cov_dat_output[gene]['raw'].to_csv(raw_save_file
                                               , index=None
                                               , sep=' '
                                               , float_format='%.5f')
            cov_dat_output[gene]['estimate'].to_csv(est_save_file
                                                    , index=None
                                                    , sep=' '
                                                    , float_format='%.5f')

            ctr += 1

            # update progress bar in 5% intervals
            if use_pbar:
                if (ctr % pbar_step_size == 0) and (ctr > 0):
                    pbar.update(5)

    # close progress bar if using one.
    if use_pbar:
        pbar.close()

    return cov_dat_output