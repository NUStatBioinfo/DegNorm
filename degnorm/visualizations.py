import matplotlib
matplotlib.use('agg')
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
import os
import pickle as pkl
from pandas import read_csv
from degnorm.utils import flatten_2d
plt.rcParams.update({'figure.max_open_warning': 0})


def get_exon_unions(x):
    """
    Helper function to take union of intersecting exons, for visualization purposes.
    For example, suppose we have two intersecting exons:

    [14563, 14600]
    [14590, 14640]

    They will be reduced into a single exon:
    [14563, 14640]

    :param x: n x 2 numpy array of exon (start, end) pairs
    :return: m x 2 numpy array of unioned exon (start, end) pairs, m <= n.
    """
    if x.shape[0] == 1:
        return x

    # sort exons according to start position.
    x = x[np.argsort(x[:, 0]), :]
    concat_exons = list()

    # determine if exons bleed into any succeeding exons.
    # this happens if exon i's end position >= exon j's start position, j > i.
    for i in range(x.shape[0] - 1):
        bleeding = np.where(x[i, 1] >= x[(i + 1):, 0])[0]
        if len(bleeding) > 0:
            concat_exons.append([i, max(bleeding) + (i + 1)])

    if concat_exons:
        drop_exons = np.unique(flatten_2d(concat_exons))
        union_exons = list()

        # take min(start), max(end) for exons with overlap.
        for j in np.unique([z[1] for z in concat_exons]):
            bleed_start = np.min(np.where([z[1] == j for z in concat_exons])[0])
            start = x[concat_exons[bleed_start][0]][0]
            end = x[j][1]
            union_exons.append([start, end])

        # drop old intersecting exons and stack up newly unioned exons.
        x = np.delete(x, obj=drop_exons, axis=0)
        x = np.vstack([x, np.vstack(union_exons)])
        x = np.unique(x, axis=0)
        x = x[np.argsort(x[:, 0]), :]

    return x


def plot_gene_coverage(ke, f, x_exon, gene
                       , chrom, sample_ids=None
                       , save_dir=None, **kwargs):
    """
    Plot a gene's DegNorm-estimated and original coverage matrices. Option to save image.

    :param ke: estimated coverage matrix
    :param f: original coverage matrix
    :param x_exon: numpy array with two columns, the first - start positions of exons comprising the gene,
    the second - end positions of the exons started in the first column.
    :param gene: str name of gene whose coverage is being plotted
    :param chrom: str name of chromosome gene is on
    :param sample_ids: list of str names of samples corresponding to coverage curves
    :param save_dir: str path to output directory where gene coverage plot should be saved as
    save_dir/<chromosome>/<gene>_coverage.png (optional, default None returns plt.Figure)
    :param kwargs: keyword arguments to pass to matplotlib.pylab.figure (e.g. figsize)
    :return: matplotlib.figure.Figure if save_dir is not specified, otherwise, a filepath to saved .png file
    """
    # quality control.
    if ke.shape != f.shape:
        raise ValueError('ke and f arrays do not have the same shape.')

    if sample_ids:
        if not len(sample_ids) == ke.shape[0]:
            raise ValueError('Number of supplied sample IDs does not match number'
                             'of samples in gene coverage matrix.')
    else:
        sample_ids = ['sample_{0}'.format(i + 1) for i in range(ke.shape[0])]

    # get union of intersecting exons to prevent confusion when plotting exon junctions.
    x_exon = get_exon_unions(x_exon)

    # establish exon positioning on chromosome and junction break points.
    rel_start = x_exon.min()
    start = rel_start

    if x_exon.shape[0] > 1:
        diffs = x_exon[:, 1] - x_exon[:, 0]
        rel_end = rel_start + np.sum(diffs)
        end = x_exon.max()

    else:
        rel_end = x_exon.max()
        end = rel_end

    # before/after plot consists of 4 subplots, establish them.
    fig = plt.figure(**kwargs)
    fig.suptitle('Gene {0} coverage -- chromosome {1}'.format(gene, chrom))
    gs = gridspec.GridSpec(2, 2,
                           width_ratios=[1, 1],
                           height_ratios=[20, 1])

    with sns.axes_style('darkgrid'):

        # upper-left subplot: original coverage curves
        ax1 = plt.subplot(gs[0])
        for i in range(f.shape[0]):
            ax1.plot(f[i, :], label=sample_ids[i])

        ax1.set_title('Original')

        # upper-right subplot: DegNorm-estimated coverage curves
        ax2 = plt.subplot(gs[1])
        for i in range(ke.shape[0]):
            ax2.plot(ke[i, :], label=sample_ids[i])

        ax2.set_title('Normalized')
        handles, labels = ax2.get_legend_handles_labels()

        for ax in [ax1, ax2]:
            ax.margins(x=0)

        # lower-left subplot: exon positioning
        ax3 = plt.subplot(gs[2])
        ax3.set_xlim(start, end)
        ax3.add_patch(Rectangle((start, 0)
                                , width=start + (rel_end - rel_start)
                                , height=1, fill=True, facecolor='red', lw=1))

        # lower-right subplot: exon positioning
        ax4 = plt.subplot(gs[3])
        ax4.set_xlim(start, end)
        ax4.add_patch(Rectangle((start, 0)
                                , width=start + (rel_end - rel_start)
                                , height=1, fill=True, facecolor='red', lw=1))

        # exon splicing y-axis has no meaning, so remove it.
        ax3.get_yaxis().set_visible(False)
        ax4.get_yaxis().set_visible(False)

        # only include start + end position of gene transcript in exon splicing map.
        for ax in [ax3, ax4]:
            ax.get_yaxis().set_visible(False)
            ax.set_xticks([start, end])
            ax.set_xticklabels([str(start), str(end)])

            for i in range(x_exon.shape[0] - 1):
                ax.axvline(x=x_exon[i, 1]
                           , ymin=0
                           , ymax=1
                           , color='w'
                           , lw=2)

    # configure legend: determine if labels expand vertically or horizontally.
    ncol = len(labels) if len(labels) < 6 else 1
    loc = 'upper right' if ncol == 1 else 'lower center'
    bbox_to_anchor = (1.1, 0.85) if ncol == 1 else None

    plt.figlegend(handles
                  , labels
                  , loc=loc
                  , title='Sample'
                  , ncol=ncol
                  , bbox_to_anchor=bbox_to_anchor)

    fig.tight_layout(rect=[0, 0.07, 1, 0.95])

    # if save directory specified, save coverage plots to file save_dir/chrom/<gene>_coverage.png and close figure.
    if not save_dir:
        return fig

    else:
        if not os.path.isdir(os.path.join(save_dir, chrom)):
            os.makedirs(os.path.join(save_dir, chrom))

        fig_path = os.path.abspath(os.path.join(save_dir, chrom, '{0}_coverage.png'.format(gene)))
        fig.savefig(fig_path
                    , dpi=150
                    , bbox_inches='tight')
        plt.close(fig)

        return fig_path


def check_for_files(data_dir, file_names):
    """
    Check that certain files exist within a DegNorm output directory.

    :param data_dir: str path to DegNorm pipeline run output directory
    :param file_names: str or list of str files output from DegNorm pipeline, contained in data_dir
    """
    if not os.path.isdir(data_dir):
        raise IOError('data_dir {0} is not a directory!'.format(data_dir))

    if not isinstance(file_names, list):
        file_names = [file_names]

    are_files = [os.path.isfile(os.path.join(data_dir, f_name)) for f_name in file_names]
    if not all(are_files):
        raise ValueError('Missing requested files. Check that {0} is a '
                         ' DegNorm output directory.'.format(data_dir))


def load_di_scores(data_dir, order=False):
    """
    Load degradation index (DI) scores into a pandas.DataFrame from a .csv file
    in a DegNorm pipeline run output directory. The index column is `gene`, and the
    genes will be ordered alphabetically.

    :param data_dir: str path to DegNorm pipeline run output directory
    :param order: Boolean should samples be ordered (increasing order) according to mean DI score?
    :return: pandas.DataFrame with `gene` as the index, with `chr`, and sample ID columns.
    """
    # load and format DI score data.
    di_file = 'degradation_index_scores.csv'
    out = check_for_files(data_dir
                          , file_names=di_file)
    rho_df = read_csv(os.path.join(data_dir, di_file)
                      , index_col='gene')

    # organize genes in alphabetical order.
    genes = rho_df.index.values
    genes.sort()
    rho_df = rho_df.loc[genes]

    # order sample ids in ascending order of mean DI score.
    ordered_sample_means = rho_df.mean(axis=0).sort_values()
    ordered_sample_ids = ordered_sample_means.index.tolist()

    return rho_df if not order else rho_df[['chr'] + ordered_sample_ids]


def get_di_heatmap(data_dir, save_dir=None, figsize=[10, 8]):
    """
    Generate an n x p (genes x samples) heatmap of degradation index (DI) scores, in ascending order of mean sample
    DI score. Red -> higher DI score, blue -> lower DI score.

    :param data_dir: str path to DegNorm pipeline run output directory
    :param figsize: list or 2-tuple of int [width, height]
    :param save_dir: str directory to save heatmap plot, and then return str filepath of saved plots.
    If None (default), return a matplotlib.figure.Figure.
    :return: See save parameter.
    """
    # load up (ordered) DI scores, drop chromosome column.
    rho_df = load_di_scores(data_dir
                            , order=True)
    rho_df.drop('chr'
                , axis=1
                , inplace=True)

    # make DI score heatmap.
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    fig.suptitle('DI score heatmap')

    sns.heatmap(rho_df
                , cmap='bwr'
                , cbar_kws={"shrink": .5})

    ax.set_xticklabels(ax.get_xticklabels()
                       , rotation=45)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    if save_dir:
        save_path = os.path.abspath(os.path.join(save_dir, 'di_heatmap.png'))
        fig.savefig(save_path
                    , dpi=200)
        plt.close(fig)

        return save_path

    return fig


def get_di_correlation(data_dir, save_dir=None, figsize=[8, 6]):
    """
    Generate a p x p (sample-wise) heatmap of the between-sample correlation matrix.

    :param data_dir: str path to DegNorm pipeline run output directory
    :param figsize: list or 2-tuple of int [width, height]
    :param save_dir: str directory to save heatmap plot, and then return str filepath of saved plots.
    If None (default), return a matplotlib.figure.Figure.
    :return: See save parameter.
    """
    # load up (ordered) DI scores, drop chromosome column.
    rho_df = load_di_scores(data_dir
                            , order=True)
    rho_df.drop('chr'
                , axis=1
                , inplace=True)

    # make DI score correlation matrix heatmap.
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    fig.suptitle('DI score correlation')

    corr = rho_df.corr()
    sns.heatmap(corr
                , xticklabels=corr.columns.values
                , yticklabels=corr.columns.values
                , cmap='YlGnBu'
                , cbar_kws={"shrink": .5})
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    if save_dir:
        save_path = os.path.abspath(os.path.join(save_dir, 'di_correlation.png'))
        fig.savefig(save_path
                    , dpi=200)
        plt.close(fig)

        return save_path

    return fig


def get_di_boxplots(data_dir, save_dir=None, figsize=[12, 8]):
    """
    Generate DI score boxplots (one for each of the p samples).

    :param data_dir: str path to DegNorm pipeline run output directory
    :param figsize: list or 2-tuple of int [width, height]
    :param save_dir: str directory to save heatmap plot, and then return str filepath of saved plots.
    If None (default), return a matplotlib.figure.Figure.
    :return: See save parameter.
    """
    # load up (ordered) DI scores, drop chromosome column.
    rho_df = load_di_scores(data_dir
                            , order=True)
    rho_df.drop('chr'
                , axis=1
                , inplace=True)

    rho_long_df = rho_df.melt(var_name='sample ID'
                              , value_name='DI score')

    # make DI score boxplots
    with sns.axes_style('darkgrid'):
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        fig.suptitle('DI scores')
        sns.boxplot(x='sample ID'
                    , y='DI score'
                    , data=rho_long_df)

        # rotate sample ID labels.
        ax.set_xticklabels(ax.get_xticklabels()
                           , rotation=30)
        ax.set_xlabel('')
        fig.tight_layout(rect=[0, 0, 1, 0.95])

        if save_dir:
            save_path = os.path.abspath(os.path.join(save_dir, 'di_boxplots.png'))
            fig.savefig(save_path
                        , dpi=200)
            plt.close(fig)

            return save_path

        return fig


def get_gene_coverage(genes, data_dir, figsize=[10, 6], save=False):
    """
    Generate gene coverage plots on demand from DegNorm output directory.

    By default, returns a list of matplotlib.figure.Figures, but save=True will
    cause gene coverage plots to be saved into corresponding chromosome directory and
    filepaths of each saved image will be returned.

    :param genes: str or list of str, gene names (case insensitive)
    :param data_dir: str path to DegNorm pipeline run output directory
    :param figsize: [width (int), height (int)] dimensions of coverage curve plots.
    :param save: Bool if True save each plot to <chromosome name>/<gene name>_coverage.png and return
    str filepaths of saved plots. If False (default) return a matplotlib.figure.Figures.
    :return: See save parameter.
    """
    # genes should be a list.
    if isinstance(genes, str):
        genes = [genes]

    # check that data_dir contains the files associated with a DegNorm output dir.
    out = check_for_files(data_dir
                         , file_names=['gene_exon_metadata.csv', 'read_counts.csv', 'degradation_index_scores.csv'])

    # read in required data: exon positioning data and sample ID's (saved in DI score data).
    exon_df = read_csv(os.path.join(data_dir, 'gene_exon_metadata.csv'))
    with open(os.path.join(data_dir, 'degradation_index_scores.csv'), 'r') as di:
        sample_ids = di.readline().strip().split(',')[2:]

    # make genes case-insensitive: cast to uppercase
    genes = [x.upper() for x in genes]
    exon_df.gene = exon_df.gene.apply(lambda x: x.upper())

    # check that all input genes are available in gene/exon metadata.
    avail_genes = exon_df.gene.unique()
    gene_diff = list(set(genes) - set(avail_genes))

    # error out if some genes are not available.
    if gene_diff:
        raise ValueError('Genes {0} were not found in DegNorm output. Check that they were run through pipeline.'
                         .format(', '.join(gene_diff)))

    # subset exon data those requested.
    exon_df = exon_df[exon_df.gene.isin(genes)]
    chroms = exon_df.chr.unique()
    figs = list()

    # iterate over unique chromosomes corresponding to genes requested.
    for chrom in chroms:

        orig_file = os.path.join(data_dir, chrom, 'coverage_matrices_{0}.pkl'.format(chrom))
        ests_file = os.path.join(data_dir, chrom, 'estimated_coverage_matrices_{0}.pkl'.format(chrom))

        # load the payload: coverage curve matrix and DegNorm-approximated coverage matrix dictionaries.
        with open(orig_file, 'rb') as orig, open(ests_file, 'rb') as ests:
            cov_dat = pkl.load(orig)
            ests_dat = pkl.load(ests)

        # cast gene keys of dictionaries to uppercase.
        cov_dat = {k.upper(): v for k, v in cov_dat.items()}
        ests_dat = {k.upper(): v for k, v in ests_dat.items()}
        exon_sub_df = exon_df[exon_df.chr == chrom]

        # determine genes in this chromosome.
        chrom_genes = exon_sub_df.gene.unique()

        # intersect desired chromosome genes with those actually run through DegNorm pipeline.
        nmfoa_genes = list()
        for gene in chrom_genes:
            if (ests_dat.get(gene) is not None) and (cov_dat.get(gene) is not None):
                nmfoa_genes.append(gene)

        chrom_genes = np.intersect1d(chrom_genes, nmfoa_genes)

        # generate plots, save if desired.
        if len(chrom_genes) > 0:

            for gene in chrom_genes:
                figs.append(plot_gene_coverage(ests_dat.get(gene)
                                               , f=cov_dat.get(gene)
                                               , x_exon=exon_sub_df[exon_sub_df.gene == gene][['start', 'end']].values
                                               , gene=gene
                                               , chrom=chrom
                                               , sample_ids=sample_ids
                                               , save_dir=data_dir if save else None
                                               , figsize=figsize))

    return figs