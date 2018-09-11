import matplotlib
matplotlib.use('agg')
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import matplotlib.pylab as plt
import seaborn as sns
from pandas import read_csv
import numpy as np
import os
import pickle as pkl


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

    # establish exon positioning on chromosome and break points.
    x_exon = x_exon[x_exon[:, 0].argsort()]
    diffs = x_exon[:, 1] - x_exon[:, 0]
    rel_start = x_exon.min()
    rel_end = rel_start + np.sum(diffs)
    start = rel_start
    end = x_exon.max()

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
    if save_dir:

        if not os.path.isdir(os.path.join(save_dir, chrom)):
            os.makedirs(os.path.join(save_dir, chrom))

        fig_path = os.path.join(save_dir, chrom, '{0}_coverage.png'.format(gene))
        fig.savefig(fig_path
                    , dpi=150
                    , bbox_inches='tight')
        plt.close(fig)

        return fig_path

    return fig


def save_chrom_coverage(coverage_file, estimates_file, exon_df,
                        sample_ids, figsize=[10, 6], output_dir='.'):
    """
    Wrapper for plot_gene_coverage: make a coverage plot for every gene in a chromosome.
    """
    for f in [coverage_file, estimates_file]:
        if not os.path.isfile(f):
            raise IOError('Could not find file {0}'.format(f))

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # load original and estimated coverage curve matrices.
    with open(coverage_file, 'rb') as f:
        orig_dat = pkl.load(f)

    with open(estimates_file, 'rb') as f:
        est_dat = pkl.load(f)

    # render and save coverage curve plots one gene at a time.
    for gene in exon_df.gene.unique():

        if (est_dat.get(gene) is not None) and (orig_dat.get(gene) is not None):
            gene_exon_df = exon_df[exon_df.gene == gene]
            out = plot_gene_coverage(est_dat.get(gene)
                                     , f=orig_dat.get(gene)
                                     , x_exon=gene_exon_df[['start', 'end']].values
                                     , gene=gene
                                     , chrom=gene_exon_df.chr.iloc[0]
                                     , sample_ids=sample_ids
                                     , save_dir=output_dir
                                     , figsize=figsize)

        else:
            continue


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
    string filenames of saved plots. If False (default) return list of matplotlib.figure.Figures.
    :return: See save parameter.
    """
    plt.rcParams.update({'figure.max_open_warning': 0})

    # genes should be a list.
    if isinstance(genes, str):
        genes = [genes]

    # quality control: make sure output_dir (and save_dir, if specified) are both real directories.
    if not os.path.isdir(data_dir):
        raise IOError('data_dir {0} is not a directory!'.format(data_dir))

    if not os.path.isfile(os.path.join(data_dir, 'gene_exon_metadata.csv')) \
        or not os.path.isfile(os.path.join(data_dir, 'read_counts.csv')):
        raise ValueError('Missing gene/exon metadata or read count file. Check that {0} is a '
                         ' DegNorm output directory.'.format(data_dir))

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