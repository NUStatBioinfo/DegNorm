from platform import platform
from re import search
import matplotlib
# if search('Linux', platform()):
matplotlib.use('agg')
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
import os
import pickle as pkl


def plot_gene_coverage(ke, f, x_exon, gene
                       , chrom, sample_ids=None, **kwargs):
    """
    Plot a gene's DegNorm-estimated and original coverage matrices.

    :param ke: estimated coverage matrix
    :param f: original coverage matrix
    :param x_exon: numpy array with two columns, the first - start positions of exons comprising the gene,
    the second - end positions of the exons started in the first column.
    :param gene: str name of gene whose coverage is being plotted
    :param chrom: str name of chromosome gene is on
    :param sample_ids: list of str names of samples corresponding to coverage curves
    :param kwargs: keyword arguments to pass to matplotlib.pylab.figure (e.g. figsize)
    :return: matplotlib.figure.Figure
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

        # upper-left subplot: estimated coverage curves
        ax1 = plt.subplot(gs[0])
        for i in range(ke.shape[0]):
            ax1.plot(ke[i, :], label=sample_ids[i])

        ax1.set_title('Normalized')

        # upper-right subplot: estimated coverage curves
        ax2 = plt.subplot(gs[1])
        for i in range(ke.shape[0]):
            ax2.plot(f[i, :], label=sample_ids[i])

        ax2.set_title('Original')
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

    plt.figlegend(handles, labels, loc='lower center', ncol=len(labels))
    fig.tight_layout(rect=[0, 0.07, 1, 0.95])
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

    with open(coverage_file, 'rb') as f:
        orig_dat = pkl.load(f)

    with open(estimates_file, 'rb') as f:
        est_dat = pkl.load(f)

    for gene in est_dat:
        tmp = exon_df[exon_df.gene == gene]
        fig = plot_gene_coverage(est_dat[gene]
                                 , f=orig_dat[gene]
                                 , x_exon=tmp[['start', 'end']].values
                                 , gene=gene
                                 , chrom=tmp.chr.iloc[0]
                                 , figsize=figsize)

        fig.savefig(os.path.join(output_dir, '{0}_coverage.png').format(gene)
                  , dpi=150)