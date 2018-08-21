import degnorm
from degnorm.visualizations import *
from pandas import DataFrame
from jinja2 import Environment, FileSystemLoader


def render_report(data_dir, genenmfoa, gene_manifest_df,
                  input_files, sample_ids, top_n_genes=5, output_dir='.'):
    """
    Render an .html summary report regarding DegNorm output and save it.

    :param data_dir: str location of DegNorm output directory containing coverage matrix plots
    generated by degnorm.visualizations.save_chrom_coverage
    :param genenmfoa: GeneNMFOA object - fitted and transformed model
    :param gene_manifest_df: pandas.DataFrame containing chromosome<->gene map. Must have 'chr' and 'gene' columns.
    :param input_files: list of str names of input .bam files
    :param sample_ids: list of str names of samples
    :param top_n_genes: int number of top-DI and min-DI genes' coverage plots to render
    :param output_dir: str location where .html report should be saved
    """
    # establish save path for report contents.
    report_dir = os.path.join(output_dir, 'report')
    if not os.path.isdir(report_dir):
        os.makedirs(report_dir)

    # ---------------------------------------------------------------------------- #
    # Construct table of input experiment files and sample IDs.
    # ---------------------------------------------------------------------------- #
    files_df = DataFrame([input_files, sample_ids]).T
    files_df.columns = ['Input file', 'Sample ID']

    # ---------------------------------------------------------------------------- #
    # Construct table with DegNorm runtime variable configuration information.
    # ---------------------------------------------------------------------------- #
    run_dat = {'NMF-OA SVD iterations': [genenmfoa.nmf_iter]
               , 'DegNorm iterations': [genenmfoa.iter]
               , 'Downsampled bases': ['N/A' if not genenmfoa.grid_points else str(genenmfoa.grid_points)]
               , 'Number of genes': [genenmfoa.n_genes]}
    degnorm_info_df = DataFrame(run_dat)

    # ---------------------------------------------------------------------------- #
    # Render (overlayed) distributions of DI scores (one histogram per sample)
    # ---------------------------------------------------------------------------- #
    fig = plt.figure(figsize=[10, 6])
    fig.suptitle('Degradation index scores, by sample')

    with sns.axes_style('darkgrid'):
        for i in range(genenmfoa.p):
            sns.distplot(genenmfoa.rho[:, i], label=sample_ids[i])

    plt.xlim(-0.05, 1.05)
    fig.legend(loc='upper right')

    sample_di_dist_plot = os.path.join(report_dir, 'di_dists_samples.png')
    fig.savefig(sample_di_dist_plot
                , dpi=200)

    # ---------------------------------------------------------------------------- #
    # Box and whiskers: (1) distribution of mean DI score per sample
    # + (2) distribution of mean DI score per gene
    # ---------------------------------------------------------------------------- #
    di_sample_means = genenmfoa.rho.mean(axis=0)
    di_gene_means = genenmfoa.rho.mean(axis=1)

    fig = plt.figure(figsize=[10, 6])
    fig.suptitle('Mean degradation index score distributions: by sample, by gene')
    gs = gridspec.GridSpec(2, 1)
    with sns.axes_style('darkgrid'):

        # upper-left subplot: estimated coverage curves
        ax1 = plt.subplot(gs[0])
        sns.boxplot(data=di_sample_means, orient='h')
        ax1.set_title('Samples')

        # upper-right subplot: estimated coverage curves
        ax2 = plt.subplot(gs[1])
        sns.boxplot(data=di_gene_means, orient='h')
        ax2.set_title('Genes')

        for ax in [ax1, ax2]:
            ax.set_yticklabels('')

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    mean_di_dist_plot = os.path.join(report_dir, 'di_mean_dists.png')
    fig.savefig(mean_di_dist_plot
                , dpi=200)

    # ---------------------------------------------------------------------------- #
    # Find saved coverage plot files corresponding to genes
    # with top_n_genes highest and lowest average DI scores over samples.
    # ---------------------------------------------------------------------------- #
    avg_dis = genenmfoa.rho.mean(axis=1)
    hi_di_idx = avg_dis.argsort()[::-1][0:top_n_genes]
    lo_di_idx = avg_dis.argsort()[0:top_n_genes]
    hi_di_imgs = list()
    lo_di_imgs = list()

    # Gather those top_n_genes' coverage plots saved to data_dir.
    for i in range(min(len(top_n_genes), genenmfoa.n_genes)):
        hi_di_gene = genenmfoa.genes[hi_di_idx[i]]
        lo_di_gene = genenmfoa.genes[lo_di_idx[i]]

        hi_di_gene_chrom = gene_manifest_df[gene_manifest_df.gene == hi_di_gene].chr.iloc[0]
        lo_di_gene_chrom = gene_manifest_df[gene_manifest_df.gene == lo_di_gene].chr.iloc[0]

        hi_di_imgs.append(os.path.join(data_dir, hi_di_gene_chrom, '{0}_coverage.png'))
        lo_di_imgs.append(os.path.join(data_dir, lo_di_gene_chrom, '{0}_coverage.png'))

    # ---------------------------------------------------------------------------- #
    # Find report template and render.
    # ---------------------------------------------------------------------------- #
    # load report template.
    template_dir = os.path.join(os.path.dirname(degnorm.__file__), 'report')
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('degnorm_report.html')

    # organize template variables.
    template_vars = {'file_info_table': files_df.to_html()
                     , 'output_dir': output_dir
                     , 'degnorm_info_table': degnorm_info_df.to_html()
                     , 'sample_di_dist_plot': sample_di_dist_plot
                     , 'mean_di_dist_plot': mean_di_dist_plot
                     , 'top_n_genes': top_n_genes
                     , 'hi_di_plots': hi_di_imgs
                     , 'lo_di_plots': lo_di_imgs}

    # render report and save.
    html_out = template.render(template_vars)
    with open(os.path.join(report_dir, 'degnorm_summary.html'), 'w') as f:
        f.write(html_out)