import pkg_resources
import subprocess
from degnorm.visualizations import *
from degnorm.utils import find_software
from pandas import DataFrame
from jinja2 import Environment, FileSystemLoader


def render_report(data_dir, genenmfoa, input_files,
                  sample_ids, top_n_genes=5, output_dir='.'):
    """
    Render an .html summary report regarding DegNorm output and save it.

    :param data_dir: str location of DegNorm output directory containing coverage matrix plots
    generated by degnorm.visualizations.save_chrom_coverage
    :param genenmfoa: GeneNMFOA object - fitted and transformed model
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
               , 'DegNorm iterations': [genenmfoa.degnorm_iter]
               , 'Downsample rate': ['1/{0}'.format(str(genenmfoa.downsample_rate))]
               , 'Number of input genes': [genenmfoa.n_genes]
               , 'Baseline selection-eligible genes': [np.sum(genenmfoa.use_baseline_selection)]}
    degnorm_info_df = DataFrame(run_dat).transpose()

    # ---------------------------------------------------------------------------- #
    # Render (overlayed) distributions of DI scores (one histogram per sample)
    # ---------------------------------------------------------------------------- #

    # figure out max DI score density.
    y_hist_max = np.ceil(np.max(np.apply_along_axis(lambda x: np.max(plt.histogram(x, density=True)[0])
                                                    , axis=0
                                                    , arr=genenmfoa.rho)))

    fig = plt.figure(figsize=[10, 6])
    fig.suptitle('Degradation index scores, by sample')

    with sns.axes_style('darkgrid'):
        if genenmfoa.rho.shape[0] > 1:

            for i in range(genenmfoa.p):
                sns.distplot(genenmfoa.rho[:, i], label=sample_ids[i])

            plt.xlim(-0.05, 1.05)
            plt.ylim(0., y_hist_max)
            fig.legend(loc='upper right')

        else:
            sns.barplot(sample_ids, y=genenmfoa.rho[0])

        sample_di_dist_plot = os.path.abspath(os.path.join(report_dir, 'di_dists_samples.png'))
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

        ax1 = plt.subplot(gs[0])
        sns.boxplot(data=di_sample_means, orient='h')
        ax1.set_title('Samples')

        ax2 = plt.subplot(gs[1])
        sns.boxplot(data=di_gene_means, orient='h')
        ax2.set_title('Genes')

        for ax in [ax1, ax2]:
            ax.set_yticklabels('')

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    mean_di_dist_plot = os.path.abspath(os.path.join(report_dir, 'di_mean_dists.png'))
    fig.savefig(mean_di_dist_plot
                , dpi=200)

    # ---------------------------------------------------------------------------- #
    # Generate and save coverage plots for genes
    # with top_n_genes highest and lowest average DI scores over samples.
    # ---------------------------------------------------------------------------- #
    avg_dis = genenmfoa.rho.mean(axis=1)
    n_genes = min(top_n_genes, genenmfoa.n_genes)
    hi_di_idx = avg_dis.argsort()[::-1][0:n_genes]
    lo_di_idx = avg_dis.argsort()[0:n_genes]

    hi_di_genes = [genenmfoa.genes[hi_di_idx[i]] for i in range(n_genes)]
    lo_di_genes = [genenmfoa.genes[lo_di_idx[i]] for i in range(n_genes)]

    hi_di_imgs = get_gene_coverage(hi_di_genes
                                   , data_dir=data_dir
                                   , figsize=[10, 6]
                                   , save=True)

    lo_di_imgs = get_gene_coverage(lo_di_genes
                                   , data_dir=data_dir
                                   , figsize=[10, 6]
                                   , save=True)

    # ---------------------------------------------------------------------------- #
    # Find report template and render.
    # ---------------------------------------------------------------------------- #

    # load report template from degnorm package resources.
    resources_dir = pkg_resources.resource_filename('degnorm', 'resources')
    env = Environment(loader=FileSystemLoader(resources_dir))
    template = env.get_template('degnorm_report.html')

    # organize template variables.
    template_vars = {'css_file': os.path.join(resources_dir, 'styles.css')
                     , 'file_info_table': files_df.to_html(index=False
                                                           , bold_rows=False)
                     , 'output_dir': os.path.abspath(output_dir)
                     , 'degnorm_info_table': degnorm_info_df.to_html(index=True
                                                                     , bold_rows=True
                                                                     , header=False)
                     , 'sample_di_dist_plot': sample_di_dist_plot
                     , 'mean_di_dist_plot': mean_di_dist_plot
                     , 'top_n_genes': top_n_genes
                     , 'hi_di_plots': [os.path.abspath(hi_di_imgs[i]) for i in range(n_genes)]
                     , 'lo_di_plots': [os.path.abspath(lo_di_imgs[i]) for i in range(n_genes)]}

    # render report and save.
    html_out = template.render(template_vars)
    html_filename = os.path.join(report_dir, 'degnorm_summary.html')
    with open(html_filename, 'w') as f:
        f.write(html_out)

    # if pandoc is installed and available in path, swap out .html report for a .pdf
    pandoc_avail = find_software('pandoc')
    if pandoc_avail:
        out = subprocess.run(['pandoc {0} -o {1}'.format(html_filename, html_filename.replace('.html', '.pdf'))]
                             , shell=True)

        os.remove(html_filename)