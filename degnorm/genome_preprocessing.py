import itertools
from degnorm.gene_groupby import group_by_genes
from degnorm.utils import *


def remove_multichrom_genes(df):
    """
    Remove from a gene annotation file the genes that show up in multiple chromosomes.


    :param df: pandas.DataFrame version of a gene annotation file that has been converted to .bed format, e.g.

    +----------------+-----------------+-----------------+-----------------+
    |     chr        |       start     |       end       |       gene      |
    +================+=================+=================+=================+
    |     chrI       |      11873      |      12227      |      DDX11L1    |
    +----------------+-----------------+-----------------+-----------------+
    |     chr6       |      17232      |      17368      |       WASH7P    |
    +----------------+-----------------+-----------------+-----------------+

    :return: df, but subsetted to genes that show up in exclusively one chromosome.
    """
    per_chrom_gene_cts = group_by_genes(df, chrom=False).chr.nunique()
    rm_genes = per_chrom_gene_cts[per_chrom_gene_cts > 1].index.tolist()

    return df[~df['gene'].isin(rm_genes)]


def par_apply_exon_overlap(e_df, g_df, chrom):
    """
    Helper function for multiprocessing pool.apply in remove_multigene_exons. Identifies
    exons that overlap with multiple genes.
    """

    # subset range DataFrames to relevant chromosome.
    e_sub_df = e_df[e_df['chr'] == chrom]
    e_sub_df.reset_index(inplace=True, drop=True)
    g_sub_df = g_df[g_df['chr'] == chrom]
    g_sub_df.reset_index(inplace=True, drop=True)

    rm_exon_list = list()

    # iterate over exons contained in this chromosome.
    for e_idx in range(e_sub_df.shape[0]):

        # for particular exon, get get exon start/end locations and exon_id
        e_start, e_end, e_id = e_sub_df.exon_start.iloc[e_idx], e_sub_df.exon_end.iloc[e_idx], \
                               e_sub_df.exon_id.iloc[e_idx]

        # subset to genes that either start or end between exon start/end
        gene_match_df = g_sub_df[(g_sub_df['gene_start'].between(e_start, e_end)) |
                                          (g_sub_df['gene_end'].between(e_start, e_end))]
        if gene_match_df.shape[0] > 1:
            rm_exon_list.append(e_id)

    return rm_exon_list


def remove_multigene_exons(e_df, g_df, n_jobs=max_cpu()):
    """
    Identify exons that span multiple genes and remove them from a population of exons
    over chromosomes in parallel.

    :param e_df: pandas.DataFrame containing details on the population of exons - namely,
    the output from the gene_groupby.exon_outline function (must have `chr`, `gene`, `exon_start`, `exon_end` cols)
    :param g_df: pandas.DataFrame containing details on the population of genes - namely,
    the output from the gene_groupby.get_gene_outline function (must have `chr`, `gene`, `gene_start`, `gene_end` cols)
    :return: pandas.DataFrame subset of exon_range_df corresponding to exons that overlap with exactly a single gene's
    position on the genome
    """
    if 'exon_id' not in e_df.columns:
        e_df['exon_id'] = range(e_df.shape[0])

    chroms = g_df.chr.unique()

    # Iterating over chromosomes in parallel.
    p = mp.Pool(processes=n_jobs)
    rm_exons = [p.apply_async(par_apply_exon_overlap, args=(e_df, g_df, chrom)) for chrom in chroms]
    rm_exons = [x.get() for x in rm_exons]
    p.close()

    # collapse 2-d list of lists into 1-d list
    rm_exons = list(itertools.chain.from_iterable(rm_exons))

    return e_df[~e_df.exon_id.isin(rm_exons)]