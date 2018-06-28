import pandas as pd


def group_by_genes(df, chrom=True):
    if chrom:
        return df.groupby(['chr', 'gene'])
    else:
        return df.groupby('gene')


def remove_multi_chrom_genes(df):
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


def group_gene_exons(df):
    """
    Pivot a DataFrame of multiple chromosomes' genes' exon start and end locations
    into a table with one row per (chromosome, gene) and an `exons` field containing
    lists of lists of exon (start, end) locations per gene.

    +----------------+-----------------+----------------------------------------------------+
    |     chr        |       gene      |                     exons                          |
    +================+=================+====================================================+
    |     chr1       |      A3GALT2    | [[33772366, 33773054], [33777652, 33777790], [..   |
    +----------------+-----------------+----------------------------------------------------+
    |     chr1       |      AADACL3    | [[12776117, 12776347], [12776117, 12776347], [...  |
    +----------------+-----------------+----------------------------------------------------+

    :param df: pandas.DataFrame containing `start` and `end` columns, e.g.
    from a .bed file
    :return: pandas.DataFrame
    """
    gb = group_by_genes(df)
    groups = list(gb.groups.keys())
    exons = list()

    for group in groups:
        group_df = gb.get_group(group)
        n = group_df.shape[0]
        exons.append([[group_df['start'].iloc[i], group_df['end'].iloc[i]] for i in range(n)])

    region_df = pd.DataFrame({'chr': [group[0] for group in groups]
                              , 'gene': [group[1] for group in groups]
                              , 'exons': exons})
    return region_df


def get_gene_outline(df):
    """
    Append a `range` column to the output of group_gene_exons function. Gives min(start) and max(end)
    location of a gene's exons on a chromosome.

    +-----------------------+-----------------------+
    |   [group_gene_exons]  |       range           |
    +=======================+=======================+
    |   [group_gene_exons]  | [33772366, 33786699]  |
    +-----------------------+-----------------------+
    |   [group_gene_exons]  | [12776117, 12788726]  |
    +-----------------------+-----------------------+

    :param df: pandas.DataFrame containing `start` and `end` columns, e.g.
    from a .bed file
    :return: pandas.DataFrame
    """
    gb = group_by_genes(df)
    gb = gb.apply(lambda x: [min(x.start), max(x.end)])
    df_wide = gb.reset_index()
    df_wide.rename(columns={0: 'range'}
                   , inplace=True)
    df_wide = df_wide.merge(group_gene_exons(df)
                            , on=['chr', 'gene'])

    return df_wide