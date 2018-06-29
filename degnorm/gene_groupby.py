import pandas as pd


def group_by_genes(df, chrom=True):
    if chrom:
        return df.groupby(['chr', 'gene'])
    else:
        return df.groupby('gene')


def group_gene_exons(df):
    """
    Pivot a DataFrame of multiple chromosomes' genes' exon start and end locations
    into a table with one row per (chromosome, gene) and an `exons` field containing
    lists of lists of exon (start, end) locations per gene.

    +----------------+-----------------+-----------------------------+----------------------------+
    |     chr        |       gene      |       exon_starts           |       exon_ends            |
    +================+=================+=============================+============================+
    |     chr1       |      A3GALT2    | [33772366, 33777652, ...]   | [33773054,33777790, ... ]  |
    +----------------+-----------------+-----------------------------+----------------------------+
    |     chr1       |      AADACL3    | [12776117, 12776117, ...]   | [12776347, 12776347, ...]  |
    +----------------+-----------------+-----------------------------+----------------------------+

    :param df: pandas.DataFrame containing `start` and `end` columns, e.g.
    from a .bed file
    :return: pandas.DataFrame
    """
    gb = group_by_genes(df)
    groups = list(gb.groups.keys())
    starts = list()
    ends = list()

    for group in groups:
        group_df = gb.get_group(group)
        n = group_df.shape[0]
        starts.append([group_df.start.iloc[i] for i in range(n)])
        ends.append([group_df.end.iloc[i] for i in range(n)])

    region_df = pd.DataFrame({'chr': [group[0] for group in groups]
                              , 'gene': [group[1] for group in groups]
                              , 'exon_starts': starts
                              , 'exon_ends': ends})
    return region_df


def get_gene_outline(df):
    """
    Append a `range` column to the output of group_gene_exons function. Gives min(start) and max(end)
    location of a gene's exons on a chromosome.

    +-----------------------+-------------------+-------------------+
    |   [group_gene_exons]  |    gene_start     |    gene_start     |
    +=======================+===================+===================+
    |   [group_gene_exons]  |     33772366      |     33779539      |
    +-----------------------+-------------------+-------------------+
    |   [group_gene_exons]  |     12776117      |     12784120      |
    +-----------------------+-------------------+-------------------+

    :param df: pandas.DataFrame containing `start` and `end` columns, e.g.
    from a .bed file
    :return: pandas.DataFrame
    """
    gb = group_by_genes(df)
    df_wide = gb.apply(lambda x: pd.Series({'gene_start': min(x.start), 'gene_end': max(x.end)}))
    df_wide.reset_index(inplace=True)
    df_wide = df_wide.merge(group_gene_exons(df)
                            , on=['chr', 'gene'])

    return df_wide


def get_exon_outline(df):
    """
    Trace the outline of individual exons in a genome annotation.

    :param df: pandas.DataFrame, output from get_gene_outline function
    :return: pandas.DataFrame with `chr`, `gene`, `exon_start`, and `exon_end` fields:

    +----------------+-----------------+---------------------+------------------+
    |     chr        |       gene      |     exon_start      |    exon_end      |
    +================+=================+=====================+==================+
    |     chr1       |      RCAN3AS    |      24828837       |     24828850     |
    +----------------+-----------------+---------------------+------------------+
    |     chr1       |      RCC1       |      28832454       |     28832596     |
    +----------------+-----------------+---------------------+------------------+
    """

    exon_df_list = list()
    for row_idx in range(df.shape[0]):
        sub_df = df.iloc[row_idx]

        exon_df_list.append(pd.DataFrame({'chr': sub_df.chr
                                             , 'gene': sub_df.gene
                                             , 'exon_start': sub_df.exon_starts
                                             , 'exon_end': sub_df.exon_ends}))

    exon_outline_df = pd.concat(exon_df_list)
    exon_outline_df.drop_duplicates(inplace=True)
    exon_outline_df['exon_id'] = range(exon_outline_df.shape[0])

    return exon_outline_df
