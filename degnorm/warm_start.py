import os
import shutil
import pickle as pkl
from numpy import intersect1d
from pandas import read_csv


def load_from_previous(degnorm_dir, new_dir):
    """
    To skip the costly DegNorm prerocessing work, e.g. finding per-gene coverage matrices,
    load and copy parsed coverage matrices, genome annotation data, read counts
    from the output directory of a previous DegNorm run into a new DegNorm output directory.

    :param degnorm_dir: str path to directory of a prior DegNorm run's output. Must contain
    gene_exon_metadata.csv and read_counts.csv files, in addition to per-chromosome
    gene coverage matrices saved in binary pickle format.
    :param new_dir: str path to new DegNorm run's output directory.
    :return: dictionary with the core data required to run the new DegNorm pipeline:
    - chrom_gene_cov_dict: 2-d dictionary {chrom: {gene: coverage matrix}}
    - read_count_df: pandas.DataFrame with chr, gene, <sample ID> fields containing per-experiment read counts
    - genes_df: pandas.DataFrame with just gene positioning within chromosomes
    - sample_ids: list of str names of sample ID's, in the same order as in read_count_df.
    """

    if not os.path.isdir(new_dir):
        raise IOError('new DegNorm output directory {0} not found.'.format(new_dir))

    exon_file = os.path.join(degnorm_dir, 'gene_exon_metadata.csv')
    read_count_file = os.path.join(degnorm_dir, 'read_counts.csv')

    # copy and load exon and read count data.
    try:
        shutil.copy(exon_file, os.path.join(new_dir, 'gene_exon_metadata.csv'))
        shutil.copy(read_count_file, os.path.join(new_dir, 'read_counts.csv'))

        exon_df = read_csv(exon_file)
        read_count_df = read_csv(read_count_file)

    except FileNotFoundError as e:
        raise e

    genes_df = exon_df[['chr', 'gene', 'gene_start', 'gene_end']].drop_duplicates().reset_index(drop=True)

    # find intersection of reads + gtf genes, keep that subset.
    intersect_genes = intersect1d(genes_df.gene, read_count_df.gene)
    genes_df = genes_df[genes_df.gene.isin(intersect_genes)]
    read_count_df = read_count_df[read_count_df.gene.isin(intersect_genes)]

    # sort transcript and read count data in same order.
    genes_df.sort_values(['chr', 'gene'], inplace=True)
    genes_df.reset_index(inplace=True, drop=True)
    read_count_df.sort_values(['chr', 'gene'], inplace=True)
    read_count_df.reset_index(inplace=True, drop=True)

    # grab required data: sample IDs, gene manifest, chromosomes encompassed in gene manifest.
    sample_ids = read_count_df.columns.tolist()[2:]
    chroms = genes_df.chr.unique().tolist()
    chrom_gene_cov_dict = dict()

    for idx in range(len(chroms)):
        chrom = chroms[idx]

        # make output sub-dir for chromosome.
        os.makedirs(os.path.join(new_dir, chrom))

        # copy coverage matrices file to new output dir.
        cov_file = os.path.join(degnorm_dir, chrom, 'coverage_matrices_{0}.pkl'.format(chrom))
        shutil.copy(cov_file
                    , dst=os.path.join(new_dir, chrom, 'coverage_matrices_{0}.pkl'.format(chrom)))

        # load the coverage matrices for the genes in this chromosome.
        with open(cov_file, 'rb') as f:
            chrom_gene_cov_dict[chrom] = pkl.load(f)

    # merge data into a dict to return.
    output = dict()
    output['chrom_gene_cov_dict'] = chrom_gene_cov_dict
    output['read_count_df'] = read_count_df
    output['genes_df'] = genes_df
    output['sample_ids'] = sample_ids

    return output

