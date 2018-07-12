import numpy as np


def relative_gene_sample_coverage(sample_df, rng):
    """
    Obtain per-base coverage for one gene for one experiment relative to the starting position
    of the gene on a chromosome.

    :param sample_df: pandas.DataFrame with one RNA-seq sample's genomecov output for a particular
    chromosome
    :param rng: list of two integers denoting start and end positions of a gene on a chromosome
    :return: np.array of read coverage counts
    """
    # outline a gene's location on chromosome including noncoding regions.
    gene_start = np.min(rng)
    gene_end = np.max(rng)

    # TODO: figure out if geneomecov outputs [start, end) or (start, end] RLE output
    gene_reads_df = sample_df[(sample_df['start'] > gene_start) & (sample_df['end'] <= gene_end)]

    # extract RLE vectors from genomecov output, make starts/ends on genome relative
    # to start of the gene.
    start_vec = gene_reads_df['start'].values - gene_start
    end_vec = gene_reads_df['end'].values - gene_start
    counts = gene_reads_df['count'].values
    n = len(start_vec)

    # initialize gene coverage vector.
    coverage_vec = np.zeros(gene_end - gene_start)

    # compute coverage over exon + intron-representative array.
    # Note: this array will need to be spliced in order to represent solely the coding regions.
    for i in range(n):
        coverage_vec[start_vec[i]:end_vec[i]] += counts[i]

    return coverage_vec
