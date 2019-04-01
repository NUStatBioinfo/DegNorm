import HTSeq
import networkx as nx
import pandas as pd
from degnorm.utils import *
from degnorm.loaders import GeneAnnotationLoader


class GeneAnnotationProcessor():

    def __init__(self, annotation_file, chroms=None, verbose=True):
        """
        Genome coverage reader for a single RNA-seq experiment.
        Goal is to assemble a dictionary (chromosome, coverage array) pairs.

        :param annotation_file: str .gtf or .gff file
        :param verbose: bool indicator should progress be written to logger?
        :param chroms: str or list of str chromosome names - subset gene annotation data to
         specified chromosomes
        """
        self.filename = annotation_file
        self.verbose = verbose
        self.chroms = chroms
        self.loader = None

        # cast any supplied chromosomes or genes into a list if not already in lists.
        if self.chroms:
            if not isinstance(self.chroms, list):
                self.chroms = [self.chroms]

    def load(self):
        """
        Loads exons from a .gtf or .gff file.
        """
        self.loader = GeneAnnotationLoader(self.filename)
        exon_df = self.loader.get_data()

        if self.chroms:
            if self.verbose:
                logging.info('Subsetting exon data to {0} chromosomes:\n'
                             '\t{1}'.format(len(self.chroms), ', '.join(self.chroms)))

            exon_df = subset_to_chrom(exon_df
                                      , chrom=self.chroms)

        if self.verbose:
            logging.info('Successfully loaded exon data -- shape: {0}'.format(exon_df.shape))

        if exon_df.empty:
            raise ValueError('Exon DataFrame is empty!')

        return exon_df

    @staticmethod
    def remove_multichrom_genes(df):
        """
        Remove from a gene annotation file the genes that show up in multiple chromosomes.

        :return: df, but subset to genes that show up in exclusively one chromosome.
        """
        grp = df.groupby(['gene'])
        per_chrom_gene_cts = grp.chr.nunique()
        rm_genes = per_chrom_gene_cts[per_chrom_gene_cts > 1].index.tolist()

        return df[~df.gene.isin(rm_genes)]

    @staticmethod
    def gene_outline(df):
        """
        Based on the exons present for a particular gene on a particular chromosome, gives the min(start) and max(end)
        location of all this gene's exons to form a gene "outline."

        +-----------+-------------+-------------------+-------------------+
        |    chr    |    gene     |    gene_start     |     gene_end      |
        +===========+=============+===================+===================+
        |   chr3    |     WASH7P  |     33772366      |     33779539      |
        +-----------+-------------+-------------------+-------------------+
        |    chr3   |  MIR6859-3  |     12776117      |     12784120      |
        +-----------+-------------+-------------------+-------------------+

        :param df: pandas.DataFrame containing 'gene', 'chr', exon 'start` and 'end' columns, e.g. from a .bed file
        :return: pandas.DataFrame with 'chr', 'gene', 'gene_start' and 'gene_end' columns
        """
        grp = df.groupby(['chr', 'gene'])
        df_wide = grp.apply(lambda x: pd.Series({'gene_start': min(x.start), 'gene_end': max(x.end)}))
        df_wide.reset_index(inplace=True)

        return df_wide

    def run(self):
        """
        Main function for GeneAnnotationProcessor. Runs transcriptome annotation processing pipeline:

         1. loads .gtf or .gff file for exon regions
         2. removes genes that occur in multiple chromosomes
         3. outlines gene start/end positions
         4. removes duplicates.

        :param chroms: list of str chromosomes with which to subset data. Use if
        only a subset of chromosomes from genome annotation file will be useful.
        :return: pandas.DataFrame with 'chr', 'gene', 'gene_start', 'gene_end',
        'start' [exon start], and 'end'[exon end] fields
        outlining exons ready to use for reads coverage computations
        """
        if self.verbose:
            logging.info('Loading genome annotation file {0}...'.format(self.filename))

        exon_df = self.load()

        if self.verbose:
            logging.info('Begin genome annotation file processing.')

        # gene / exon processing
        exon_df = self.remove_multichrom_genes(exon_df)
        exon_df.drop_duplicates(inplace=True)
        gene_df = self.gene_outline(exon_df)
        exon_df = exon_df.merge(gene_df
                                , on=['chr', 'gene'])
        exon_df.drop_duplicates(inplace=True)

        if self.verbose:
            logging.info('Processing successful. Final shape -- {0}'.format(exon_df.shape))

        return exon_df


def get_gene_overlap_structure(gene_df):
    """
    Build gene intersection matrix (a type of adjacency matrix), splitting genes into groups
    of mutually overlapping genes (i.e. groups of genes that are reachable on paths within
    adjacency matrix) and isolated genes that have no overlap with others.

    Example: Let gene_df be the pandas.DataFrame

    +-----------+-------------+-------------------+-------------------+
    |    chr    |    gene     |    gene_start     |     gene_end      |
    +===========+=============+===================+===================+
    |   chr3    |     WASH7P  |     100           |     200           |
    +-----------+-------------+-------------------+-------------------+
    |    chr3   |  MIR6859-3  |     150           |     230           |
    +-----------+-------------+-------------------+-------------------+
    |    chr3   |     RDVC    |     215           |     280           |
    +-----------+-------------+-------------------+-------------------+
    |    chr3   |     EZH2    |     600           |      822          |
    +-----------+-------------+-------------------+-------------------+

    then get_gene_overlap_structure will return
    {'isolated_genes': ['EZH2'], 'overlap_genes': ['WASH7P', 'MIR6895-3', 'RDVC']}


    :param gene_df: pandas.DataFrame containing a chromosome's genes' location data, has at least these columns:
    `gene` (str name of gene), `gene_start` (int leftmost base position of gene), `gene_end` (int rightmost
    base position of gene)
    :return: dict of two elements, 'overlap_genes' which is a list of lists (sublists are groups of overlapping genes,
    e.g. if gene A overlaps gene B and gene B overlaps gene C - no assurance that gene A overlaps gene C - then
    genes A, B, and C form a group). Second element is 'isolated genes', a list of genes that have no overlap
    with others.
    """

    genes = gene_df.gene.values
    n_genes = len(genes)

    # initialize gene overlap adjacency matrix.
    adj_mat = np.zeros([n_genes, n_genes])

    gene_starts = gene_df.gene_start.values
    gene_ends = gene_df.gene_end.values
    gas = HTSeq.GenomicArrayOfSets(['chrom']
                                   , stranded=False)

    # build gene locator gas. Use 0-indexing.
    for i in range(n_genes):
        iv = HTSeq.GenomicInterval('chrom', gene_starts[i] - 1, gene_ends[i], '.')
        gas[iv] += str(i)

    # iterate over genes, find other genes that have overlap.
    # for genes with overlap, store their exon region bounds.
    for i in range(n_genes):
        gene_start, gene_end = gene_starts[i] - 1, gene_ends[i]

        # search gas for overlapping genes.
        gas_intersect = [(st[0], sorted(st[1])) for
                         st in gas[HTSeq.GenomicInterval('chrom', gene_start, gene_end, '.')].steps()]

        # parse results from gas search, identify intersecting genes.
        gene_intersect = list()
        for ii in range(len(gas_intersect)):
            cross_genes = gas_intersect[ii][1]
            if cross_genes:
                gene_intersect.extend(cross_genes)

        # get indices of intersecting genes, update row of adjacency matrix.
        gene_intersect = [int(x) for x in list(set(gene_intersect))]
        if gene_intersect:
            adj_mat[i, gene_intersect] = 1

    # build network graph from adjacency (overlap) matrix.
    graph = nx.from_numpy_matrix(adj_mat)

    # now, parse graph:
    # 1. starting with one gene, find all genes reachable from this gene. This is one
    # adjacency group. If no other genes are reachable, gene is isolated.
    # 2. Find set diff of adjacency groups with remaining genes, start search step 1. with
    # any of these remaining genes.
    # (continue iterating 1. and 2. until all genes have been grouped or labeled isolated.)
    search_gene_idx = 0
    progress = 0
    isolated_genes = list()
    overlap_genes = list()  # overlap_genes will be list of lists (sublists are groups of overlapping genes)
    gene_ids = list(range(n_genes))

    while progress < n_genes:
        # get gene id's for all genes reachable from search_gene_idx
        reachable_gene_idx = list(nx.single_source_shortest_path(graph, search_gene_idx).keys())
        progress += len(reachable_gene_idx)

        # if gene is only reachable to itself -->> no overlap, append name to list of isolated genes.
        if len(reachable_gene_idx) == 1:
            isolated_genes.append(genes[reachable_gene_idx[0]])

        # if gene is in group of overlapping genes -->> store entire this set of overlapping genes.
        else:
            overlap_genes.append(genes[[x for x in reachable_gene_idx]].tolist())

        gene_ids = list(set(gene_ids) - set(reachable_gene_idx))
        if gene_ids:
            search_gene_idx = gene_ids[0]
        else:
            break

    return {'overlap_genes': overlap_genes
            , 'isolated_genes': isolated_genes}


# if __name__ == '__main__':
#     import os
#     import time
#
#     gtf_file = os.path.join(os.getenv('HOME'), 'nu', 'jiping_research', 'data', 'rna_seq', 'genes.gtf')
#     processor = GeneAnnotationProcessor(gtf_file)
#
#     start = time.time()
#     exon_df = processor.run()
#     end = time.time()
#
#     print('human genome .gtf file parsing time: {0}'.format(end - start))
#     print('exon_df.shape: {0}'.format(exon_df.shape))
#     print(exon_df.head())