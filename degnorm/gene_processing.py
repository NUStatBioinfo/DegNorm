import itertools
import pandas as pd
from degnorm.utils import *
from degnorm.loaders import GeneAnnotationLoader


class GeneAnnotationProcessor():

    def __init__(self, annotation_file, n_jobs=max_cpu(),
                 genes=None, chroms=None, verbose=True):
        """
        Genome coverage reader for a single RNA-seq experiment.
        Goal is to assemble a dictionary (chromosome, coverage array) pairs.

        :param annotation_file: str .gtf or .gff file
        :param n_jobs: int number of CPUs to use for determining genome coverage. Default
        is number of CPUs on machine - 1
        :param verbose: bool indicator should progress be written to logger?
        :param genes: str or list of str gene names - subset gene annotation data to specified genes
        :param chroms: str or list of str chromosome names - subset gene annotation data to
         specified chromosomes
        """
        self.filename = annotation_file
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.genes = genes
        self.chroms = chroms
        self.loader = None

        # cast any supplied chromosomes or genes into a list if not already in lists.
        if self.chroms:
            if not isinstance(self.chroms, list):
                self.chroms = [self.chroms]

        if self.genes:
            if not isinstance(self.genes, list):
                self.genes = [self.genes]

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

            exon_df = subset_to_chrom(exon_df, chrom=self.chroms)

        if self.genes:
            if self.verbose:
                logging.info('Subsetting exon data to {0} genes.'
                             .format(len(self.genes)))

            exon_df = exon_df[exon_df.gene.isin(self.genes)]

        if self.verbose:
            logging.info('Successfully loaded exon data -- shape: {0}'.format(exon_df.shape))

        if exon_df.empty:
            raise ValueError('Exon DataFrame is empty! Perhaps your chromosome or gene subsetting'
                             'is not going according to plan.')

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
    def get_gene_outline(df):
        """
        Based on the exons present for a particular gene on a particular chromosome, gives the min(start) and max(end)
        location of all this gene's exons to form a gene "outline."

        +-----------+-------------+-------------------+-------------------+
        |    chr    |    gene     |    gene_start     |    gene_start     |
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

    @staticmethod
    def find_overlapping_exons(exon_df, gene_df, chrom):
        """
        Helper function for multiprocessing pool.apply in remove_multigene_exons. Identifies
        exons within one chromosome that overlap with multiple genes.

        :param exon_df: pandas.DataFrame with fields 'chr', 'gene', 'start', 'end' outlining the transcriptome
        on a chromosome; there are multiple rows (so, multiple exons) per (chr, gene) combination
        :param gene_df: pandas.DataFrame with fields 'chr', 'gene', 'gene_start', 'gene_end', and 'exon_id'
         outlining the positions of genes on a chromosome
        :param chrom: str name of chromosome
        :return: list of exon identifiers for exons that intersect with > 1 genes
        """
        # subset range DataFrames to relevant chromosome.
        exon_chrom_df = subset_to_chrom(exon_df, chrom=chrom, reindex=True)
        exon_chrom_df = exon_chrom_df[['start', 'end', 'exon_id']]
        gene_chrom_df = subset_to_chrom(gene_df, chrom=chrom, reindex=True)

        # store list of exons that need to be removed.
        rm_exon_list = list()

        # iterate over exons contained in this chromosome.
        for exon_idx in range(exon_chrom_df.shape[0]):

            # for particular exon, get get exon start/end locations and exon_id
            exon_start, exon_end, exon_id = exon_chrom_df.iloc[exon_idx].values

            # subset to genes that have overlap with this exon
            gene_match_df = gene_chrom_df[(gene_chrom_df.gene_start.between(exon_start, exon_end)) |
                                          (gene_chrom_df.gene_end.between(exon_start, exon_end)) |
                                          ((gene_chrom_df.gene_start > exon_start) & (gene_chrom_df.gene_end < exon_end))]

            if gene_match_df.shape[0] > 1:
                rm_exon_list.append(exon_id)

        return rm_exon_list

    def remove_multigene_exons(self, exon_df, gene_df):
        """
        Identify exons that span multiple genes and remove them from a population of exons
        over chromosomes in parallel.

        :param exon_df: pandas.DataFrame with fields 'chr', 'gene', 'start', 'end' outlining the transcriptome
        on a chromosome; there are multiple rows (so, multiple exons) per (chr, gene) combination
        :param gene_df: pandas.DataFrame with fields 'chr', 'gene', 'gene_start', 'gene_end', and 'exon_id'
         outlining the positions of genes on a chromosome
        :return: pandas.DataFrame subset of exon_df corresponding to exons that overlap with exactly a single gene's
        position on the genome
        """
        if 'exon_id' not in exon_df.columns:
            exon_df['exon_id'] = range(exon_df.shape[0])

        chroms = gene_df.chr.unique()

        # Iterating over chromosomes in parallel.
        p = mp.Pool(processes=self.n_jobs)
        rm_exons = [p.apply_async(self.find_overlapping_exons, args=(exon_df, gene_df, chrom)) for chrom in chroms]
        p.close()
        rm_exons = [x.get() for x in rm_exons]

        # collapse 2-d list of lists into 1-d list
        rm_exons = list(itertools.chain.from_iterable(rm_exons))

        return exon_df[~exon_df.exon_id.isin(rm_exons)]

    def run(self):
        """
        Main function for GeneAnnotationProcessor. Runs transcriptome annotation processing pipeline:

         1. loads .gtf or .gff file for exon regions
         2. removes genes that occur in multiple chromosomes
         3. eliminates exons that intersect with > 1 genes on a single chromosome

        :param chroms: list of str chromosomes with which to subset data. Use if
        only a subset of chromosomes from genome annotation file will be useful.
        :return: pandas.DataFrame with 'chr', 'gene', 'gene_start', 'gene_end', 'exon_start', 'exon_end' fields
        outlining exons ready to use for reads coverage computations
        """
        if self.verbose:
            logging.info('Loading genome annotation file {0} into pandas.DataFrame'.format(self.filename))

        exon_df = self.load()

        logging.info('Removing genes that occur in multiple chromosomes.')
        exon_df = self.remove_multichrom_genes(exon_df)
        exon_df.drop_duplicates(inplace=True)

        if self.verbose:
            logging.info('Successfully removed multiple-chromosome genes.')
            logging.info('Removing exons that occur in multiple genes.')

        gene_df = self.get_gene_outline(exon_df)
        n_exons = exon_df.shape[0]
        exon_df = self.remove_multigene_exons(exon_df=exon_df
                                              , gene_df=gene_df)

        if self.verbose:
            logging.info('Multiple-chromosome genes were removed.')

        exon_df = exon_df.merge(gene_df, on=['chr', 'gene'])
        exon_df.drop_duplicates(inplace=True)

        if self.verbose:
            logging.info('{0} multiple-gene exons were removed. Final shape -- {1}'
                         .format(n_exons - exon_df.shape[0], exon_df.shape))
            logging.info('Genome annotation file processing pipeline was successful.')

        return exon_df


# if __name__ == '__main__':
#     import os
#     gtf_file = os.path.join(os.getenv('HOME'), 'nu', 'jiping_research', 'data', 'rna_seq', 'genes.gtf')
#     processor = GeneAnnotationProcessor(gtf_file, n_jobs=2, chroms=['chr1', 'chr2'])
#     exon_df = processor.run()