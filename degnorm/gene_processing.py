import HTSeq
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
            raise ValueError('Exon DataFrame is empty! Perhaps your chromosome or gene subsetting '
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

    @staticmethod
    def remove_multigene_exons(ex_df):
        """
        Helper function for parallel execution in remove_multigene_exons. Identifies
        exons within one chromosome that overlap with multiple genes and removes them.

        Input: exons with gene where each row is an exon. Chromosome, gene metadata attached
        to each exon. Output: one -> many relationship between genes and exons.

        +-----------+---------+------------+------------+
        |    chr    |  start  |    end     |    gene    |
        +===========+=========+============+============+
        |   chr3    |  99909  |    12003   |    WASH7P  |
        +-----------+---------+------------+------------+
        |   chr3    |  14560  |    14612   |    WASH7P  |
        +-----------+---------+------------+------------+

        :param ex_df: pandas.DataFrame with fields 'chr', 'gene', 'start', 'end' outlining exons in a
        chromosome; genes with >1 exon are represented across multiple rows.
        :param gene_df: pandas.DataFrame with fields 'chr', 'gene', 'gene_start', 'gene_end', and 'exon_id'
         outlining the positions of genes on a chromosome
        :param chrom: str name of chromosome
        :return: pandas.DataFrame with as many or fewer exons as original exon_df. Now each exon only
        corresponds to one gene.
        """
        # ensure column order.
        ex_df = ex_df[['chr', 'start', 'end', 'gene']]

        # ensure exons are indexed continuously.
        ex_df.reset_index(inplace=True
                          , drop=True)

        # Step 1: initialize GenomicArrayOfSets object. See https://bit.ly/2SPko6J for docs.
        gas = HTSeq.GenomicArrayOfSets('auto'
                                       , stranded=False)

        # Step 2: Add exon features to gas.
        for i in range(ex_df.shape[0]):
            chrom, exon_start, exon_end, gene = ex_df.iloc[i].values
            exon_intrv = HTSeq.GenomicInterval(chrom, exon_start, exon_end + 1, ".")
            gas[exon_intrv] += str(i)

        # Step 3. Find exons in exon_df that overlap with other exons.
        overlap_exons = list()
        for chrom in ex_df.chr.unique():
            # determine effective end of chromosome.
            chrom_end = ex_df[ex_df.chr == chrom].end.max() + 1

            # extract indices of chromosome's overlapping exons.
            chrom_steps = [(st[0], sorted(st[1])) for st in
                           gas[HTSeq.GenomicInterval(chrom, 0, chrom_end + 1, '.')].steps()]
            chrom_ol_exons = [chrom_steps[i][1] for i in range(len(chrom_steps)) if len(chrom_steps[i][1]) > 1]

            # cast exon id strings to int, extend set of overlap_exons.
            chrom_ol_exons = [[int(x) for x in y] for y in chrom_ol_exons]
            overlap_exons += chrom_ol_exons

        # Step 4: among exons that overlap with other exons, identify ones that overlap on different genes.
        # E.g., two exons may overlap but on the same gene, or, may overlap on different genes.
        rm_exons = list()
        for exons in overlap_exons:
            ex_sub_df = ex_df.iloc[exons]

            # identify exons forming part of an overlap region. If > 1 gene, append to exons for removal.
            if len(ex_sub_df.gene.unique()) > 1:
                rm_exons += exons

        # if there are any exons to remove, remove them.
        if len(rm_exons) > 0:
            ex_df.drop(rm_exons
                       , inplace=True)
        
        return ex_df

    def run(self):
        """
        Main function for GeneAnnotationProcessor. Runs transcriptome annotation processing pipeline:

         1. loads .gtf or .gff file for exon regions
         2. removes genes that occur in multiple chromosomes
         3. eliminates exons that intersect with > 1 genes on a single chromosome

        :param chroms: list of str chromosomes with which to subset data. Use if
        only a subset of chromosomes from genome annotation file will be useful.
        :return: pandas.DataFrame with 'chr', 'gene', 'gene_start', 'gene_end',
        'start' [exon start], and 'end'[exon end] fields
        outlining exons ready to use for reads coverage computations
        """
        if self.verbose:
            logging.info('Reading genome annotation file {0}.'.format(self.filename))

        exon_df = self.load()

        logging.info('Removing genes that occur in multiple chromosomes.')
        exon_df = self.remove_multichrom_genes(exon_df)
        exon_df.drop_duplicates(inplace=True)

        if self.verbose:
            logging.info('Successfully removed multiple-chromosome genes.')
            logging.info('Removing exons that occur in multiple genes.')

        gene_df = self.gene_outline(exon_df)
        n_exons = exon_df.shape[0]
        exon_df = self.remove_multigene_exons(exon_df)

        if (exon_df.shape[0] < n_exons) and (self.verbose):
            logging.info('{0} exons with multi-gene overlap were removed.'.format(n_exons - exon_df.shape[0]))

        exon_df = exon_df.merge(gene_df
                                , on=['chr', 'gene'])
        exon_df.drop_duplicates(inplace=True)

        if self.verbose:
            logging.info('Genome annotation file processing successful. Final shape -- {0}'.format(exon_df.shape))

        return exon_df


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