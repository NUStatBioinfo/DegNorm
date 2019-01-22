import pickle as pkl
import HTSeq
from pandas import DataFrame, concat
from collections import OrderedDict
from degnorm.utils import *
from degnorm.loaders import BamLoader
from joblib import Parallel, delayed
from scipy import sparse


def cigar_segment_bounds(cigar, start):
    """
    Determine the start and end positions on a chromosome of a non-no-matching part of an
    RNA-seq read based on a read's cigar string.

    cigar string meaning: http://bioinformatics.cvr.ac.uk/blog/tag/cigar-string/

    Example:
        '50M25N50M' with start = 100 -> [100, 149, 175, 224]. Note that start and end integers
        are inclusive, i.e. all positions at or between 100 and 149 and at or between 175 and 224
        are covered by reads.

    :param cigar: str a read's cigar string, e.g. "49M165N51M"
    :param start: int a read's start position on a chromosome
    :return: list of integers representing cigar match start, end points, in order of matching subsequences
    """
    # if CIGAR string is a single full match (i.e. "<positive integer>M")
    # extract length of the match, return match segment.
    full_match = re.match('(\d+)M$', cigar)
    if full_match is not None:
        extension = int(cigar[:(full_match.span()[-1] - 1)]) - 1

        return [start, start + extension]

    # break up cigar string into list of 2-tuples (letter indicative of match/no match, run length integer).
    cigar_split = [(v, int(k)) for k, v in re.findall(r'(\d+)([A-Z]?)', cigar)]

    # output storage.
    match_idx_list = list()

    for idx in range(len(cigar_split)):
        segment = cigar_split[idx]

        if segment[0] == 'M':
            extension = segment[1] - 1  # end of a match run is inclusive.
            augment = True
            match_idx_list += [start, start + extension]  # append a match run to output.

        else:
            if augment:
                extension = segment[1] + 1
                augment = False
            else:
                extension = segment[1]

        start += extension

    return match_idx_list


def fill_in_bounds(bounds_vec, endpoint=False):
    """
    Fill in the outline of contiguous integer regions with integers. For example:
    [10, 13, 20, 24] -> [10, 11, 12, 20, 21, 22, 23]

    :param bounds_vec: list of int or 1-d numpy array of int, outline of contiguous integer regions. Must
    have even number of elements.
    :param endpoint: bool should odd-indexed elements of bounds_vec (the region endpoints) be included
    in the filled in output? Same as numpy.linspace `endpoint` parameter.
    :return: 1-d numpy array of int
    """
    n = len(bounds_vec)

    if n % 2 != 0:
        raise ValueError('bounds_vec must have even number of values.')

    if not endpoint:
        filled_in = np.concatenate([np.arange(bounds_vec[j - 1], bounds_vec[j])
                                    for j in np.arange(1, n, step=2)])
    else:
        filled_in = np.concatenate([np.arange(bounds_vec[j - 1], bounds_vec[j]+1)
                                    for j in np.arange(1, n, step=2)])

    return filled_in


class BamReadsProcessor():
    def __init__(self, bam_file, index_file, chroms=None, n_jobs=max_cpu(),
                 output_dir=None, unique_alignment=False, verbose=True):
        """
        Transcript coverage and read counts processor, for a single alignment file (.bam).
        The main method for this class is coverage_read_counts, which computes coverage arrays and read counts
        for each gene, on in parallel over chromosomes.

        :param bam_file: str .bam filename
        :param index_file: str corresponding .bai (.bam index file) filename
        :param output_dir: str path to DegNorm output directory where coverage array files will be saved.
        If not specified, will use directory where RNA Seq experiment file is located.
        :param chroms: list of str names of chromosomes to load.
        :param n_jobs: int number of CPUs to use for determining genome coverage. Default
        is number of CPUs on machine - 1.
        :param unique_alignment: bool indicator - drop reads with NH:i:<x> flag where x > 1.
        :param verbose: bool indicator should progress be written to logger?
        """
        self.filename = bam_file

        # determine directory where we can dump coverage, reads files.
        if not output_dir:
            output_dir = os.path.join(os.path.dirname(self.filename), 'tmp')

        file_basename = '.'.join(os.path.basename(self.filename).split('.')[:-1])
        self.index_filename = index_file
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.sample_id = file_basename
        self.save_dir = os.path.join(output_dir, self.sample_id)
        self.header = None
        self.paired = None
        self.chroms = chroms
        self.unique_alignment = unique_alignment
        self.loader = BamLoader(self.filename, self.index_filename)
        self.get_header()
        self.determine_if_paired()

    def get_header(self):
        """
        Parse the header of a .bam file and extract the chromosomes and corresponding lengths
        of the chromosomes with reads contained in the .bam file. E.g.

        +----------------+-----------------+
        |     chr        |      length     |
        +================+=================+
        |     chr1       |     23445432    |
        +----------------+-----------------+
        |     chr2       |     31192127    |
        +----------------+-----------------+
        """

        # open .bam file connection.
        bam_file = self.loader.get_data()
        chrom_len_dict = dict()

        # parse header contained within .bam file.
        header_dict = bam_file.header.as_dict()['SQ']

        # close .bam file connection.
        bam_file.close()

        for header_line in header_dict:
            chrom_len_dict[header_line.get('SN')] = header_line.get('LN')

        # cast header as a pandas.DataFrame.
        self.header = DataFrame(list(chrom_len_dict.items())
                                , columns=['chr', 'length'])

        # based on supplied chromosome set and chromosomes in header, take intersection.
        if self.chroms is not None:
            self.chroms = np.intersect1d(self.chroms, self.header.chr.unique()).tolist()

        # if no chromosomes specified, assume load includes every chromosome in header.
        else:
            self.chroms = self.header.chr.unique().tolist()

    def determine_if_paired(self):
        """
        Determine if a .bam file is from a paired read experiment or from a single-end read experiment
        by studying the pattern of the first 500 reads' query names. Checks for "<query_name>.1" "<query_name>.2"
        pattern, indicating paired reads.
        """
        self.paired = False
        bam_file = self.loader.get_data()

        # pull first 300 queries' query names.
        ctr = 0
        qnames = list()
        for read in bam_file.fetch(self.chroms[0]):
            qnames.append(read.query_name)
            ctr += 1

            if ctr > 300:
                break

        # close .bam file connection.
        bam_file.close()

        # check if first 300 queries match the pattern of a query string from a paired read experiment.
        pair_indices = set(list(map(lambda x: x.split('.')[-1], qnames)))
        if pair_indices == {'1', '2'}:
            self.paired = True

    def load_chromosome_reads(self, chrom):
        """
        Load the reads from a .bam file for one particular chromosome.

        +-------------------+-------------+--------------+-------------------------------------+
        |        qname      |      pos    |     cigar    |  qname_unpaired [for paired reads]  |
        +===================+=============+==============+=====================================+
        | SRR873838.292.1   |   46189662  |      101M    |              SRR873838.292          |
        +-------------------+-------------+--------------+-------------------------------------+
        | SRR873838.292.2   |   46189763  |  77M255N24M  |              SRR873838.292          |
        +-------------------+-------------+----------------------------------------------------+

        :param chroms: list of str names of chromosomes to load.
        :return: pandas.DataFrame. If for paired reads file, additionally comes with qname_unpaired column and
        is sorted by qname_unpaired
        """
        reads = list()
        read_attributes = ['query_name', 'pos', 'cigarstring']
        bam_file = self.loader.get_data()

        for read in bam_file.fetch(chrom):

            # if working only with unique alignment reads, skip read if NH tag is > 1.
            if self.unique_alignment:
                if read.has_tag('NH'):
                    if read.get_tag('NH') > 1:
                        continue

            # if reading paired reads and the read is paired,
            # then grab the attributes of interest from the read in question.
            if self.paired:
                # pysam encodes RNEXT field as integer: -1 for "*" and 15 for "="
                if read.rnext != -1:
                    reads.append([getattr(read, attr) for attr in read_attributes])

            # otherwise (single-end reads) just grab them.
            else:
                reads.append([getattr(read, attr) for attr in read_attributes])

        # close .bam file connection.
        bam_file.close()

        # transform .bam lines into pandas.DataFrame
        df = DataFrame(reads
                       , columns=['qname', 'pos', 'cigar'])
        df['pos'] = df['pos'].astype('int')

        # remove reads list (large memory).
        del reads
        gc.collect()

        # if working with paired data: sort reads by unpaired query name.
        if self.paired:
            df['qname_unpaired'] = df.qname.apply(lambda x: '.'.join(x.split('.')[:-1]))
            df.sort_values('qname_unpaired', inplace=True)

        return df

    @staticmethod
    def get_gene_intersections(gene_df, exon_df):
        """
        Remove reads that meet the following criteria for being useless:
         1. Reads without any overlap with a single gene
         2. Reads that overlap with multiple genes, as such reads come from an ambiguous gene transcript.

        :param gene_df: pandas.DataFrame containing a chromosome's genes' location data, has at least these columns:
        `gene` (str name of gene), `gene_start` (int leftmost base position of gene), `gene_end` (int rightmost
        base position of gene)
        :param exon_df: pandas.DataFrame containing a chromosome's genes' exon composition, has at least these columns:
        `gene` (str name of gene), `start` (int leftmost exon base position), `end` (int rightmost exon base position),
        each row is an exon, one gene can map to multiple exons.
        :return: dict of OrderedDicts. Each key is a gene, each subkey is another gene that intersects with the first
        key. Values of subkeys are 1-d numpy arrays of exon boundaries [E1 start, E1 end, E2 start, E2 end, ...]
        """
        # initialize output storage.
        gene_intersection_dat = dict()

        # initialize GenomicArrayOfSets to hold gene locations.
        dat = exon_df[['gene', 'start', 'end']].values
        gas = HTSeq.GenomicArrayOfSets(['chrom'], stranded=False)

        # build gene locator gas.
        for i in range(dat.shape[0]):
            iv = HTSeq.GenomicInterval('chrom', dat[i, 1], dat[i, 2], '.')
            gas[iv] += dat[i, 0]

        # iterate over genes, find other genes that have overlap.
        # for genes with overlap, store their exon region bounds.
        for i in range(gene_df.shape[0]):
            gene, gene_start, gene_end = gene_df[['gene', 'gene_start', 'gene_end']].iloc[i].values
            gene_intersection_dat[gene] = OrderedDict()

            # search gas for overlapping genes.
            gas_intersect = [(st[0], sorted(st[1])) for
                             st in gas[HTSeq.GenomicInterval('chrom', gene_start, gene_end, '.')].steps()]

            # parse results from gas search.
            gene_intersect = list()
            for ii in range(len(gas_intersect)):
                cross_genes = gas_intersect[ii][1]
                if cross_genes:
                    gene_intersect.extend(cross_genes)

            # reorder union of intersecting genes so that gene of interest leads.
            gene_intersect = list(set(gene_intersect))
            gene_intersect.remove(gene)
            gene_intersect = [gene] + gene_intersect

            # for overlapping genes, identify corresponding exon boundaries.
            for int_gene in gene_intersect:
                my_exon_df = exon_df[exon_df.gene == int_gene]
                e_starts, e_ends = np.sort(my_exon_df.start.values), np.sort(my_exon_df.end.values)
                gene_intersection_dat[gene][int_gene] = np.concatenate(
                    [[e_starts[j], e_ends[j]] for j in range(len(e_starts))])

        return gene_intersection_dat

    @staticmethod
    def determine_full_inclusion(read_idx, gene_idx_list):
        """
        Determine which genes' exon regions fully include a read's base positions.
        For example:
        read_idx: [3, 4, 5]
        gene_idx_list: [[1, 2, 3, 4], [2, 3, 4, 5, 6], [11, 14, 15, 16]]
        -> array([1])

        :param read_idx: list of int or 1-d numpy array of int, read's base positions
        :param gene_idx_list: list of list of int or 1-d numpy array of int, a set of genes' exon positions,
        one sublist per gene.
        :return: 1-d numpy array with integer indices of genes in gene_idx_list that fully include read's bases
        """
        gene_capture = np.where(list(map(lambda y: len(np.setdiff1d(read_idx, y)) == 0, gene_idx_list)))[0]

        return gene_capture

    def chromosome_coverage_read_counts(self, chrom_gene_df, chrom_exon_df, chrom):
        """
        Determine per-chromosome reads coverage and per-gene read counts from an RNA-seq experiment in
        a way that properly considers ambiguous reads - if a (paired) read falls entirely within the
        exonic regions of a *single* gene, only then does read contribute to read count and coverage.
        The cigar scores from single and paired reads are parsed according to cigar_segment_bounds.

        1. Saves compressed coverage array to self.save_dir with file name 'sample_[sample_id]_[chrom].npz'
        2. Saves read counts to self.save_dir with file name 'read_counts_[sample_id]_[chrom].csv'

        It is assumed that overlapping gene exons have been removed from exon set.
        See gene_processing.GeneAnnotationProcessor.remove_multigene_exons method.

        :param chrom_gene_df: pandas.DataFrame with `chr`, `gene`, `gene_start`, and `gene_end` columns
        that delineate the start and end position of a gene's transcript on a chromosome, must be
        subset to the chromosome in study.
        :param chrom_exon_df: pandas.DataFrame with `chr`, `gene`, `start`, `end` columns that delineate
        the start and end positions of exons on a gene.
        :param chrom: str chromosome name
        :return: str full file path to where coverage array is saved in a compressed .npz file.
        """
        # First, load this chromosome's reads.
        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} begin loading reads from {2}'
                         .format(self.sample_id, chrom, self.filename))

        reads_df = self.load_chromosome_reads(chrom)

        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} reads successfully loaded -- shape: {2}'
                         .format(self.sample_id, chrom, reads_df.shape))

        # append end position to reads based on cigar score.
        reads_df['end_pos'] = reads_df['pos'] + reads_df['cigar'].apply(
            lambda x: sum([int(k) for k, v in re.findall(r'(\d+)([A-Z]?)', x)]))

        # assign row number to read ID column, then use as index.
        reads_df['read_id'] = range(reads_df.shape[0])
        reads_df.set_index('read_id'
                           , drop=False
                           , inplace=True)

        # If working with paired reads,
        # ensure that we've sequestered paired reads (eliminate any query names only occurring once).
        if self.paired:
            qname_counts = reads_df.qname_unpaired.value_counts()
            paired_occ_reads = qname_counts[qname_counts == 2].index.values.tolist()
            reads_df = reads_df[reads_df.qname_unpaired.isin(paired_occ_reads)]

        # get gene intersection data to guide ambiguous read detection.
        gene_intersection_dat = self.get_gene_intersections(chrom_gene_df
                                                            , exon_df=chrom_exon_df)

        # start getting coverage: iterate over genes.
        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} begin computing coverage for {2} genes.'
                         .format(self.sample_id, chrom, chrom_gene_df.shape[0]))

        read_counts = np.zeros(chrom_gene_df.shape[0])
        gene_cov_dat = dict()
        for i in range(chrom_gene_df.shape[0]):
            # extract this gene's start, end positions.
            gene, gene_start, gene_end = chrom_gene_df[['gene', 'gene_start', 'gene_end']].iloc[i].values

            # initialize gene coverage array (will be pared down to exon regions later).
            cov_vec = np.zeros(gene_end - gene_start + 1)

            # create GenomicArrayOfSets for gene-read search.
            gas = HTSeq.GenomicArrayOfSets(['chrom'], stranded=False)
            for intersect_gene in gene_intersection_dat[gene]:
                exon_bounds = gene_intersection_dat[gene][intersect_gene]

                # append to gas the exonic regions comprising gene.
                for ii in np.arange(1, stop=len(exon_bounds), step=2):
                    iv = HTSeq.GenomicInterval('chrom', exon_bounds[ii - 1], exon_bounds[ii], '.')
                    gas[iv] += intersect_gene

            # obtain full set of gene exon positions for coverage dicing later on.
            transcript_idx = fill_in_bounds(gene_intersection_dat[gene][gene])

            # storage for reads to drop.
            drop_reads = list()
            read_count = 0

            # subset reads to those that start and end within scope of gene.
            dat = reads_df[((reads_df.pos >= gene_start) & (reads_df.end_pos <= gene_end))][['cigar', 'pos', 'read_id']].values

            # for paired reads, perform special parsing of CIGAR strings to avoid double-counting of overlap regions.
            if self.paired:
                for ii in np.arange(1, dat.shape[0], 2):

                    # obtain read region bounds.
                    bounds_1 = cigar_segment_bounds(dat[ii - 1, 0]
                                                    , start=dat[ii - 1, 1])
                    bounds_2 = cigar_segment_bounds(dat[ii, 0]
                                                    , start=dat[ii, 1])

                    # leverage nature of alignments of paired reads to find disjoint coverage ranges.
                    min_bounds_1, max_bounds_1 = min(bounds_1), max(bounds_1)
                    min_bounds_2, max_bounds_2 = min(bounds_2), max(bounds_2)

                    if max_bounds_2 >= max_bounds_1:
                        bounds_2 = [max_bounds_1 + 1 if j <= max_bounds_1 else j for j in bounds_2]
                    else:
                        # bounds_2 = list(set([min_bounds_1 - 1 if j >= min_bounds_1 else j for j in bounds_2]))
                        bounds_2 = [min_bounds_1 - 1 if j >= min_bounds_1 else j for j in bounds_2]
                        bounds_2.sort()

                    # aggregate read pair's bounds.
                    bounds = bounds_1 + bounds_2

                    # search for read's matching regions within gene exon gas:
                    # For each read matching region, find names of genes whose exons fully capture the region.
                    read_gene = None  # storage for the single capturing gene.
                    for j in np.arange(1, stop=len(bounds), step=2):

                        gas_steps = [(st[0], sorted(st[1])) for st in
                                     gas[HTSeq.GenomicInterval('chrom', bounds[j - 1], bounds[j] + 1, '.')].steps()]

                        # if match region broken up over multiple intervals, there is no hope for a full capture.
                        # Otherwise, need to check how many distinct genes' fully capture this match region with
                        # their exons. Note: if one gene has contiguous exons, e.g. [100, 130][130, 150], GAS
                        # considers this one GenomicInterval.
                        if len(gas_steps) == 1:
                            capture_genes = gas_steps[0][1]
                            n_capture_genes = len(capture_genes)

                            # only if number of fully capturing genes is exactly 1 for this match region,
                            # we might consider it.
                            if n_capture_genes == 1:
                                # if analyzing the first match region, the read_gene is set to the capturing gene.
                                if j == 1:
                                    read_gene = capture_genes[0]
                                # if analyzing another match region, only consider read_gene if it's equal to
                                # the first capturing gene.
                                else:
                                    if read_gene != capture_genes[0]:
                                        read_gene = None

                    # Ambiguous read determination logic:
                    use_read = False
                    drop_read = False

                    # if paired reads lie fully within only 1 gene, keep it.
                    if read_gene is not None:
                        # if that capturing gene is the gene in question, use the read.
                        # Otherwise, do not use read but do not drop read.
                        use_read = read_gene == gene

                    # if 0 or 2+ genes capture read, drop read and do not use.
                    else:
                        drop_read = True

                    if use_read:
                        cov_vec[fill_in_bounds(bounds, endpoint=True) - gene_start] += 1
                        read_count += 1

                    # append read pair to list of reads to drop if need be.
                    if drop_read:
                        drop_reads.extend([dat[ii - 1, 2], dat[ii, 2]])

            # for single-read RNA-Seq experiments, we do not need such special consideration.
            else:
                for ii in np.arange(dat.shape[0]):

                    # obtain read regions bounds.
                    bounds = cigar_segment_bounds(dat[ii, 0]
                                                  , start=dat[ii, 1])

                    # search for read's matching regions within gene exon gas:
                    # For each read matching region, find names of genes whose exons fully capture the region.
                    read_gene = None  # storage for the single capturing gene.
                    for j in np.arange(1, stop=len(bounds), step=2):

                        gas_steps = [(st[0], sorted(st[1])) for st in
                                     gas[HTSeq.GenomicInterval('chrom', bounds[j - 1], bounds[j] + 1, '.')].steps()]

                        # if match region broken up over multiple intervals, there is no hope for a full capture.
                        # Otherwise, need to check how many distinct genes' fully capture this match region with
                        # their exons. Note: if one gene has contiguous exons, e.g. [100, 130][130, 150], GAS
                        # considers this one GenomicInterval.
                        if len(gas_steps) == 1:
                            capture_genes = gas_steps[0][1]
                            n_capture_genes = len(capture_genes)

                            # only if number of fully capturing genes is exactly 1 for this match region,
                            # we might consider it.
                            if n_capture_genes == 1:
                                # if analyzing the first match region, the read_gene is set to the capturing gene.
                                if j == 1:
                                    read_gene = capture_genes[0]
                                # if analyzing another match region, only consider read_gene if it's equal to
                                # the first capturing gene.
                                else:
                                    if read_gene != capture_genes[0]:
                                        read_gene = None

                    # Ambiguous read determination logic:
                    use_read = False
                    drop_read = False

                    # if read lies fully within only 1 gene, keep it.
                    if read_gene is not None:
                        # if that capturing gene is the gene in question, use the read.
                        # Otherwise, do not use read but do not drop read.
                        use_read = read_gene == gene

                    # if 0 or 2+ genes capture read, drop read and do not use.
                    else:
                        drop_read = True

                    if use_read:
                        cov_vec[fill_in_bounds(bounds, endpoint=True) - gene_start] += 1
                        read_count += 1

                    # add read to list of reads to drop if need be.
                    if drop_read:
                        drop_reads.extend(dat[ii, 2])

            # store exonic regions of gene coverage array as sample coverage matrix.
            gene_cov_dat[gene] = cov_vec[transcript_idx - gene_start]

            # store read count.
            read_counts[i] = read_count

            # do not consider the current gene later on in the intersecting genes' reads parsing.
            for other_gene in list(gene_intersection_dat[gene].keys())[1:]:
                if gene in gene_intersection_dat[other_gene]:
                    del gene_intersection_dat[other_gene][gene]

            # drop ambiguous reads from larger set of chromosome reads,
            # should speed up gene-read searches in the future.
            if drop_reads:
                reads_df.drop(drop_reads
                              , inplace=True)

        # free up some memory -- delete large memory objects.
        del reads_df, dat, cov_vec, gene_intersection_dat, drop_reads
        gc.collect()

        # save gene coverage dictionary to disk as pickle file.
        pkl_file = os.path.join(self.save_dir, 'coverage_' + self.sample_id + '_' + chrom + '.pkl')

        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} saving {2} coverage arrays to {3}'
                         .format(self.sample_id, chrom, len(gene_cov_dat), pkl_file))

        with open(pkl_file, 'wb') as f:
            pkl.dump(gene_cov_dat, f)

        # free up memory, delete gene coverage data.
        del gene_cov_dat
        gc.collect()

        # build read counts DataFrame.
        read_count_df = DataFrame({'gene': chrom_gene_df.gene.values
                                      , self.sample_id: read_counts})

        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} mean per-gene read count: {2:.4}'
                         .format(self.sample_id, chrom, read_count_df[self.sample_id].mean()))

        # write sample chromosome read counts to .csv for joining later.
        read_count_file = os.path.join(self.save_dir, 'read_counts_' + self.sample_id + '_' + chrom + '.csv')
        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} saving read counts {2}'
                         .format(self.sample_id, chrom, read_count_file))
        read_count_df.to_csv(read_count_file
                             , index=False)

        # return per-chromosome coverage array filename and per-chromosome read counts.
        return pkl_file, read_count_file

    def coverage_read_counts(self, gene_df, exon_df):
        """
        Main function for computing coverage arrays in parallel over chromosomes.

        :param gene_df: pandas.DataFrame with `chr`, `gene`, `gene_start`, and `gene_end` columns
        that delineate the start and end position of a gene's transcript on a chromosome. See
        GeneAnnotationProcessor.
        :return: list of str file paths of compressed .npz files containing coverage arrays.
        """
        # create directory in DegNorm output dir where sample coverage vecs are saved.
        if not os.path.exists(self.save_dir):
            os.makedirs(self.save_dir)

        if self.verbose:
            logging.info('SAMPLE {0}: determining read coverage and read counts for {1} chromosomes.\n'
                         'Output will be saved to directory {2}.'
                         .format(self.sample_id, len(self.chroms), self.save_dir))

        # distribute work with joblib.Parallel:
        par_output = Parallel(n_jobs=min(self.n_jobs, len(self.chroms))
                              , verbose=0
                              , backend='threading')(delayed(self.chromosome_coverage_read_counts)(
            chrom_gene_df=subset_to_chrom(gene_df, chrom=chrom),
            chrom_exon_df=subset_to_chrom(exon_df, chrom=chrom),
            chrom=chrom)
            for chrom in self.chroms)

        # parse output from parallel workers.
        cov_filepaths = [x[0] for x in par_output]
        read_count_filepaths = [x[1] for x in par_output]

        return cov_filepaths, read_count_filepaths


if __name__ == '__main__':
    from degnorm.gene_processing import GeneAnnotationProcessor

    data_path = '/Users/fineiskid/nu/jiping_research/DegNorm/degnorm/tests/data/'
    bam_file = os.path.join(data_path, 'hg_small_1.bam')
    bai_file = os.path.join(data_path, 'hg_small_1.bai')
    gtf_file = os.path.join(data_path, 'chr1_small.gtf')

    gtf_processor = GeneAnnotationProcessor(gtf_file)
    exon_df = gtf_processor.run()
    gene_df = exon_df[['chr', 'gene', 'gene_start', 'gene_end']].drop_duplicates().reset_index(drop=True)

    output_dir = '/Users/fineiskid/nu/jiping_research/degnorm_test_files'

    reader = BamReadsProcessor(bam_file=bam_file
                               , index_file=bai_file
                               , chroms=['chr1']
                               , n_jobs=1
                               , output_dir=output_dir
                               , verbose=True)

    sample_id = reader.sample_id
    print('found sample id {0}'.format(sample_id))

    cov_files, read_count_files = reader.coverage_read_counts(gene_df
                                                              , exon_df=exon_df)
