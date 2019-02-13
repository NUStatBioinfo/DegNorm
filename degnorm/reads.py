import pickle as pkl
import HTSeq
import networkx as nx
from pandas import DataFrame, concat, IntervalIndex
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
    def get_gene_overlap_structure(gene_df):
        """
        Build gene intersection matrix (a type of adjacency matrix), splitting genes into groups
        of mutually overlapping genes (i.e. groups of genes that are reachable on paths within
        adjacency matrix) and isolated genes that have no overlap with others.

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

        # assess how many genes we have.
        n_genes = chrom_gene_df.shape[0]

        # ---------------------------------------------------------------------- #
        # Step 1. Load chromosome's reads and index them.
        # ---------------------------------------------------------------------- #
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

        # easy wins: drop reads whose start position is < minimum start position of a gene,
        # and drop reads whose end position is > maximum start position of a gene
        min_gene_start, max_gene_end = chrom_gene_df.gene_start.min(), chrom_gene_df.gene_end.max()
        reads_df = reads_df[(reads_df.pos >= min_gene_start) & (reads_df.end_pos <= max_gene_end)]

        # If working with paired reads,
        # ensure that we've sequestered paired reads (eliminate any query names only occurring once).
        if self.paired:
            qname_counts = reads_df.qname_unpaired.value_counts()
            paired_occ_reads = qname_counts[qname_counts == 2].index.values.tolist()
            reads_df = reads_df[reads_df.qname_unpaired.isin(paired_occ_reads)]

        # ---------------------------------------------------------------------- #
        # Step 2. Drop reads that don't fully fall within union of all exons.
        # ---------------------------------------------------------------------- #
        chrom_len = self.header[self.header.chr == chrom].length.iloc[0]
        tscript_vec = np.ones([chrom_len]
                              , dtype=int)  # costly memory-wise!

        # build binary 0/1 exon/intron indicator vector.
        # Need to account for exon data being 1-indexed, tscript_vec is 0-indexed.
        exon_starts = chrom_exon_df.start.values - 1
        exon_ends = chrom_exon_df.end.values - 1
        for i in range(len(exon_starts)):
            tscript_vec[exon_starts[i]:exon_ends[i]] = 0

        del exon_starts, exon_ends
        gc.collect()

        # use values array, faster access.
        dat = reads_df[['cigar', 'pos', 'read_id']].values

        # store read_ids of reads to drop.
        drop_reads = list()

        # store read match region bounds, so that we only parse CIGAR strings once.
        read_bounds = list()

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
                    bounds_2 = [min_bounds_1 - 1 if j >= min_bounds_1 else j for j in bounds_2]
                    bounds_2.sort()

                # aggregate read pair's bounds.
                bounds = bounds_1 + bounds_2

                # iterate over match regions. If a single region is not fully contained
                # within exon regions, drop the pair.
                drop_read = False
                for j in np.arange(1, len(bounds), step=2):

                    # check whether matching regions on tscript_vec are fully contained within exonic regions.
                    # note that bounds are 1-indexed, tscript_vec is 0-indexed.
                    if np.sum(tscript_vec[(bounds[j - 1] - 1):bounds[j]]) > 0:
                        drop_read = True

                # append read id to set of read indices to drop (if appropriate).
                if drop_read:
                    drop_reads.extend([dat[ii - 1, 2], dat[ii, 2]])

                # otherwise, append match region bounds list. Note: endpoints of regions are inclusive.
                else:
                    read_bounds.append(bounds)

        # for single-read RNA-Seq experiments, we do not need such special consideration.
        else:
            for ii in np.arange(dat.shape[0]):
                # obtain read regions bounds.
                bounds = cigar_segment_bounds(dat[ii, 0]
                                              , start=dat[ii, 1])

                # iterate over match regions. If a single region is not fully contained
                # within exon regions, drop the read.
                drop_read = False
                for j in np.arange(1, len(bounds), step=2):

                    if np.sum(tscript_vec[(bounds[j - 1] - 1):bounds[j]]) > 0:
                        drop_read = True

                # append read id to set of read indices to drop (if appropriate).
                if drop_read:
                    drop_reads.append(dat[ii, 2])

                # otherwise, append match region bounds list. Note: endpoints of regions are inclusive.
                else:
                    read_bounds.append(bounds)

        # drop reads that don't fully intersect exonic regions.
        if drop_reads:
            reads_df.drop(drop_reads
                          , inplace=True)

        if self.paired:
            # if paired reads, don't actually need .1 and .2 constituent reads anymore.
            # So to save time + memory, take every other read.
            reads_df = reads_df.iloc[np.arange(1, reads_df.shape[0], step=2)]

        # add parsed match region bounds to reads!
        reads_df['bounds'] = read_bounds

        # drop the long chromosome-length vector, dropped read ids, reads_df values array.
        del tscript_vec, drop_reads, dat
        gc.collect()

        # ---------------------------------------------------------------------- #
        # Step 3. Compute coverage, reads across groups of mutually overlapping genes.
        # (This is costly from a time perspective. Should constitute
        #  coverage, read count calculations for ~ 10-20% of genes.)
        # ---------------------------------------------------------------------- #

        # initialize read count dictionary.
        read_count_dict = {gene: 0 for gene in chrom_gene_df.gene}

        # get gene intersection structure to guide ambiguous read detection.
        gene_overlap_dat = self.get_gene_overlap_structure(chrom_gene_df)
        n_isolated_genes = len(gene_overlap_dat['isolated_genes'])

        # display summary statistics around rate of gene intersection.
        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} -- {2} / {3} genes have overlap with others.\n'
                         'Begin computing coverage, read count for overlapping genes.'
                         .format(self.sample_id, chrom, n_genes - n_isolated_genes, n_genes))

        # initialize chromosome coverage array.
        cov_vec = np.zeros([chrom_len]
                           , dtype=int)

        # for genes in a group of overlapping genes, compute read coverage + count.
        if gene_overlap_dat['overlap_genes']:

            # iterate over groups of overlapping genes.
            for intersect_genes in gene_overlap_dat['overlap_genes']:

                intersect_gene_df = chrom_gene_df[chrom_gene_df.gene.isin(intersect_genes)]
                intersect_gene_start = intersect_gene_df.gene_start.min()
                intersect_gene_end = intersect_gene_df.gene_end.max()

                # obtain exon regions for each gene in intersection group.
                transcript_idx = list()
                for igene in intersect_genes:
                    igene_exon_df = chrom_exon_df[chrom_exon_df.gene == igene]
                    e_starts, e_ends = np.sort(igene_exon_df.start.values), np.sort(igene_exon_df.end.values)
                    exon_bounds = np.concatenate([[e_starts[j], e_ends[j]] for j in range(len(e_starts))])
                    transcript_idx.append(fill_in_bounds(exon_bounds
                                                         , endpoint=True))

                # storage for reads to drop.
                drop_reads = list()

                # subset reads to those that start and end within scope of this bloc of intersecting genes.
                intersect_reads_dat = reads_df[(reads_df.pos >= intersect_gene_start) &
                                               (reads_df.end_pos <= intersect_gene_end)][['bounds', 'read_id']].values()

                # for single-read RNA-Seq experiments, we do not need such special consideration.
                for ii in np.arange(intersect_reads_dat.shape[0]):

                    # obtain read regions bounds.
                    bounds, read_id = intersect_reads_dat[ii, :].values

                    # obtain read positions, shift by -1 so we zero-index reads.
                    read_idx = fill_in_bounds(bounds
                                              , endpoint=True) - 1

                    # find genes that fully include this read.
                    caught_genes = self.determine_full_inclusion(read_idx
                                                                 , gene_idx_list=transcript_idx)

                    # Ambiguous read determination logic:
                    # - if paired reads lie fully within 0 or 2+ genes, do not use the reads pair and drop them.
                    # - if read lies fully within a single gene:
                    #   - do not drop it.
                    #   - if the caught gene is the current gene being analyzed, use the read. Otherwise, do not.
                    drop_read = False
                    read_gene = None

                    # if more than one gene fully captures read, drop and do not use.
                    if len(caught_genes) != 1:
                        use_read = False
                        drop_read = True

                    # if only one gene captures read, use the read + identify capturing gene for incrementing count.
                    else:
                        use_read = caught_genes[0] == 0
                        read_gene = intersect_genes[caught_genes[0]]

                    # if only full intersection is with with a single gene, increment coverage and read count
                    # for that gene.
                    if use_read:
                        cov_vec[read_idx] += 1
                        read_count_dict[read_gene] += 1

                    # if need be, add read to list of reads to be dropped.
                    if drop_read:
                        drop_reads.append(read_id)

                # drop ambiguous reads from larger set of chromosome reads,
                # should speed up gene-read searches in the future.
                if drop_reads:
                    reads_df.drop(drop_reads
                                  , inplace=True)

            # free up some memory -- delete groups of intersecting genes.
            del gene_overlap_dat['overlap_genes'], intersect_reads_dat
            gc.collect()

        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} -- overlapping gene coverage, read counting successful.'
                         .format(self.sample_id, chrom))

        # ---------------------------------------------------------------------- #
        # Step 4. Compute coverage, reads for invidual isolated genes.
        # ---------------------------------------------------------------------- #
        if n_isolated_genes > 0:

            if self.verbose:
                logging.info('SAMPLE {0}: CHROMOSOME {1} -- Begin computing coverage, read count for isolated genes.'
                             .format(self.sample_id, chrom))

            # reduce chrom_gene_df to remaining genes
            chrom_gene_df = chrom_gene_df[chrom_gene_df.gene.isin(gene_overlap_dat['isolated_genes'])]

            # run same inclusion/exclusion transcript test but on the isolated genes.
            tscript_vec = np.ones([chrom_len]
                                  , dtype=int)

            # identify regions of chromosome covered by isolated genes.
            gene_starts = chrom_gene_df.gene_start.values - 1
            gene_ends = chrom_gene_df.gene_end.values - 1
            for i in range(len(gene_starts)):
                tscript_vec[gene_starts[i]:gene_ends[i]] = 0

            # identify reads that do not fall within an isolated gene's (start, end).
            drop_reads = list()
            dat = reads_df[['pos', 'end_pos', 'read_id']].values
            for i in range(dat.shape[0]):
                read_start, read_end, read_id = dat[i, :]

                if np.sum(tscript_vec[(read_start - 1):(read_end - 1)]) > 0:
                    drop_reads.append(read_id)

            # drop large memory structures.
            del dat, gene_starts, gene_ends, tscript_vec

            if drop_reads:
                print('Dropping {0} reads that do not lie completely within isolated genes'.format(len(drop_reads)))
                reads_df.drop(drop_reads
                              , inplace=True)

            del drop_reads
            gc.collect()

            # only continue if we have any reads intersecting isolated genes (a precaution).
            if not reads_df.empty:

                # add IntervalIndex index to chromosome gene data
                chrom_gene_df.index = IntervalIndex.from_arrays(chrom_gene_df.gene_start
                                                                , right=chrom_gene_df.gene_end
                                                                , closed='both')

                # "join" genes on reads data! (so that reads are tied to genes, for read counting.)
                reads_df['gene'] = chrom_gene_df.loc[reads_df.pos].gene.values

                # loop over reads for isolated genes, incrementing read count and coverage.
                dat = reads_df[['bounds', 'gene']].values
                for i in range(dat.shape[0]):
                    bounds, gene = dat[i, :]

                    # obtain read positions, shift by -1 so we zero-index reads.
                    read_idx = fill_in_bounds(bounds
                                              , endpoint=True) - 1

                    # increment coverage and read count.
                    cov_vec[read_idx] += 1
                    read_count_dict[gene] += 1

                del dat

            # drop remaining large data structures.
            del chrom_gene_df, chrom_exon_df, reads_df
            gc.collect()

            if self.verbose:
                logging.info('SAMPLE {0}: CHROMOSOME {1} -- isolated gene coverage, read counting successful.'
                             .format(self.sample_id, chrom))

        # ---------------------------------------------------------------------- #
        # Step 5. Save output.
        # chromosome overage vector ->> compressed csr numpy array
        # chromosome read counts ->> .csv file
        # ---------------------------------------------------------------------- #
        cov_file = os.path.join(self.save_dir, 'coverage_' + self.sample_id + '_' + chrom + '.npz')
        count_file = os.path.join(self.save_dir, 'read_counts_' + self.sample_id + '_' + chrom + '.csv')

        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} saving csr-compressed coverage array to {2}'
                         .format(self.sample_id, chrom, cov_file))

        # save coverage vector as a compressed-sparse row matrix.
        sparse.save_npz(cov_file
                        , matrix=sparse.csr_matrix(cov_vec))

        del cov_vec

        # construct read count DataFrame from read count dictionary.
        read_count_df = DataFrame({'gene': list(read_count_dict.keys())
                                   , self.sample_id: list(read_count_dict.values())})

        del read_count_dict

        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} mean per-gene read count: {2:.4}'
                         .format(self.sample_id, chrom, read_count_df[self.sample_id].mean()))

        # write sample chromosome read counts to .csv for joining later.
        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} saving read counts {2}'
                         .format(self.sample_id, chrom, count_file))
        read_count_df.to_csv(count_file
                             , index=False)

        gc.collect()

        # return per-chromosome coverage array filename and per-chromosome read counts.
        return cov_file, count_file

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


# if __name__ == '__main__':
#     from degnorm.gene_processing import GeneAnnotationProcessor
#
#     data_path = '/Users/fineiskid/nu/jiping_research/DegNorm/degnorm/tests/data/'
#     bam_file = os.path.join(data_path, 'hg_small_1.bam')
#     bai_file = os.path.join(data_path, 'hg_small_1.bai')
#     gtf_file = os.path.join(data_path, 'chr1_small.gtf')
#
#     gtf_processor = GeneAnnotationProcessor(gtf_file)
#     exon_df = gtf_processor.run()
#     gene_df = exon_df[['chr', 'gene', 'gene_start', 'gene_end']].drop_duplicates().reset_index(drop=True)
#
#     output_dir = '/Users/fineiskid/nu/jiping_research/degnorm_test_files'
#
#     reader = BamReadsProcessor(bam_file=bam_file
#                                , index_file=bai_file
#                                , chroms=['chr1']
#                                , n_jobs=1
#                                , output_dir=output_dir
#                                , verbose=True)
#
#     sample_id = reader.sample_id
#     print('found sample id {0}'.format(sample_id))
#
#     cov_files, read_count_files = reader.coverage_read_counts(gene_df
#                                                               , exon_df=exon_df)
