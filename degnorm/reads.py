from pandas import DataFrame, IntervalIndex, read_sql_query
from degnorm.utils import *
from degnorm.loaders import BamLoader
from joblib import Parallel, delayed
from scipy import sparse
from sqlite3 import connect
import pickle as pkl


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
    full_match = re.match(r'(\d+)M$', cigar)
    if full_match is not None:
        extension = int(cigar[:(full_match.span()[-1] - 1)]) - 1

        return [start, start + extension]

    # break up cigar string into list of 2-tuples (letter indicative of match/no match, run length integer).
    cigar_split = [(v, int(k)) for k, v in re.findall(r'(\d+)([A-Z]?)', cigar)]

    # initialize parse params.
    # Allow for "hard clipping" where aligned read can start with non-matching region (https://bit.ly/2K6TJ5Y)
    augment = False
    any_match = False

    # output storage.
    match_idx_list = list()

    for idx in range(len(cigar_split)):
        segment = cigar_split[idx]

        if segment[0] == 'M':
            any_match = True
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

    # if no matching regions found, throw error.
    if not any_match:
        raise ValueError('CIGAR string {0} has no matching region.'.format(cigar))

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
        raise ValueError('bounds_vec = {0}, must have even number of values!'.format(bounds_vec))

    if not endpoint:
        filled_in = np.concatenate([np.arange(bounds_vec[j - 1], bounds_vec[j])
                                    for j in np.arange(1, n, step=2)])
    else:
        filled_in = np.concatenate([np.arange(bounds_vec[j - 1], bounds_vec[j]+1)
                                    for j in np.arange(1, n, step=2)])

    return filled_in


class BamReadsProcessor:

    def __init__(self, bam_file, index_file, chroms=None, n_jobs=1,
                 output_dir=None, unique_alignment=True, verbose=True):
        """
        Transcript coverage and read counts processor, for a single alignment file (.bam).
        The main method for this class is coverage_read_counts, which computes coverage arrays and read counts
        for each gene, on in parallel over chromosomes.

        :param bam_file: str .bam filename
        :param index_file: str corresponding .bai (.bam index file) filename
        :param output_dir: str path to DegNorm output directory where coverage array files will be saved.
        If not specified, will use directory where RNA Seq experiment file is located.
        :param chroms: list of str names of chromosomes to load.
        :param n_jobs: int number of threads to use while determining genome coverage.
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

        # tell user whether or not sample has been detected as either paired or single-end reads.
        if self.verbose:
            logging.info('SAMPLE {0} -- sample contains {1} reads'
                         .format(self.sample_id, 'paired' if self.paired else 'single-end'))


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
    def determine_full_inclusion(read_bounds, gene_exon_bounds):
        """
        Determine which genes fully encapsulate a set of read matching regions.

        For example: a given read has 2 matching regions that are captured by only the first gene's exons.

        read_bounds = [10, 32, 45, 90]
        gene_exon_bounds = [[[8, 40], [44, 100]], [[2, 20], [60, 400]]]

        ->>

        [0] # gene 0

        :param read_bounds: list or 1-d array of even length alternating between positions of read
        matching region starts, matching region ends
        :param gene_exon_bounds: list of list of lists, each sublist is a list of [exon start, exon end] subsublists,
        one sublist per gene.
        :return: list with integer indices of gene_exon_bounds corresponding to genes that fully capture read bounds.
        """
        full_capture_idx = list()

        # iterate over genes.
        for gene_idx in range(len(gene_exon_bounds)):
            exon_bounds = gene_exon_bounds[gene_idx]
            full_capture = True

            # iterate over read matching regions.
            for j in np.arange(1, len(read_bounds), step=2):
                seg_capture = False
                start = read_bounds[j - 1]
                end = read_bounds[j]

                # check gene's exon regions to see if they capture read match region.
                # Stop checking exon regions the moment there's one full match.
                for exon_bound in exon_bounds:
                    if start >= exon_bound[0] and end <= exon_bound[1]:
                        seg_capture = True
                        break

                # if any segment not captured, stop checking match segments
                # and exclude gene from set of fully capturing genes.
                if not seg_capture:
                    full_capture = False
                    break

            if full_capture:
                full_capture_idx.append(gene_idx)

        return full_capture_idx

    def chromosome_coverage_read_counts(self, gene_overlap_dat, chrom_gene_df, chrom_exon_df, chrom):
        """
        Determine per-chromosome reads coverage and per-gene read counts from an RNA-seq experiment in
        a way that properly considers ambiguous reads - if a (paired) read falls entirely within the
        exonic regions of a *single* gene, only then does read contribute to read count and coverage.
        The cigar scores from single and paired reads are parsed according to cigar_segment_bounds.

        1. Saves compressed coverage array to self.save_dir with file name 'sample_[sample_id]_[chrom].npz' for
        genes with no overlap with any other gene.
        2. Saves gene coverage arrays (only their exonic regions) to serialized pickle files.
        3. Saves read counts to self.save_dir with file name 'read_counts_[sample_id]_[chrom].csv'

        :param chrom_gene_df: pandas.DataFrame with `chr`, `gene`, `gene_start`, and `gene_end` columns
        that delineate the start and end position of a gene's transcript on a chromosome, must be
        subset to the chromosome in study.
        :param gene_overlap_dat: dictionary with keys 'isolated_genes' and 'overlap_genes' detailing
        groups of genes that do not overlap with others and then groups of genes that share any overlap.
        See gene_processing.get_gene_overlap_structure function.
        :param chrom_exon_df: pandas.DataFrame with `chr`, `gene`, `start`, `end` columns that delineate
        the start and end positions of exons on a gene.
        :param chrom: str chromosome name
        :return: None. Coverage and read count files are written to self.save_dir.
        """
        # First, load this chromosome's reads.
        if self.verbose:
            logging.info('SAMPLE {0}, CHR {1} -- begin loading reads from {2}'
                         .format(self.sample_id, chrom, self.filename))

        # assess how many genes we have.
        n_genes = chrom_gene_df.shape[0]

        # gene_overlap_dat data check: ensure that number isolated genes + number overlapping genes
        # equals number of genes in genes DataFrame.
        n_isolated_genes, n_overlap_genes = 0, 0
        if gene_overlap_dat['isolated_genes']:
            n_isolated_genes = len(gene_overlap_dat['isolated_genes'])

        if gene_overlap_dat['overlap_genes']:
            n_overlap_genes = np.sum([len(x) for x in gene_overlap_dat['overlap_genes']])

        if n_isolated_genes + n_overlap_genes != n_genes:
            raise ValueError('number of genes contained in gene_overlap_dat does not match that of chrom_gene_df.')

        # initialize read counts.
        read_count_dict = {gene: 0 for gene in chrom_gene_df.gene}

        # ---------------------------------------------------------------------- #
        # Step 1. Load chromosome's reads and index them.
        # ---------------------------------------------------------------------- #
        reads_df = self.load_chromosome_reads(chrom)

        if self.verbose:
            logging.info('SAMPLE {0}, CHR {1} -- reads successfully loaded. shape = {2}'
                         .format(self.sample_id, chrom, reads_df.shape))

        # append end position to reads based on cigar score.
        reads_df['end_pos'] = reads_df['pos'] + reads_df['cigar'].apply(
            lambda x: sum([int(k) for k, v in re.findall(r'(\d+)([A-Z]?)', x)]))

        # assign row number to read ID column, then use as index.
        reads_df['read_id'] = range(reads_df.shape[0])
        reads_df.set_index('read_id'
                           , drop=False
                           , inplace=True)

        # easy win: drop reads whose start position is < minimum start position of a gene,
        # and drop reads whose end position is > maximum start position of a gene
        min_gene_start, max_gene_end = chrom_gene_df.gene_start.min() - 1, chrom_gene_df.gene_end.max() - 1
        reads_df = reads_df[(reads_df.pos >= (min_gene_start)) & (reads_df.end_pos <= (max_gene_end))]

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
                              , dtype=int)  # large vector, will delete after using.

        # build binary 0/1 exon/intron indicator vector.
        # Need to account for exon data being 1-indexed, tscript_vec is 0-indexed, but
        # exon end positions are inclusive.
        exon_starts = chrom_exon_df.start.values - 1
        exon_ends = chrom_exon_df.end.values
        for i in range(len(exon_starts)):
            tscript_vec[exon_starts[i]:exon_ends[i]] = 0

        del exon_starts, exon_ends
        gc.collect()

        # store read_ids of reads to drop, and initialize dropped read count.
        drop_reads = list()

        # store read match region bounds, so that we only parse CIGAR strings once.
        read_bounds = list()

        # use values array, faster access.
        dat = reads_df[['cigar', 'pos', 'read_id']].values

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
                    # note that right-bounds are inclusive.
                    if np.sum(tscript_vec[(bounds[j - 1]):(bounds[j] + 1)]) > 0:
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

                    if np.sum(tscript_vec[(bounds[j - 1]):(bounds[j] + 1)]) > 0:
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

        # delete objs, attempt to save on memory.
        del tscript_vec, drop_reads, dat, read_bounds
        gc.collect()

        # ---------------------------------------------------------------------- #
        # Step 3. Compute coverage, reads across groups of mutually overlapping genes.
        # (This is costly from a time perspective. Should constitute
        #  coverage, read count calculations for ~ 10-20% of genes.)
        # ---------------------------------------------------------------------- #

        # display summary statistics around rate of gene intersection.
        if self.verbose:
            logging.info('SAMPLE {0}, CHR {1} -- overlap genes = {2} / {3}.'
                         .format(self.sample_id, chrom, n_overlap_genes, n_genes))
            logging.info('SAMPLE {0}, CHR {1} -- begin overlap gene group reads processing.'
                         .format(self.sample_id, chrom))

        # for genes in a group of overlapping genes, compute read coverage + count.
        if n_overlap_genes > 0:

            ol_cov_dict = dict()

            # iterate over groups of overlapping genes.
            for ol_genes in gene_overlap_dat['overlap_genes']:

                ol_gene_df = chrom_gene_df[chrom_gene_df.gene.isin(ol_genes)]
                ol_gene_group_start = ol_gene_df.gene_start.min() - 1
                ol_gene_group_end = ol_gene_df.gene_end.max() - 1

                ol_gene_starts = list()
                gene_exon_bounds = list()
                transcript_idx = list()

                # obtain exon regions for each gene in overlap group.
                # Exon starts/ends are 1-indexed, change them to be 0-indexed.
                for ol_gene in ol_genes:
                    ol_gene_exon_df = chrom_exon_df[chrom_exon_df.gene == ol_gene]

                    # store gene starts for constructing per-gene coverage vectors.
                    # 0-index gene starts/ends.
                    ol_gene_start = ol_gene_exon_df.gene_start.iloc[0] - 1
                    ol_gene_end = ol_gene_exon_df.gene_end.iloc[0] - 1
                    ol_gene_starts.append(ol_gene_start)

                    # initialize gene coverage vector for each gene in overlap group.
                    ol_cov_dict[ol_gene] = np.zeros([ol_gene_end - ol_gene_start + 1]
                                                    , dtype=int)

                    # save gene exon positioning, for determining which reads captured by which genes.
                    # 0-index exon positions, and include gene end positioning.
                    e_starts, e_ends = np.sort(ol_gene_exon_df.start.values) - 1, np.sort(ol_gene_exon_df.end.values)
                    gene_exon_bounds += [[[e_starts[j], e_ends[j]] for j in range(len(e_starts))]]      # list of list of lists, includes exon end pos.
                    transcript_idx.append(np.unique(fill_in_bounds(flatten_2d(gene_exon_bounds[-1]))))  # transcript vector is 0-indexed, includes exon end pos.

                # drop things we don't need any more.
                del ol_gene_df, ol_gene_exon_df, e_starts, e_ends

                # storage for reads to drop.
                drop_reads = list()

                # subset reads to those that start and end within scope of this bloc of overlapping genes.
                ol_reads_dat = reads_df[(reads_df.pos >= (ol_gene_group_start)) &
                                        (reads_df.end_pos <= (ol_gene_group_end))][['bounds', 'read_id']].values

                # for single-read RNA-Seq experiments, we do not need such special consideration.
                for i in range(ol_reads_dat.shape[0]):

                    # obtain read regions bounds.
                    read_bounds, read_id = ol_reads_dat[i, :]

                    # find genes that fully include this read. Everything is 0-indexed.
                    caught_genes = self.determine_full_inclusion(read_bounds
                                                                 , gene_exon_bounds=gene_exon_bounds)

                    # Ambiguous read determination logic:
                    # - if paired reads lie fully within 0 or 2+ genes, do not use the reads pair and drop them.
                    # - if read lies fully within a single gene:
                    #    - do not drop it.
                    #    - if the caught gene is the current gene being analyzed, use the read. O/w do not.
                    n_caught_genes = len(caught_genes)

                    # if only one gene captures read, use the read and identify capturing gene for
                    # incrementing count, but drop it from consideration later (it's been accounted for).
                    # if only full intersection is with with a single gene, increment coverage and read count
                    # for that gene, and drop read.
                    # Note: need to restart coverage calculations relative to gene's start position.
                    if n_caught_genes == 1:
                        drop_read = True
                        read_gene = ol_genes[caught_genes[0]]
                        read_gene_start = ol_gene_starts[caught_genes[0]]
                        read_idx = fill_in_bounds(read_bounds
                                                  , endpoint=True) - read_gene_start - 1
                        ol_cov_dict[read_gene][read_idx] += 1
                        read_count_dict[read_gene] += 1

                    # if no gene fully captures the read, do not use read *but do not drop it*,
                    # for the possibility that some isolated gene captures the read later on.
                    elif n_caught_genes == 0:
                        drop_read = False

                    # if > 1 gene fully captures the read,
                    # do not use read and drop it from consideration.
                    else:
                        drop_read = True

                    # if need be, add read to list of reads to be dropped.
                    if drop_read:
                        drop_reads.append(read_id)

                # drop ambiguous reads from larger set of chromosome reads,
                # should speed up gene-read searches in the future.
                if drop_reads:
                    reads_df.drop(drop_reads
                                  , inplace=True)

                    del drop_reads

                # pare down coverage vectors for genes in overlap grou to their concatenated exon regions.
                for i in range(len(ol_genes)):
                    ol_gene = ol_genes[i]
                    ol_cov_dict[ol_gene] = ol_cov_dict[ol_gene][transcript_idx[i] - ol_gene_starts[i]]

            # ---------------------------------------------------------------------- #
            # Step 3.5: save overlapping genes' coverage vectors.
            # overlapping gene coverage vector dict ->> pkl file.
            # ---------------------------------------------------------------------- #
            ol_cov_file = os.path.join(self.save_dir, 'overlap_coverage_' + self.sample_id + '_' + chrom + '.pkl')
            if self.verbose:
                logging.info('SAMPLE {0}, CHR {1} -- saving overlapping gene coverage vectors.'
                             .format(self.sample_id, chrom))

            # dump overlapping genes' coverage matrices.
            with open(ol_cov_file, 'wb') as f:
                pkl.dump(ol_cov_dict, f)

            # free up some memory -- delete groups of intersecting genes, etc.
            del ol_reads_dat, ol_cov_dict, transcript_idx, gene_exon_bounds
            gc.collect()

            if self.verbose:
                logging.info('SAMPLE {0}, CHR {1} -- overlapping gene reads processing successful.'
                             .format(self.sample_id, chrom))

        # ---------------------------------------------------------------------- #
        # Step 4. Compute coverage, reads for individual isolated genes.
        # ---------------------------------------------------------------------- #
        if n_isolated_genes > 0:

            if self.verbose:
                logging.info('SAMPLE {0}, CHR {1} -- begin isolated gene reads processing.'
                             .format(self.sample_id, chrom))

            # reduce chrom_gene_df to remaining genes
            chrom_gene_df = chrom_gene_df[chrom_gene_df.gene.isin(gene_overlap_dat['isolated_genes'])]

            # run same inclusion/exclusion transcript test but on the isolated genes.
            tscript_vec = np.ones([chrom_len]
                                  , dtype=int)

            # identify regions of chromosome covered by isolated genes.
            # change gene starts/ends to 0-indexed to match 0-indexed tscript_vec array, but
            # gene ends are inclusive.
            gene_starts = chrom_gene_df.gene_start.values - 1
            gene_ends = chrom_gene_df.gene_end.values
            for i in range(len(gene_starts)):
                tscript_vec[gene_starts[i]:gene_ends[i]] = 0

            # identify reads that do not fall within an isolated gene's (start, end).
            drop_reads = list()
            dat = reads_df[['pos', 'end_pos', 'read_id']].values
            for i in range(dat.shape[0]):
                read_start, read_end, read_id = dat[i, :]

                # remember to include read end position. reads are 0-indexed.
                if np.sum(tscript_vec[read_start:(read_end + 1)]) > 0:
                    drop_reads.append(read_id)

            # drop memory hogs.
            del dat, gene_starts, gene_ends, tscript_vec

            # drop reads that do not lie completely within area covered by isolated genes.
            if drop_reads:
                reads_df.drop(drop_reads
                              , inplace=True)

            del drop_reads
            gc.collect()

            # (a precaution) only continue if we have any reads intersecting isolated genes.
            if not reads_df.empty:

                # initialize chromosome coverage array.
                cov_vec = np.zeros([chrom_len]
                                   , dtype=int)

                # ---------------------------------------------------------------------- #
                # Step 4.5.1: join genes on reads data
                # so that each read is tied to a gene, for read counting purposes.
                # ---------------------------------------------------------------------- #

                # add IntervalIndex index to chromosome gene data. 0-index gene_starts, gene_ends because
                # reads are 0-indexed.
                chrom_gene_df[['gene_start', 'gene_end']] = chrom_gene_df[['gene_start', 'gene_end']].values - 1
                chrom_gene_df.index = IntervalIndex.from_arrays(chrom_gene_df.gene_start
                                                                , right=chrom_gene_df.gene_end
                                                                , closed='both')

                try:
                    reads_df['gene'] = chrom_gene_df.loc[reads_df.pos].gene.values

                # if there remains at least one read that doesn't land within a gene span,
                # need to run the interval join in SQL (TODO: diagnose why this happens)
                # as pandas doesn't offer a solution. See https://bit.ly/2Wg7Kjl.
                except KeyError:
                    # Make a sqlite db in memory.
                    conn = connect(':memory:')

                    # write the required data to SQL.
                    reads_df[['read_id', 'pos']].to_sql('reads'
                                                        , con=conn
                                                        , index=False)
                    chrom_gene_df[['gene_start', 'gene_end']].to_sql('genes'
                                                                     , con=conn
                                                                     , index=False)

                    # query read IDs for reads with position within [gene_start, gene_end].
                    qry = '''
                    select  
                        read_id
                    from
                        reads join genes on
                        pos between gene_start and gene_end
                    '''

                    # execute valid read ID query.
                    valid_read_ids = read_sql_query(qry
                                                    , con=conn).read_id
                    conn.close()

                    # subset reads to reads w/ valid read ID, then join with interval index again.
                    reads_df = reads_df[reads_df.read_id.isin(valid_read_ids)]
                    reads_df['gene'] = chrom_gene_df.loc[reads_df.pos].gene.values

                    del valid_read_ids

                # loop over reads for isolated genes, incrementing read count and coverage.
                dat = reads_df[['bounds', 'gene']].values
                for i in range(dat.shape[0]):
                    bounds, gene = dat[i, :]

                    # reads are already 0-indexed.
                    read_idx = fill_in_bounds(bounds
                                              , endpoint=True)

                    # increment coverage and read count.
                    cov_vec[read_idx] += 1
                    read_count_dict[gene] += 1

                # ---------------------------------------------------------------------- #
                # Step 4.5.2: save chromosome coverage vector.
                # chromosome overage vector ->> compressed csr numpy array
                # ---------------------------------------------------------------------- #
                chrom_cov_file = os.path.join(self.save_dir, 'chrom_coverage_' + self.sample_id + '_' + chrom + '.npz')

                if self.verbose:
                    logging.info('SAMPLE {0}, CHR {1} -- saving csr-compressed chrom coverage array.'
                                 .format(self.sample_id, chrom))

                # save coverage vector as a compressed-sparse row matrix.
                sparse.save_npz(chrom_cov_file
                                , matrix=sparse.csr_matrix(cov_vec))

                # drop large data objects.
                del cov_vec, dat, reads_df

            # drop remaining large data data objects.
            del chrom_gene_df, chrom_exon_df
            gc.collect()

            if self.verbose:
                logging.info('SAMPLE {0}, CHR {1} -- isolated gene reads processing successful.'
                             .format(self.sample_id, chrom))

        # ---------------------------------------------------------------------- #
        # Step 5. Save read counts.
        # chromosome read counts ->> .csv file
        # ---------------------------------------------------------------------- #
        count_file = os.path.join(self.save_dir, 'read_counts_' + self.sample_id + '_' + chrom + '.csv')

        # construct read count DataFrame from read count dictionary.
        read_count_df = DataFrame({'gene': list(read_count_dict.keys())
                                   , self.sample_id: list(read_count_dict.values())})

        del read_count_dict
        gc.collect()

        if self.verbose:
            logging.info('SAMPLE {0}, CHR {1} -- mean per-gene read count: {2:.4}'
                         .format(self.sample_id, chrom, read_count_df[self.sample_id].mean()))
            logging.info('SAMPLE {0}, CHR {1} -- saving read counts.'
                         .format(self.sample_id, chrom))

        # save sample's chromosome read counts to .csv for joining later.
        read_count_df.to_csv(count_file
                             , index=False)

    def coverage_read_counts(self, gene_overlap_dict, gene_df, exon_df):
        """
        Main function for computing coverage arrays in parallel over chromosomes.

        :param gene_overlap_dat: dictionary, keys are chromosomes, values are sub-dicts
         with output from gene_processing.get_gene_overlap_structure function.
        :param gene_df: pandas.DataFrame with `chr`, `gene`, `gene_start`, and `gene_end` columns
        that delineate the start and end position of a gene's transcript on a chromosome. See
        GeneAnnotationProcessor.
        :return: list of str file paths of compressed .npz files containing coverage arrays.
        """
        # create directory in DegNorm output dir where sample coverage vecs are saved.
        if not os.path.exists(self.save_dir):
            os.makedirs(self.save_dir)

        if self.verbose:
            logging.info('SAMPLE {0}: begin computing coverage, read counts for {1} chromosomes...'
                         .format(self.sample_id, len(self.chroms)))

        # distribute work across chromosomes with joblib.Parallel.
        out = Parallel(n_jobs=min(self.n_jobs, len(self.chroms))
                       , verbose=0
                       , backend='threading')(delayed(self.chromosome_coverage_read_counts)(
            gene_overlap_dat=gene_overlap_dict.get(chrom),
            chrom_gene_df=subset_to_chrom(gene_df, chrom=chrom),
            chrom_exon_df=subset_to_chrom(exon_df, chrom=chrom),
            chrom=chrom)
            for chrom in self.chroms)
