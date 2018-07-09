import re
from degnorm.utils import *
from degnorm.loaders import SamLoader, BamLoader


class ReadsCoverageParser():

    def __init__(self, sam_file=None, bam_file=None, n_jobs=max_cpu(), verbose=False):
        """
        Genome coverage reader for a single RNA-seq experiment.
        Goal is to assemble a dictionary (chromosome, coverage array) pairs.

        :param sam_file: chr .sam filename (optional if .bam file is specified)
        :param bam_file: chr .bam filename (option if .sam file is specified)
        :param n_jobs: int number of CPUs to use for determining genome coverage. Default
        is number of CPUs on machine - 1
        """
        if bam_file and sam_file:
            raise ValueError('cannot specify both a .sam and a .bam file')

        # determine if .sam or .bam file
        self.is_bam = False
        if bam_file:
            self.is_bam = True
            self.filename = bam_file

        else:
            self.filename = sam_file

        self.n_jobs = n_jobs
        self.verbose = verbose

    def load(self):
        if self.is_bam:
            df_dict = BamLoader(self.filename).get_data()
        else:
            df_dict = SamLoader(self.filename).get_data()

        return df_dict[self.filename]

    @staticmethod
    def _cigar_segments(cigar, start):
        """
        Determine the start and end positions on a chromosome of a non-no-matching part of an
        RNA-seq read based on a read's cigar string. All regions

        cigar string meanings: http://bioinformatics.cvr.ac.uk/blog/tag/cigar-string/

        :param cigar: chr a read's cigar string, e.g. "49M165N51M"
        :param start: int a read's start position on a chromosome
        :return: list of lists; each sub-list has two elements, the start and end position (end of match is inclusive)
        of a matching region of the read.
        """
        cigar_split = [(v, int(k)) for k, v in re.findall(r'(\d+)([A-Z][a-z]?)', cigar)]

        shift_one = True
        segments = list()
        for segment in cigar_split:

            extension = segment[1] + 1 if shift_one else segment[1]
            shift_one = False

            if segment[0] == 'M':
                segments.append([start, start + extension])

            start += extension

        return segments

    @staticmethod
    def _fill_segments(lst2d):
        """
        For each sub-list in a 2-d list structure (a list of lists), fill in
        the integers between the first and second elements of the sub-list.

        For example:
        [[4, 7], [9, 13]] -> [(4, 5, 6), (9, 10, 11, 12)]

        :param lst2d: list of lists; each sub-list should have two integers defining
        a start and an end position, respectively
        :return: list of numpy arrays; each array is a range spanning the integers
        between the first and second integer the sub-lists contained in lst2d
        """
        lst2d_expanded = [np.arange(lst1d[0], lst1d[1]) for lst1d in lst2d]

        return lst2d_expanded

    def cigar_processor(self, cigars, starts):
        """
        Determine the locations covered by a single read or set of two paired reads based on the
        cigar score(s) -- use for computing complete reads coverage.

        For a single read or set of two paired reads, take the unique combined set of
        matched regions of a chromosome.

        :param cigars: list or numpy array of str cigar strings
        :param starts: list or numpy array of int read start positions. Must have same lengths as `cigars`.
        :return: 1-d numpy array of int positions covered by the
        """
        n = len(cigars)
        if n > 2:
            raise ValueError('len(cigars) > 2 ...... cigars = {0}'.format(', '.join(cigars)))

        for idx in range(len(cigars)):
            append_me = flatten_2d(self._fill_segments(self._cigar_segments(cigars[idx]
                                                                            , start=starts[idx])))
            if idx == 0:
                read_idx_arr = append_me
            else:
                read_idx_arr = np.append(read_idx_arr, append_me)

        return np.unique(read_idx_arr)

    def chromosome_reads_coverage(self, df, chrom=None):
        """
        Determine per-chromosome reads coverage from an RNA-seq experiment. The cigar scores from
        single and paired reads are parsed according to cigar_processor

        :param df: pandas.DataFrame loaded .sam file for a single chromosome
        :param chrom: str chromosome name. If not supplied, presume df is already chromosome-subsetted.
        :return: tuple (chrom, 1-d numpy array -- entire chromosome's reads coverage)
        """
        if chrom:
            df_chrom = subset_to_chrom(df, chrom=chrom)
        else:
            df_chrom = df

        coverage_arr = np.zeros(df_chrom.end_pos.max() + 1)

        gb = df_chrom.groupby(['qname_unpaired'])
        for grp in gb:
            covered_idx = self.cigar_processor(grp[1].cigar.values
                                               , starts=grp[1].pos.values)
            coverage_arr[covered_idx] += 1

        if self.verbose:
            logging.info('---- CHROMOSOME {0} ----'.format(chrom))
            logging.info('Max read coverage -- {0}'.format(str(np.max(coverage_arr))))
            logging.info('Mean read coverage -- {0}'.format(str(np.mean(coverage_arr))))
            logging.info('% of chromosome covered -- {0}'.format(str(np.mean(coverage_arr > 0))))

        return (chrom, coverage_arr)

    def coverage(self):
        if self.verbose:
            logging.info('Loading file {0} into pandas.DataFrame'.format(self.filename))

        df = self.load()

        # hack: drop nonunique paired reads. TODO: look into why there are nonunique qname pairs.
        df.drop_duplicates(subset=['chr', 'qname']
                           , inplace=True)

        if self.verbose:
            logging.info('Successfully loaded {0} file. Shape -- {1}'.format(
                '.bam' if self.is_bam else '.sam', df.shape))

        chroms = df.chr.unique()

        if self.verbose:
            logging.info('Determining position coverage for {0} chromosomes...'.format(len(chroms)))

        p = mp.Pool(processes=self.n_jobs)
        coverages = [p.apply(self.chromosome_reads_coverage, args=(df, chrom)) for chrom in chroms]
        p.close()

        chrom_cov_dict = dict()
        for tup in coverages:
            chrom_cov_dict[tup[0]] = tup[1]

        return chrom_cov_dict

