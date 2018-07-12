import re
from degnorm.utils import *
from degnorm.loaders import SamLoader, BamLoader


class ReadsCoverageParser():
    def __init__(self, sam_file=None, bam_file=None, n_jobs=max_cpu(), tmp_dir=None, verbose=False):
        """
        Genome coverage reader for a single RNA-seq experiment.
        Goal is to assemble a dictionary (chromosome, coverage array) pairs.

        :param sam_file: chr .sam filename (optional if .bam file is specified)
        :param bam_file: chr .bam filename (option if .sam file is specified)
        :param tmp_dir: chr path to directory where cigar parse .txt files will be dumped for Fortran IO.
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

        # determine where to dump .txt files from cigar string parse.
        if tmp_dir:
            self.tmp_dir = tmp_dir

        else:
            self.tmp_dir = os.path.join(os.path.dirname(self.filename), 'tmp')

        self.n_jobs = n_jobs
        self.verbose = verbose

    def load(self):
        """
        Load a .sam or .bam file, obtain reads data and header.
        """
        if self.is_bam:
            df_dict = BamLoader(self.filename).get_data()
        else:
            df_dict = SamLoader(self.filename).get_data()

        df = df_dict['data']
        header_df = df_dict['header']

        # hack: drop nonunique paired reads. TODO: look into why there are nonunique qname pairs.
        df.drop_duplicates(subset=['chr', 'qname']
                           , inplace=True)

        return df, header_df

    @staticmethod
    def _cigar_segment_bounds(cigar, start):
        """
        Determine the start and end positions on a chromosome of a non-no-matching part of an
        RNA-seq read based on a read's cigar string.

        cigar string meaning: http://bioinformatics.cvr.ac.uk/blog/tag/cigar-string/

        :param cigar: str a read's cigar string, e.g. "49M165N51M"
        :param start: int a read's start position on a chromosome
        :return: tuple
            list: a list of integers representing cigar match start,end points, e.g.
            50M25N50M starting from 100 -> [100, 149, 175, 224]. Note that start and end integers
            are inclusive, i.e. all positions at or between 100 and 149 and at or between 175 and 224
            are covered by reads.
            int: number of matching segments found within the cigar string. E.g. "100M" is one matching segment.
        """
        if cigar == '100M':
            return [start, start + 99], 1

        cigar_split = [(v, int(k)) for k, v in re.findall(r'(\d+)([A-Z]?)', cigar)]

        num_match_segs = 0
        match_idx_list = list()

        for idx in range(len(cigar_split)):
            segment = cigar_split[idx]

            if segment[0] == 'M':
                extension = segment[1] - 1
                augment = True
                match_idx_list += [start, start + extension]
                num_match_segs += 1

            else:
                if augment:
                    extension = segment[1] + 1
                    augment = False
                else:
                    extension = segment[1]

            start += extension

        return match_idx_list, num_match_segs

    def chromosome_cigar_segment_dump(self, df, chrom=None, chrom_len=0):
        """
        Parse cigar strings from a loaded .sam or .bam file using self._cigar_segment_strings
        and dump them to a .txt file so that they can be parsed into coverage arrays in another application.

        File will [path to .sam or .bam file]/[.sam or .bam filename]_[chromosome]_cigar.txt

        Example structure of file contents:

        chr7 3445092
        g1 104 204 110 210 5679 5685 0 0 0 0
        g2 543 3212 9801 9890 10101 10201 0 0

        :param df: pandas.DataFrame loaded .sam or .bam file for a single chromosome
        :param chrom: str chromosome name. If not supplied, presume df is already a chromosome subset
        :param chrom_len: int length of chromosome from reference genome
        :param output_path: str path to output directory to dump .txt files. If None, use
        directory where .sam file came from.
        :return: None (write file to disk)
        """
        if chrom:
            df_chrom = subset_to_chrom(df, chrom=chrom)
        else:
            df_chrom = df

        dat = df_chrom[['cigar', 'pos']].values
        n_reads = dat.shape[0] + 1

        cig_bounds_list = list()
        for i in np.arange(1, n_reads, 2):
            bounds_1, n_match_1 = self._cigar_segment_bounds(dat[i - 1, 0], start=dat[i - 1, 1])
            bounds_2, n_match_2 = self._cigar_segment_bounds(dat[i, 0], start=dat[i, 1])

            # leverage nature of alignments of paired reads to find disjoint coverage ranges.
            min_bounds_1, max_bounds_1 = min(bounds_1), max(bounds_1)
            min_bounds_2, max_bounds_2 = min(bounds_2), max(bounds_2)

            if max_bounds_2 > max_bounds_1:
                bounds_2 = [max_bounds_1 + 1 if j <= max_bounds_1 else j for j in bounds_2]
            else:
                bounds_2 = list(set([min_bounds_1 - 1 if j >= min_bounds_1 else j for j in bounds_2]))
                bounds_2.sort()

            cig_bounds_list.append(bounds_1 + bounds_2)

        # determine maximum number of cigar match segments per read pair,
        # pad cigar strings with 0's for easing fortran IO
        max_col_width = max(map(lambda x: len(x), cig_bounds_list))
        for i in range(len(cig_bounds_list)):
            bounds = cig_bounds_list[i]
            n_bounds = len(bounds)

            pad = ''
            if n_bounds < max_col_width:
                pad = ' ' + ' '.join(['0'] * (max_col_width - n_bounds))

            cig_bounds_list[i] = 'g' + str(i + 1) + ' ' + ' '.join([str(j) for j in bounds]) + pad

        # build output file name.
        file_basename = '.'.join(os.path.basename(self.filename).split('.')[:-1])
        output_filepath = os.path.join(self.tmp_dir
                                       , file_basename + '_' + chrom + '_cigar.txt')
        if self.verbose:
            logging.info('CHROMOSOME {0} ---- writing cigar segment strings to {1}'.format(chrom, output_filepath))

        with open(output_filepath, 'w') as out:
            out.write(chrom + ' ' + str(chrom_len) + '\n')
            out.write('\n'.join(cig_bounds_list))

    def dump_cigar_parse(self):
        """
        Main function for parsing cigar strings into start,end string segments and dumping them to disk.
        """
        if self.verbose:
            logging.info('Loading file {0} into pandas.DataFrame'.format(self.filename))

        df, header_df = self.load()

        if self.verbose:
            logging.info('Successfully loaded {0} file. Total reads: {1}'.format(
                '.bam' if self.is_bam else '.sam', df.shape[0]))

        chroms = df.chr.unique()
        header_df = header_df[header_df['chr'].isin(chroms)]

        if not os.path.isdir(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        if self.verbose:
            logging.info('Begin writing cigar segment strings for {0} chromosomes...'.format(len(chroms)))

        p = mp.Pool(processes=self.n_jobs)
        out = [p.apply_async(self.chromosome_cigar_segment_dump
                             , args=(df, header_df.chr.iloc[i], header_df.length.iloc[i])) for i in
               range(header_df.shape[0])]
        p.close()


    # --------------------------------------------------------------------------------------------------- #
    #                                           DEPRECATED                                                #
    # --------------------------------------------------------------------------------------------------- #

    # @staticmethod
    # def _cigar_segments(cigar, start):
    #     """
    #     Determine the start and end positions on a chromosome of a non-no-matching part of an
    #     RNA-seq read based on a read's cigar string. All regions
    #
    #     cigar string meaning: http://bioinformatics.cvr.ac.uk/blog/tag/cigar-string/
    #
    #     :param cigar: str a read's cigar string, e.g. "49M165N51M"
    #     :param start: int a read's start position on a chromosome
    #     :return: list of lists; each sub-list has two elements, the start and end position (end of match is inclusive)
    #     of a matching region of the read.
    #     """
    #     if cigar == '100M':
    #         return [start, start + 101]
    #
    #     cigar_split = [(v, int(k)) for k, v in re.findall(r'(\d+)([A-Z]?)', cigar)]
    #
    #     shift_one = True
    #     segments = list()
    #     for segment in cigar_split:
    #
    #         extension = segment[1] + 1 if shift_one else segment[1]
    #         shift_one = False
    #
    #         if segment[0] == 'M':
    #             segments.append([start, start + extension])
    #
    #         start += extension
    #
    #     return segments

    # @staticmethod
    # def _fill_segments(lst2d):
    #     """
    #     For each sub-list in a 2-d list structure (a list of lists), fill in
    #     the integers between the first and second elements of the sub-list.
    #
    #     For example:
    #     [[4, 7], [9, 13]] -> [(4, 5, 6), (9, 10, 11, 12)]
    #
    #     :param lst2d: list of lists; each sub-list should have two integers defining
    #     a start and an end position, respectively
    #     :return: list of numpy arrays; each array is a range spanning the integers
    #     between the first and second integer the sub-lists contained in lst2d
    #     """
    #     lst2d_expanded = [np.arange(lst1d[0], lst1d[1]) for lst1d in lst2d]
    #
    #     return lst2d_expanded

    # def cigar_segment_processor(self, cigars, starts):
    #     """
    #     Determine the locations covered by a single read or set of two paired reads based on the
    #     cigar score(s) -- use for computing complete reads coverage.
    #
    #     For a single read or set of two paired reads, take the unique combined set of
    #     matched regions of a chromosome.
    #
    #     :param cigars: list or numpy array of str cigar strings
    #     :param starts: list or numpy array of int read start positions. Must have same lengths as `cigars`.
    #     :return: 1-d numpy array of int positions covered by the
    #     """
    #     n = len(cigars)
    #
    #     for idx in range(n):
    #         append_me = flatten_2d(self._fill_segments(self._cigar_segments(cigars[idx]
    #                                                                         , start=starts[idx])))
    #         if idx == 0:
    #             read_idx_arr = append_me
    #         else:
    #             read_idx_arr = np.append(read_idx_arr, append_me)
    #
    #     return np.unique(read_idx_arr)

    # def chromosome_reads_coverage(self, df, chrom=None):
    #     """
    #     Determine per-chromosome reads coverage from an RNA-seq experiment. The cigar scores from
    #     single and paired reads are parsed according to cigar_processor
    #
    #     :param df: pandas.DataFrame loaded .sam or .bam file for a single chromosome
    #     :param chrom: str chromosome name. If not supplied, presume df is already chromosome-subsetted.
    #     :return: tuple (chrom, 1-d numpy array -- entire chromosome's reads coverage)
    #     """
    #     if chrom:
    #         df_chrom = subset_to_chrom(df, chrom=chrom)
    #     else:
    #         df_chrom = df
    #
    #     coverage_arr = np.zeros(df_chrom.end_pos.max() + 1)
    #
    #     gb = df_chrom.groupby(['qname_unpaired'])
    #     for grp in gb:
    #         covered_idx = self.cigar_segment_processor(grp[1].cigar.values
    #                                                    , starts=grp[1].pos.values)
    #         coverage_arr[covered_idx] += 1
    #
    #     if self.verbose:
    #         logging.info('---- CHROMOSOME {0} ----'.format(chrom))
    #         logging.info('Max read coverage -- {0}'.format(str(np.max(coverage_arr))))
    #         logging.info('Mean read coverage -- {0}'.format(str(np.mean(coverage_arr))))
    #         logging.info('% of chromosome covered -- {0}'.format(str(np.mean(coverage_arr > 0))))
    #
    #     return (chrom, coverage_arr)
    #
    # def coverage(self):
    #     """
    #     Main function for computing coverage arrays in parallel over chromosomes.
    #     """
    #     if self.verbose:
    #         logging.info('Loading file {0} into pandas.DataFrame'.format(self.filename))
    #
    #     df, _ = self.load()
    #
    #     if self.verbose:
    #         logging.info('Successfully loaded {0} file. Total reads -- {1}'.format(
    #             '.bam' if self.is_bam else '.sam', df.shape[0]))
    #
    #     chroms = df.chr.unique()
    #
    #     if self.verbose:
    #         logging.info('Determining position coverage for {0} chromosomes...'.format(len(chroms)))
    #
    #     p = mp.Pool(processes=self.n_jobs)
    #     cov_tups = [p.apply_async(self.chromosome_reads_coverage, args=(df, chrom)) for chrom in chroms]
    #     p.close()
    #     cov_tups = [x.get() for x in cov_tups]
    #
    #     chrom_cov_dict = dict()
    #     for tup in cov_tups:
    #         chrom_cov_dict[tup[0]] = tup[1]
    #
    #     return chrom_cov_dict
