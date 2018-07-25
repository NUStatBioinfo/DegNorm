import re
from degnorm.utils import *
from degnorm.loaders import SamLoader, BamLoader


class ReadsCoverageProcessor():
    def __init__(self, sam_file=None, bam_file=None, n_jobs=max_cpu(), tmp_dir=None, verbose=True):
        """
        Genome coverage reader for a single RNA-seq experiment.
        Goal is to assemble a dictionary (chromosome, coverage array) pairs.

        :param sam_file: str .sam filename (optional if .bam file is specified)
        :param bam_file: str .bam filename (option if .sam file is specified)
        :param tmp_dir: str path to directory where coverage array files are saved.
        :param n_jobs: int number of CPUs to use for determining genome coverage. Default
        is number of CPUs on machine - 1
        :param verbose: bool indicator should progress be written to logger?
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
        if not tmp_dir:
            tmp_dir = os.path.join(os.path.dirname(self.filename), 'tmp')

        file_basename = '.'.join(os.path.basename(self.filename).split('.')[:-1])
        self.tmp_dir = os.path.join(tmp_dir, file_basename)
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.sample_id = file_basename
        self.loader = None

    def load(self):
        """
        Load a .sam or .bam file, obtain reads data and header.
        """
        if self.is_bam:
            self.loader = BamLoader(self.filename)
        else:
            self.loader = SamLoader(self.filename)

        df_dict = self.loader.get_data()
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

    def chromosome_coverage(self, df, chrom=None, chrom_len=0):
        """
        Determine per-chromosome reads coverage from an RNA-seq experiment. The cigar scores from
        single and paired reads are parsed according to _cigar_segment_bounds.

        Saves compressed coverage array to self.tmp_dir with file name 'sample_[sample_id]_[chrom].npz'

        :param df: pandas.DataFrame loaded .sam or .bam file for a single chromosome
        :param chrom: str chromosome name. If not supplied, presume df is already a chromosome subset
        :param chrom_len: int length of chromosome from reference genome
        directory where .sam file came from.
        :return: str full file path to where coverage array is saved in a compressed .npz file.
        """
        if chrom:
            df_chrom = subset_to_chrom(df, chrom=chrom)
        else:
            df_chrom = df

        # grab cigar string and read starting position. Initialize coverage array.
        dat = df_chrom[['cigar', 'pos']].values
        cov_vec = np.zeros([chrom_len])

        for i in np.arange(1, dat.shape[0], 2):
            bounds_1, _ = self._cigar_segment_bounds(dat[i - 1, 0], start=dat[i - 1, 1])
            bounds_2, _ = self._cigar_segment_bounds(dat[i, 0], start=dat[i, 1])

            # leverage nature of alignments of paired reads to find disjoint coverage ranges.
            min_bounds_1, max_bounds_1 = min(bounds_1), max(bounds_1)
            min_bounds_2, max_bounds_2 = min(bounds_2), max(bounds_2)

            if max_bounds_2 > max_bounds_1:
                bounds_2 = [max_bounds_1 + 1 if j <= max_bounds_1 else j for j in bounds_2]
            else:
                bounds_2 = list(set([min_bounds_1 - 1 if j >= min_bounds_1 else j for j in bounds_2]))
                bounds_2.sort()

            bounds = bounds_1 + bounds_2
            for j in np.arange(1, len(bounds), 2):
                cov_vec[(bounds[j - 1]) : (bounds[j] + 1)] += 1

        if self.verbose:
            logging.info('CHROMOSOME {0} -- length: {1}'.format(chrom, len(cov_vec)))
            logging.info('CHROMOSOME {0} -- max read coverage: {1}'.format(chrom, str(np.max(cov_vec))))
            logging.info('CHROMOSOME {0} -- mean read coverage: {1}'.format(chrom, str(np.mean(cov_vec))))
            logging.info('CHROMOSOME {0} -- % of chromosome covered: {1}'.format(chrom, str(np.mean(cov_vec > 0))))

        # create output directory if it does not exist, and then make output file name.
        out_file = os.path.join(self.tmp_dir, 'sample_' + self.sample_id + '_' + chrom + '.npz')
        if not os.path.isdir(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        logging.info('CHROMOSOME {0} -- saving coverage array to {1}'.format(chrom, out_file))
        np.savez_compressed(out_file
                            , cov=cov_vec)

        return out_file

    def coverage(self):
        """
        Main function for computing coverage arrays in parallel over chromosomes.

        :return: list of str file paths of compressed .npz files containing coverage arrays.
        """
        if self.verbose:
            logging.info('Loading file {0} into pandas.DataFrame'.format(self.filename))

        # load .sam or .bam file's reads + header.
        df, header_df = self.load()

        if self.verbose:
            logging.info('Successfully loaded {0} file. Total reads -- {1}'.format(
                '.bam' if self.is_bam else '.sam', df.shape[0]))

        # determine chromosomes whose coverage will be computed.
        chroms = df.chr.unique()
        header_df = header_df[header_df['chr'].isin(chroms)]

        if self.verbose:
            logging.info('Determining coverage for {0} chromosomes...\n'
                         '{1}'.format(len(chroms), ', '.join(chroms)))

        # run .chromosome_coverage in parallel over chromosomes.
        p = mp.Pool(processes=self.n_jobs)
        cov_filepaths = [p.apply_async(self.chromosome_coverage
                                       , args=(df, header_df.chr.iloc[i], header_df.length.iloc[i])) for i in
                         range(header_df.shape[0])]
        p.close()
        cov_filepaths = [x.get() for x in cov_filepaths]

        return cov_filepaths


# ==================================================================================================================== #
#                                                  Deprecated                                                          #
# ==================================================================================================================== #

# def chromosome_coverage_fortran(self, df, chrom=None, chrom_len=0):
    #     """
    #     Parse cigar strings from a loaded .sam or .bam file using self._cigar_segment_strings
    #     and dump them to a .txt file so that they can be parsed into coverage arrays in another application.
    #
    #     File will [path to .sam or .bam file]/[.sam or .bam filename]_[chromosome]_cigar.txt
    #
    #     Example structure of file contents:
    #
    #     chr7 3445092
    #     g1 104 204 110 210 5679 5685 0 0 0 0
    #     g2 543 3212 9801 9890 10101 10201 0 0
    #
    #     :param df: pandas.DataFrame loaded .sam or .bam file for a single chromosome
    #     :param chrom: str chromosome name. If not supplied, presume df is already a chromosome subset
    #     :param chrom_len: int length of chromosome from reference genome
    #     :param output_path: str path to output directory to dump .txt files. If None, use
    #     directory where .sam file came from.
    #     :return: None (write file to disk)
    #     """
    #     if chrom:
    #         df_chrom = subset_to_chrom(df, chrom=chrom)
    #     else:
    #         df_chrom = df
    #
    #     dat = df_chrom[['cigar', 'pos']].values
    #     n_reads = dat.shape[0] + 1
    #
    #     cig_bounds_list = list()
    #     for i in np.arange(1, n_reads, 2):
    #         bounds_1, n_match_1 = self._cigar_segment_bounds(dat[i - 1, 0], start=dat[i - 1, 1])
    #         bounds_2, n_match_2 = self._cigar_segment_bounds(dat[i, 0], start=dat[i, 1])
    #
    #         # leverage nature of alignments of paired reads to find disjoint coverage ranges.
    #         min_bounds_1, max_bounds_1 = min(bounds_1), max(bounds_1)
    #         min_bounds_2, max_bounds_2 = min(bounds_2), max(bounds_2)
    #
    #         if max_bounds_2 > max_bounds_1:
    #             bounds_2 = [max_bounds_1 + 1 if j <= max_bounds_1 else j for j in bounds_2]
    #         else:
    #             bounds_2 = list(set([min_bounds_1 - 1 if j >= min_bounds_1 else j for j in bounds_2]))
    #             bounds_2.sort()
    #
    #         cig_bounds_list.append(bounds_1 + bounds_2)
    #
    #     # determine maximum number of cigar match segments per read pair,
    #     # pad cigar strings with 0's for easing fortran IO
    #     max_col_width = max(map(lambda x: len(x), cig_bounds_list))
    #     for i in range(len(cig_bounds_list)):
    #         bounds = cig_bounds_list[i]
    #         n_bounds = len(bounds)
    #
    #         pad = ''
    #         if n_bounds < max_col_width:
    #             pad = ' ' + ' '.join(['0'] * (max_col_width - n_bounds))
    #
    #         cig_bounds_list[i] = 'g' + str(i + 1) + ' ' + ' '.join([str(j) for j in bounds]) + pad
    #
    #     # build output file name.
    #     file_basename = '.'.join(os.path.basename(self.filename).split('.')[:-1])
    #     output_filepath = os.path.join(self.tmp_dir
    #                                    , file_basename + '_' + chrom + '_cigar.txt')
    #     if self.verbose:
    #         logging.info('CHROMOSOME {0} ---- writing cigar segment strings to {1}'.format(chrom, output_filepath))
    #
    #     with open(output_filepath, 'w') as out:
    #         out.write(chrom + ' ' + str(chrom_len) + '\n')
    #         out.write('\n'.join(cig_bounds_list))

    # def coverage_fortran(self):
    #     """
    #     Main function for parsing cigar strings into start,end string segments and dumping them to disk.
    #     """
    #     if self.verbose:
    #         logging.info('Loading file {0} into pandas.DataFrame'.format(self.filename))
    #
    #     df, header_df = self.load()
    #
    #     if self.verbose:
    #         logging.info('Successfully loaded {0} file. Total reads: {1}'.format(
    #             '.bam' if self.is_bam else '.sam', df.shape[0]))
    #
    #     chroms = df.chr.unique()
    #     header_df = header_df[header_df['chr'].isin(chroms)]
    #
    #     if not os.path.isdir(self.tmp_dir):
    #         os.makedirs(self.tmp_dir)
    #
    #     if self.verbose:
    #         logging.info('Begin writing cigar segment strings for {0} chromosomes...'.format(len(chroms)))
    #
    #     p = mp.Pool(processes=self.n_jobs)
    #     out = [p.apply_async(self.chromosome_coverage_fortran
    #                          , args=(df, header_df.chr.iloc[i], header_df.length.iloc[i])) for i in
    #            range(header_df.shape[0])]
    #     p.close()