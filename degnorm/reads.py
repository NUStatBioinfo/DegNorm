import re
from pandas import DataFrame, concat
from degnorm.utils import *
from degnorm.loaders import SamLoader
from joblib import Parallel, delayed
from scipy import sparse


class ReadsProcessor():
    def __init__(self, sam_file=None, chroms=None, n_jobs=max_cpu(),
                 output_dir=None, unique_alignment=False, verbose=True):
        """
        Genome coverage reader for a single RNA-seq experiment, contained in a .sam file.
        Goal is to assemble a dictionary (chromosome, coverage array) pairs.

        :param sam_file: str .sam filename
        :param output_dir: str path to DegNorm output directory where coverage array files will be saved.
        If not specified, will use directory where RNA Seq experiment file is located.
        :param chroms: list of str names of chromosomes to load.
        :param n_jobs: int number of CPUs to use for determining genome coverage. Default
        is number of CPUs on machine - 1.
        :param unique_alignment: bool indicator - drop reads with NH:i:<x> flag where x > 1.
        :param verbose: bool indicator should progress be written to logger?
        """
        self.filename = sam_file

        # determine where to dump coverage .npz files.
        if not output_dir:
            output_dir = os.path.join(os.path.dirname(self.filename), 'tmp')

        file_basename = '.'.join(os.path.basename(self.filename).split('.')[:-1])
        self.save_dir = os.path.join(output_dir, file_basename)
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.sample_id = file_basename
        self.loader = None
        self.data = None
        self.header = None
        self.paired = None
        self.chroms = chroms
        self.unique_alignment = unique_alignment

    def load(self):
        """
        Load a .sam, obtain reads data and header.
        """
        self.loader = SamLoader(self.filename)

        # load up .sam file.
        df_dict = self.loader.get_data(chrom=self.chroms
                                       , unique_alignment=self.unique_alignment)
        df = df_dict['data']

        # .sam file must have valid header.
        try:
            header_df = df_dict['header']

        except KeyError:
            logging.error('No header was found in {0}!'.format(self.filename))
            raise

        # hack: drop nonunique paired reads. TODO: look into why there are nonunique (chr, qname) pairs.
        df.drop_duplicates(subset=['chr', 'qname']
                           , inplace=True)

        self.data = df
        self.header = header_df
        self.paired = df_dict['paired']

    @staticmethod
    def _cigar_segment_bounds(cigar, start):
        """
        Determine the start and end positions on a chromosome of a non-no-matching part of an
        RNA-seq read based on a read's cigar string.

        cigar string meaning: http://bioinformatics.cvr.ac.uk/blog/tag/cigar-string/

        :param cigar: str a read's cigar string, e.g. "49M165N51M"
        :param start: int a read's start position on a chromosome
        :return: list of integers representing cigar match start, end points, e.g.
            50M25N50M starting from 100 -> [100, 149, 175, 224]. Note that start and end integers
            are inclusive, i.e. all positions at or between 100 and 149 and at or between 175 and 224
            are covered by reads.
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

    def chromosome_coverage_read_counts(self, reads_sub_df, gene_sub_df, chrom, chrom_len=0):
        """
        Determine per-chromosome reads coverage and per-gene read counts from an RNA-seq experiment.
        The cigar scores from single and paired reads are parsed according to _cigar_segment_bounds.

        Saves compressed coverage array to self.save_dir with file name 'sample_[sample_id]_[chrom].npz'

        :param reads_sub_df: pandas.DataFrame of ordered paired reads data with `cigar` and `pos` fields,
        must be subset to the chromosome in study.
        :param gene_df: pandas.DataFrame with `chr`, `gene`, `gene_start`, and `gene_end` columns
        that delineate the start and end position of a gene's transcript on a chromosome, must be
        subset to the chromosome in study.
        :param chrom: str chromosome name
        :param chrom_len: int length of chromosome from reference genome
        :return: str full file path to where coverage array is saved in a compressed .npz file.
        """
        # first, if working with paired reads,
        # ensure that we've sequestered paired reads (eliminate any query names only occurring once).
        if self.paired:
            qname_counts = reads_sub_df.qname_unpaired.value_counts()
            paired_occ_reads = qname_counts[qname_counts == 2].index.values.tolist()
            reads_sub_df = reads_sub_df[reads_sub_df.qname_unpaired.isin(paired_occ_reads)]

        # grab cigar string and read starting position. Initialize coverage array.
        dat = reads_sub_df[['cigar', 'pos']].values
        cov_vec = np.zeros([chrom_len])

        # for paired reads, perform special parsing of CIGAR strings to avoid double-counting of overlap regions.
        if self.paired:
            for i in np.arange(1, dat.shape[0], 2):
                bounds_1 = self._cigar_segment_bounds(dat[i - 1, 0], start=dat[i - 1, 1])
                bounds_2 = self._cigar_segment_bounds(dat[i, 0], start=dat[i, 1])

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
                    cov_vec[(bounds[j - 1]):(bounds[j] + 1)] += 1

        # for single-read RNA-Seq experiments, we do not need such special consideration.
        else:
            for i in np.arange(0, dat.shape[0]):
                bounds = self._cigar_segment_bounds(dat[i - 1, 0], start=dat[i - 1, 1])

                for j in np.arange(1, len(bounds), 2):
                    cov_vec[(bounds[j - 1]):(bounds[j] + 1)] += 1

        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} length: {2}'.format(self.sample_id, chrom, len(cov_vec)))
            logging.info('SAMPLE {0}: CHROMOSOME {1} max read coverage: {2:.4}'
                         .format(self.sample_id, chrom, str(np.max(cov_vec))))
            logging.info('SAMPLE {0}: CHROMOSOME {1} mean read coverage: {2:.4}'
                         .format(self.sample_id, chrom, str(np.mean(cov_vec))))
            logging.info('SAMPLE {0}: CHROMOSOME {1} % of chromosome covered: {2:.4}'
                         .format(self.sample_id, chrom, str(np.mean(cov_vec > 0))))

        # create output directory if it does not exist, and then make output file name.
        out_file = os.path.join(self.save_dir, 'sample_' + self.sample_id + '_' + chrom + '.npz')

        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} saving csr-compressed coverage array to {2}'
                         .format(self.sample_id, chrom, out_file))

        # save coverage vector as a compressed-sparse row matrix.
        sparse.save_npz(out_file
                        , matrix=sparse.csr_matrix(cov_vec))

        # free up memory, delete coverage vector (it's saved to disk).
        del cov_vec
        gc.collect()

        # finally, parse gene read counts.
        n_genes = gene_sub_df.shape[0]
        counts = np.zeros(n_genes)
        dat = gene_sub_df[['gene_start', 'gene_end']].values

        # iterate over genes, count number of reads (entirely) falling between gene_start and gene_end.
        qname_col = 'qname_unpaired' if self.paired else 'qname'
        for i in range(n_genes):
            # entire read must fall w/in gene start/end
            counts_df = reads_sub_df[reads_sub_df.pos.between(dat[i, 0], dat[i, 1])]
            counts[i] = counts_df[qname_col].unique().shape[0]

        # turn read counts into a DataFrame so we can join on genes later.
        read_count_df = DataFrame({'chr': chrom
                                      , 'gene': gene_sub_df.gene.values
                                      , 'read_count': counts})

        # quality control.
        if read_count_df.empty:
            raise ValueError('Missing read counts!')

        if self.verbose:
            logging.info('SAMPLE {0}: CHROMOSOME {1} mean per-gene read count: {2:.4}'
                         .format(self.sample_id, chrom, read_count_df.read_count.mean()))

        return out_file, read_count_df

    def coverage_read_counts(self, gene_df):
        """
        Main function for computing coverage arrays in parallel over chromosomes.

        :param gene_df: pandas.DataFrame with `chr`, `gene`, `gene_start`, and `gene_end` columns
        that delineate the start and end position of a gene's transcript on a chromosome. See
        GeneAnnotationProcessor.
        :return: list of str file paths of compressed .npz files containing coverage arrays.
        """
        if self.verbose:
            logging.info('Begin reading file {0}...'.format(self.filename))

        # load .sam file's reads + header.
        self.load()

        if self.verbose:
            logging.info('Successfully read file {0}. Total transcript reads -- {1}'
                         .format(self.filename, self.data.shape[0]))

        # determine chromosomes whose coverage will be computed.
        chroms = self.data.chr.unique()
        header_df = self.header[self.header.chr.isin(chroms)]

        # create directory in DegNorm output dir where sample coverage vecs are saved.
        if not os.path.exists(self.save_dir):
            os.makedirs(self.save_dir)

        if self.verbose:
            logging.info('SAMPLE {0}: determining read coverage and read counts for {1} chromosomes.\n'
                         'Saving output to directory {2}'
                         .format(self.sample_id, len(chroms), self.save_dir))

        # distribute work with joblib.Parallel:
        par_output = Parallel(n_jobs=min(self.n_jobs, len(chroms))
                              , verbose=0
                              , backend='threading')(delayed(self.chromosome_coverage_read_counts)(
            reads_sub_df=subset_to_chrom(self.data, chrom=chrom),
            gene_sub_df=subset_to_chrom(gene_df, chrom=chrom),
            chrom=chrom,
            chrom_len=header_df[header_df.chr == chrom].length.iloc[0])
            for chrom in chroms)

        # parse output from parallel workers.
        # par_output = [x.get() for x in par_output]
        cov_filepaths = [x[0] for x in par_output]
        read_count_dfs = [x[1] for x in par_output]

        return cov_filepaths, concat(read_count_dfs)