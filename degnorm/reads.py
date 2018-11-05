from pandas import DataFrame, concat
from degnorm.utils import *
from degnorm.loaders import BamLoader
from joblib import Parallel, delayed
from scipy import sparse


def cigar_segment_bounds(cigar, start):
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


class BamReadsProcessor():
    def __init__(self, bam_file, index_file, chroms=None, n_jobs=max_cpu(),
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
        self.filename = bam_file

        # determine where to dump coverage .npz files.
        if not output_dir:
            output_dir = os.path.join(os.path.dirname(self.filename), 'tmp')

        file_basename = '.'.join(os.path.basename(self.filename).split('.')[:-1])
        self.save_dir = os.path.join(output_dir, file_basename)
        self.index_filename = index_file
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.sample_id = file_basename
        self.header = None
        self.paired = None
        self.chroms = chroms
        self.unique_alignment = unique_alignment
        self.loader = BamLoader(self.filename, self.index_filename)
        self.get_header()
        self.determine_if_paired()

    def get_header(self):

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
            self.chroms = np.intersect1d(self.chroms, self.header.chr.unique())

        # if no chromosomes specified, assume load includes every chromosome in header.
        else:
            self.chroms = self.header.chr.unique()

    def determine_if_paired(self):

        self.paired = False
        bam_file = self.loader.get_data()

        # pull first 300 queries' query names.
        ctr = 0
        qnames = list()
        for read in bam_file.fetch(self.header.chr.iloc[0]):
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

        read_attributes = ['query_name', 'pos', 'cigarstring']
        reads = list()
        bam_file = self.loader.get_data()

        for read in bam_file.fetch(chrom):

            # if working only with unique alignment reads, skip read if NH tag is > 1.
            if self.unique_alignment:
                if read.has_tag('NH'):
                    if read.get_tag('NH') > 1:
                        continue

            # if reading paired reads and the read is paired,
            # or if reading single end reads and the read is not paired,
            # grab the attributes of interest from the read.
            if ((read.is_paired) and (self.paired)) or ((not read.is_paired) and (not self.paired)):
                reads.append([getattr(read, attr) for attr in read_attributes])

        # close .bam file connection.
        bam_file.close()

        # transform .sam lines into pandas.DataFrame
        df = DataFrame(reads
                       , columns=['qname', 'pos', 'cigar'])

        # remove reads list (large memory).
        del reads
        gc.collect()

        # if working with paired data: sort reads by unpaired query name.
        if self.paired:
            df['qname_unpaired'] = df.qname.apply(lambda x: '.'.join(x.split('.')[:-1]))
            df.sort_values('qname_unpaired', inplace=True)

        return df

    def chromosome_coverage_read_counts(self, gene_sub_df, chrom):
        """
        Determine per-chromosome reads coverage and per-gene read counts from an RNA-seq experiment.
        The cigar scores from single and paired reads are parsed according to cigar_segment_bounds.

        Saves compressed coverage array to self.save_dir with file name 'sample_[sample_id]_[chrom].npz'

        :param gene_df: pandas.DataFrame with `chr`, `gene`, `gene_start`, and `gene_end` columns
        that delineate the start and end position of a gene's transcript on a chromosome, must be
        subset to the chromosome in study.
        :param chrom: str chromosome name
        :return: str full file path to where coverage array is saved in a compressed .npz file.
        """

        # First, parse this chromosome's reads.
        reads_df = self.load_chromosome_reads(chrom)

        # Second, if working with paired reads,
        # ensure that we've sequestered paired reads (eliminate any query names only occurring once).
        if self.paired:
            qname_counts = reads_df.qname_unpaired.value_counts()
            paired_occ_reads = qname_counts[qname_counts == 2].index.values.tolist()
            reads_df = reads_df[reads_df.qname_unpaired.isin(paired_occ_reads)]

        # grab cigar string and read starting position. Initialize coverage array.
        dat = reads_df[['cigar', 'pos']].values
        cov_vec = np.zeros([self.header[self.header.chr == chrom].length[0]])

        # for paired reads, perform special parsing of CIGAR strings to avoid double-counting of overlap regions.
        if self.paired:
            for i in np.arange(1, dat.shape[0], 2):
                bounds_1 = cigar_segment_bounds(dat[i - 1, 0], start=dat[i - 1, 1])
                bounds_2 = cigar_segment_bounds(dat[i, 0], start=dat[i, 1])

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
            for i in np.arange(dat.shape[0]):
                bounds = cigar_segment_bounds(dat[i, 0], start=dat[i, 1])

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
            counts_df = reads_df[reads_df.pos.between(dat[i, 0], dat[i, 1])]
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

        # return per-chromosome coverage array filename and per-chromosome read counts.
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
            gene_sub_df=subset_to_chrom(gene_df, chrom=chrom),
            chrom=chrom)
            for chrom in self.chroms)

        # parse output from parallel workers.
        cov_filepaths = [x[0] for x in par_output]
        read_count_dfs = [x[1] for x in par_output]

        return cov_filepaths, concat(read_count_dfs)


if __name__ == '__main__':
    # bam_file = '/Users/fineiskid/Downloads/ff_small.bam'
    # bai_file = '/Users/fineiskid/Downloads/ff_small.bai'
    # bam_file = '/Users/fineiskid/Downloads/hg_small_1.bam'
    # bai_file = '/Users/fineiskid/Downloads/hg_small_1.bai'
    bam_file = '/Users/fineiskid/Downloads/hg_small_2.bam'
    bai_file = '/Users/fineiskid/Downloads/hg_small_2.bai'

    from degnorm.gene_processing import GeneAnnotationProcessor
    gap = GeneAnnotationProcessor('/Users/fineiskid/Downloads/chr1_small.gtf'
                                  , n_jobs=1
                                  , verbose=True)
    exon_df = gap.run()
    genes_df = exon_df[['chr', 'gene', 'gene_start', 'gene_end']].drop_duplicates().reset_index(drop=True)

    bammer = BamReadsProcessor(bam_file=bam_file
                               , index_file=bai_file)
    dat = bammer.coverage_read_counts(genes_df)
    print(dat[0])
    print(dat[1].head())
    print(bammer.paired)