import os
import re
from pandas import DataFrame, read_table


class Loader():
    def __init__(self, filetypes):
        """
        Generic file loader class.

        :param filetypes: str or list of file extensions, only load files that .endswith(filetype). Using a list
        would be appropriate if there are multiple allowable file formats, e.g. .gtf and .gff files, because
        they share mostly the same file format.
        """
        self.filetypes = filetypes if isinstance(filetypes, list) else [filetypes]
        self.filename = None

    def get_file(self, to_load):
        """
        Obtain list of str realpaths to files to load.

        :param to_load: str the realpath to a file ending with self.filetype file types.
        """
        if isinstance(to_load, str):
            if os.path.exists(to_load):

                loadable = False
                for ft in self.filetypes:
                    if to_load.endswith(ft):
                        self.filename = to_load
                        loadable = True

                if not loadable:
                    raise ValueError('file {0} does not end with {0}'.format(to_load, ', '.join(self.filetypes)))
            else:
                raise IOError('file {0} not found'.format(to_load))

        else:
            raise ValueError('{0} data type not understood'.format(to_load))

    def get_data(self):
        raise NotImplementedError('get_data not yet implemented for {0}'.format(self.__class__.__name__))


class SamLoader(Loader):

    def __init__(self, to_load):
        """
        .sam file loader

        :param to_load: str the realpath to a .sam file.
        """
        Loader.__init__(self, '.sam')
        self.get_file(to_load)

    def get_data(self):
        """
        Extract .sam files listed in self.to_load into a dictionary of pandas.DataFrames.
        Each DataFrame has the following fields:

        - qname
        - chr
        - pos
        - cigar
        - pnext
        - tlen
        - qname_unpaired: qname field but without '.1' or '.2' extension indicating the
        paired-order of the read in a paired read
        - end_pos: rightmost position of aligned reads

        See .sam file specification: http://samtools.github.io/hts-specs/SAMv1.pdf

        :return: dictionary {.sam file: pandas.DataFrame}
        """
        df_dict = dict()

        with open(self.filename, 'r') as f:
            lines = f.readlines()

        # identify line number where header lines stop
        header_idx = 0
        is_header_line = lines[header_idx][0] == '@'

        # store header info: chromosome lengths
        chrom_len_dict = dict()

        while is_header_line:
            header_idx += 1
            line = lines[header_idx]
            is_header_line = line[0] == '@'

            if line[0:3] == '@SQ':

                for split in line.split('\t'):
                    content = split.split(':')[-1]

                    if split[0:2] == 'SN':
                        chrom = content
                    elif split[0:2] == 'LN':
                        chrom_len = int(content)

                chrom_len_dict[chrom] = chrom_len

        # turn header into a pandas.DataFrame
        header_df = DataFrame(list(chrom_len_dict.items())
                              , columns = ['chr', 'length'])

        if not header_df.empty:
            df_dict['header'] = header_df

        # extract relevant columns from .sam file lines.
        reads = list()
        for x in lines[header_idx:]:
            splt = x.split('\t')
            reads.append([splt[i] for i in [0, 2, 3, 5, 6]])

        # transform .sam lines into pandas.DataFrame
        colnames = ['qname', 'chr', 'pos', 'cigar', 'rnext']
        df = DataFrame(reads
                       , columns=colnames)

        # subset .sam file to paired reads using rnext column.
        df = df[df['rnext'] == '=']

        # typecasting: change int fields from str.
        int_cols = ['pos']
        for col in int_cols:
            df[col] = df[col].astype('int')

        # preprocessing: un-specify paired reads together so that a pair of alignments
        # share the same un-paired QNAME
        df['qname_unpaired'] = df.qname.apply(lambda x: '.'.join(x.split('.')[:-1]))

        # sort the reads so that pairs are grouped together.
        df.sort_values('qname_unpaired', inplace=True)

        df_dict['data'] = df.reset_index(drop=True)

        return df_dict


class GeneAnnotationLoader(Loader):

    def __init__(self, to_load):
        """
        .gtf or .gff file loader

        More about .gtf or .gff fields: https://useast.ensembl.org/info/website/upload/gff.html

        :param to_load: str the realpath to a .gtf file.
        """
        Loader.__init__(self, ['.gtf', '.gff'])
        self.get_file(to_load)

    @staticmethod
    def _attribute_to_gene(attribute, exprs):
        """
        Parse a .gtf/.gff attribute string for a gene_id or gene_name.

        For example:

        self._attribute_to_gene('gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018";'
                                , exprs = [re.compile('gene_id')])
        'DDX11L1'

        :param attribute: str attribute string from .gtf/.gff file -- "A semicolon-separated list of tag-value pairs,
        providing additional information about each feature."
        :param exprs: list of compiled regex expressions to use to find a gene_id or gene_name
        :return: str a gene_id or gene_name parsed out of attribute string
        """
        splt = attribute.split(';')
        gene = False
        expr_idx = 0
        while not gene:
            gene_matches = list(filter(exprs[expr_idx].match, splt))
            if gene_matches:
                gene = gene_matches[0].replace(exprs[expr_idx].pattern, '').strip(' "')

            expr_idx += 1

        return gene

    def get_data(self):
        """
        Load a .gtf or .gff file, extract exon regions, and subset to a pandas.DataFrame that looks like this:

        +----------------+-----------------+-----------------+--------------------+
        |     chr        |       start     |       end       |         gene       |
        +================+=================+=================+====================+
        |     chrI       |      11873      |      12227      |     MIR1302-11     |
        +----------------+-----------------+-----------------+--------------------+
        |     chr6       |      17232      |      17368      |     LINC00266-3    |
        +----------------+-----------------+-----------------+--------------------+

        :return: pandas.DataFrame for exon annotated regions with fields 'chr', 'start', 'end', 'gene'
        """
        # load file only if there are at least 9 columns.
        try:
            df = read_table(self.filename
                            , sep='\t'
                            , header=None
                            , usecols=list(range(9)))
        except ValueError as e:
            raise ValueError('file {0} must have the 9 mandatory .gtf/.gff columns.'
                             'Read more at https://useast.ensembl.org/info/website/upload/gff.html')

        cols = ['chr', 'source', 'feature', 'start',
                'end', 'score', 'strand', 'frame', 'attribute']
        df.columns = cols

        # subset annotation to just exons.
        df = df[df.feature.apply(lambda x: x.lower()) == 'exon']

        # parse out gene identifiers from attribute strings.
        find_me = [re.compile('gene_id'), re.compile('gene_name')]
        df['gene'] = df.attribute.apply(lambda x: self._attribute_to_gene(x, exprs=find_me))

        # subset to the data we'll actually need, turning data into a .bed file format.
        return df[['chr', 'start', 'end', 'gene']].drop_duplicates()


# if __name__ == '__main__':
#     gtf_file = os.path.join(os.getenv('HOME'), 'nu', 'jiping_research', 'data', 'rna_seq', 'genes.gtf')
#     loader = GeneAnnotationLoader(gtf_file)
#     df = loader.get_data()
#     print(df.head(20))