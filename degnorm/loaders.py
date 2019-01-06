import pysam
from degnorm.utils import *
from pandas import read_table


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
                raise FileNotFoundError('file {0} not found'.format(to_load))

        else:
            raise ValueError('{0} data type not understood'.format(to_load))

    def get_data(self):
        raise NotImplementedError('get_data not yet implemented for {0}'.format(self.__class__.__name__))


class BamLoader(Loader):

    def __init__(self, to_load, index_file):
        """
        .bam file loader. Uses pysam utility.

        :param to_load: str the realpath to a .bam file.
        """

        Loader.__init__(self, '.bam')
        self.get_file(to_load)

        # check that supporting .bai index file exists.
        if not os.path.isfile(index_file):
            raise FileNotFoundError('{0} .bam index file not found')
        else:
            if not index_file.endswith('.bai'):
                raise ValueError('.bam index file does not have correct .bai file extension.')

        self.index_file = index_file

    def get_data(self):
        bamfile = pysam.AlignmentFile(self.filename
                                      , mode='rb'
                                      , index_filename=self.index_file)

        return bamfile


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
        return df[['chr', 'start', 'end', 'gene']].drop_duplicates().reset_index(drop=True)