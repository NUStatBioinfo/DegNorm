import os
import re
import argparse
import numpy as np
from pandas import DataFrame


class Loader():
    def __init__(self, filetype):
        """
        Generic file loader class.

        :param filetype: str file extension, only load files that .endswith(filetype)
        """
        self.filetype = filetype
        self.files = list()

    def get_file(self, to_load):
        """
        Obtain list of str realpaths to files to load.

        :param to_load: str the realpath to a `filetype` file.
        """
        if isinstance(to_load, str):
            if os.path.exists(to_load):
                if to_load.endswith(self.filetype):
                    self.filename = to_load
                else:
                    raise ValueError('to_load file {0} does not end with {0}'.format(to_load, self.filetype))
            else:
                raise IOError('to_load file {0} not found'.format(to_load))

        else:
            raise ValueError('to_load data type not understood')

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

        df_dict['data'] = df

        return df_dict


class BamLoader(Loader):

    def __init__(self, to_load):
        """
        .bam file loader

        :param to_load: str the realpath to a .bam file.
        """
        Loader.__init__(self, '.bam')
        self.get_file(to_load)


class GtfLoader(Loader):

    def __init__(self, to_load):
        """
        .gtf file loader

        :param to_load: str the realpath to a .gtf file.
        """
        Loader.__init__(self, '.gtf')
        self.get_file(to_load)

    def get_data(self):
