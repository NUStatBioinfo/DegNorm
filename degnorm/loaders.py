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

    def get_files(self, to_load):
        """
        Obtain list of str realpaths to files to load.

        :param to_load: str or list; if str, can be path to directory containing `filetype` files or
        the realpath to a `filetype` file. If list, a list of str realpaths to files ending in `filetype` extension.
        :return:
        """
        if isinstance(to_load, str):

            if os.path.isdir(to_load):
                dir_files = os.listdir(to_load)
                for f in dir_files:
                    if f.endswith(self.filetype):
                        self.files.append(os.path.join(to_load, f))

            else:
                if to_load.endswith(self.filetype):
                    self.files = [to_load]
                else:
                    raise ValueError('to_load file does not end with {0}'.format(to_load))

        elif isinstance(to_load, list):
            for f in dir_files:
                if f.endswith(self.filetype):
                    self.files.append(os.path.join(to_load, f))

        else:
            raise ValueError('to_load type not understood')

    def get_data(self):
        raise NotImplementedError('get_data not yet implemented for {0}'.format(self.__class__.__name__))


class SamLoader(Loader):

    def __init__(self, to_load):
        """
        .sam file loader

        :param to_load: str or list; if str, can be path to directory containing .sam files or
        the realpath to a .sam file. If list, a list of str .sam filenames.
        """
        Loader.__init__(self, '.sam')
        self.get_files(to_load)

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

        for fname in self.files:

            with open(fname, 'r') as f:
                lines = f.readlines()

            # identify line number where header lines stop.
            header_idx = 0
            is_header_line = lines[header_idx][0] == '@'
            while is_header_line:
                header_idx += 1
                is_header_line = lines[header_idx][0] == '@'

            # extract relevant columns from .sam file lines.
            reads = list()
            for x in lines[header_idx:]:
                splt = x.split('\t')
                reads.append([splt[i] for i in [0, 2, 3, 4, 5, 7, 8]])

            # transform .sam lines into pandas.DataFrame
            colnames = ['qname', 'chr', 'pos', 'mapq', 'cigar', 'pnext', 'tlen']
            df = DataFrame(reads
                           , columns=colnames)

            # typecasting: change int fields from str.
            int_cols = ['pos', 'pnext', 'tlen']
            for col in int_cols:
                df[col] = df[col].astype('int')

            float_cols = ['mapq']
            for col in float_cols:
                df[col] = df[col].astype('float')

            # preprocessing: un-specify paired reads together so that a pair of alignments
            # share the same un-paired QNAME
            df['qname_unpaired'] = df.qname.apply(lambda x: '.'.join(x.split('.')[:-1]))

            # preprocessing: identify max position of aligned reads. For this, need to identify
            # length of each alignment from cigar string and add that to read start position.
            df['end_pos'] = df['pos'] + df.cigar.apply(
                lambda x: np.sum([int(seg_len) for seg_len in re.findall(r'(\d+)', x)]))

            df_dict[fname] = df

        return df_dict


class BamLoader(Loader):

    def __init__(self, to_load):
        """
        .bam file loader

        :param to_load: str or list; if str, can be path to directory containing .bam files or
        the realpath to a .bam file. If list, a list of str .bam filenames.
        """
        Loader.__init__(self, '.bam')
        self.get_files(to_load)