import logging
import os
import shutil
import gzip

import numpy as np
import pandas as pd
import yaml


class getFileHandle(object):
    def __init__(self, filename, mode='rt'):
        self.filename = filename
        self.mode = mode

    def __enter__(self):
        if self.get_file_format(self.filename) in ["csv", 'plain-text', "yaml"]:
            self.handle = open(self.filename, self.mode)
        elif self.get_file_format(self.filename) == "gzip":
            self.handle = gzip.open(self.filename, self.mode)
        elif self.get_file_format(self.filename) == "h5":
            self.handle = pd.HDFStore(self.filename, self.mode)
        return self.handle

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handle.close()

    def get_file_format(self, filepath):
        if filepath.endswith('.tmp'):
            filepath = filepath[:-4]

        _, ext = os.path.splitext(filepath)

        if ext == ".csv":
            return "csv"
        elif ext == ".gz":
            return "gzip"
        elif ext == ".h5" or ext == ".hdf5":
            return "h5"
        elif ext == ".yaml":
            return "yaml"
        else:
            logging.getLogger("single_cell.helpers").warn(
                "Couldn't detect output format for {}. extension {}".format(filepath, ext)
            )
            return "plain-text"


class CsvParseError(Exception):
    pass

class CsvInputError(Exception):
    pass

class CsvWriterError(Exception):
    pass


def get_dtypes_from_df(df, na_rep='NA'):
    pandas_to_std_types = {
        "bool": "bool",
        "int64": "int",
        "float64": "float",
        "object": "str",
        "category": "str",
    }

    typeinfo = {}
    for column, dtype in df.dtypes.items():
        if column in ['chr', 'chrom', 'chromosome']:
            typeinfo[column] = 'str'
        else:
            if df.empty:
                typeinfo[column] = na_rep
            else:
                typeinfo[column] = pandas_to_std_types[str(dtype)]

    return typeinfo


class CsvInput(object):
    def __init__(self, filepath, na_rep='NA'):
        """
        csv file and all related metadata
        :param filepath: path to csv
        :type filepath: str
        :param na_rep: replace na with this
        :type na_rep: str
        """
        self.filepath = filepath
        self.compression = self.__get_compression_type_pandas()
        self.na_rep = na_rep

        if os.path.exists(self.yaml_file):
            metadata = self.__parse_metadata()
        else:
            metadata = self.generate_metadata()

        self.header, self.sep, self.dtypes, self.columns = metadata

    @property
    def yaml_file(self):
        return self.filepath + '.yaml'

    def __detect_sep_from_header(self, header):
        """
        detect whether file is tab or comma separated from header
        :param header: header line
        :type header: str
        :return: separator
        :rtype: str
        """
        if '\t' in header and ',' in header:
            raise CsvParseError("Unable to detect separator from {}".format(header))

        if '\t' not in header and ',' not in header:
            raise CsvParseError("Unable to detect separator from {}".format(header))

        if '\t' in header:
            return '\t'
        elif ',' in header:
            return ','

    def __detect_sep_from_file(self):
        with getFileHandle(self.filepath) as reader:
            header = reader.readline().strip()
            return self.__detect_sep_from_header(header)

    def __get_compression_type_pandas(self):
        filepath = self.filepath
        if filepath.endswith('.tmp'):
            filepath = filepath[:-4]

        _, ext = os.path.splitext(filepath)

        if ext == ".csv":
            return None
        elif ext == ".gz":
            return "gzip"
        elif ext == ".h5" or ext == ".hdf5":
            raise CsvInputError("HDF is not supported")
        else:
            logging.getLogger("single_cell.utils.csv").warn(
                "Couldn't detect output format for {}. extension {}".format(filepath, ext)
            )
            return None

    def __is_empty(self):
        with open(self.filepath) as f:
            first_line = f.readline()
            if len(first_line) == 0:
                logging.getLogger("single_cell.utils.csv").warn("empty csv file")
                return True

    def __generate_dtypes(self, columns=None, sep=','):
        data = pd.read_csv(
            self.filepath, compression=self.compression, chunksize=10 ** 6,
            sep=sep
        )
        data = next(data)

        if columns:
            data.columns = columns

        typeinfo = get_dtypes_from_df(data)
        return typeinfo

    def __parse_metadata(self):
        with getFileHandle(self.filepath + '.yaml') as yamlfile:
            yamldata = yaml.safe_load(yamlfile)

        header = yamldata['header']
        sep = yamldata.get('sep', ',')

        dtypes = {}
        columns = []
        for coldata in yamldata['columns']:
            colname = coldata['name']

            dtypes[colname] = coldata['dtype']

            columns.append(colname)

        return header, sep, dtypes, columns

    def generate_metadata(self):
        with getFileHandle(self.filepath) as inputfile:
            header = inputfile.readline().strip()
            sep = self.__detect_sep_from_header(header)
            columns = header.split(sep)
            header = True
            dtypes = self.__generate_dtypes(sep=sep)
            return header, sep, dtypes, columns

    def __verify_data(self, df, usecols):
        cols = self.columns
        if usecols is not None:
            cols = []
            for col in usecols:
                if col not in self.columns:
                    raise ValueError(f'requested column {col} not in columns')
            for col in self.columns:
                if col in usecols:
                    cols.append(col)
        if not list(df.columns.values) == cols:
            raise CsvParseError("metadata mismatch in {}".format(self.filepath))

    def read_csv(self, chunksize=None, dtypes_override=None, usecols=None):
        def return_gen(df_iterator):
            for df in df_iterator:
                self.__verify_data(df, usecols)
                yield df

        dtypes = {k: v for k, v in self.dtypes.items() if v != "NA"}

        if dtypes_override is not None:
            for name, dtype in dtypes_override.items():
                if name in dtypes:
                    dtypes[name] = dtype

        if self.header:
            header = 0
            names = None

        else:
            header = None
            names = self.columns

        try:
            data = pd.read_csv(
                self.filepath, compression=self.compression, chunksize=chunksize,
                sep=self.sep, header=header, names=names, usecols=usecols, dtype='str')
        except pd.errors.EmptyDataError:
            data = pd.DataFrame(columns=self.columns)

        for column_name, dtype in dtypes.items():
            if usecols is not None and column_name not in usecols:
                continue
            try:
                if dtype == 'int':
                    data[column_name] = data[column_name].astype(float).astype('Int64')
                elif dtype == 'bool':
                    values = set(data[column_name].unique())
                    if len(values - {False, True, np.nan}) != 0:
                        raise ValueError(f'{column_name} is expected to be bool, has values {values}')
                else:
                    data[column_name] = data[column_name].astype(dtype)
            except TypeError:
                logging.exception(f'unable to convert {column_name} to {dtype}')
                raise

        if chunksize:
            return return_gen(data)
        else:
            self.__verify_data(data, usecols)
            return data


class CsvOutput(object):
    def __init__(
            self, filepath, header=True, sep=',', columns=None, dtypes=None,
            na_rep='NA'
    ):
        self.filepath = filepath
        self.header = header
        self.columns = columns
        self.sep = sep
        self.dtypes = dtypes if dtypes else {}
        self.na_rep = na_rep

    @property
    def yaml_file(self):
        return self.filepath + '.yaml'

    @property
    def header_line(self):
        return self.sep.join(self.columns) + '\n'

    def __get_compression_type_pandas(self):
        filepath = self.filepath
        if filepath.endswith('.tmp'):
            filepath = filepath[:-4]

        _, ext = os.path.splitext(filepath)

        if ext == ".csv":
            return None
        elif ext == ".gz":
            return "gzip"
        elif ext == ".h5" or ext == ".hdf5":
            raise CsvInputError("HDF is not supported")
        else:
            logging.getLogger("single_cell.utils.csv").warn(
                "Couldn't detect output format for {}. extension {}".format(filepath, ext)
            )
            return None

    def __write_yaml(self):

        yamldata = {'header': self.header, 'sep': self.sep, 'columns': []}

        for column in self.columns:
            data = {'name': column, 'dtype': self.dtypes[column]}
            yamldata['columns'].append(data)

        with getFileHandle(self.yaml_file, 'wt') as f:
            yaml.safe_dump(yamldata, f, default_flow_style=False)

    def write_df(self, df):
        if self.columns:
            if not self.columns == df.columns.values:
                raise CsvWriterError("Writer initialized with wrong col names")

        header = df.columns.values if self.header else False
        compression = self.__get_compression_type_pandas()
        df.to_csv(
            self.filepath, sep=self.sep, na_rep=self.na_rep,
            index=False, header=header, compression=compression
        )

        self.columns = list(df.columns.values)
        self.dtypes = get_dtypes_from_df(df)

        self.__write_yaml()

    def write_csv_data(self, reader, writer):
        reader_gzip = type(reader) == gzip.GzipFile
        writer_gzip = type(writer) == gzip.GzipFile

        if reader_gzip == writer_gzip:
            same_format = True
        else:
            same_format = False

        if same_format:
            shutil.copyfileobj(
                reader, writer, length=16 * 1024 * 1024
            )
        else:
            for line in reader:
                writer.write(line)


    def concatenate_files(self, infiles):
        header = self.header_line if self.header else None

        with getFileHandle(self.filepath, 'wt') as writer:
            if header:
                writer.write(header)
            for infile in infiles:
                with getFileHandle(infile) as reader:
                    self.write_csv_data(reader, writer)

        self.__write_yaml()

    def write_headerless_csv(self, infile):
        with getFileHandle(self.filepath, 'wt') as writer:
            with getFileHandle(infile) as reader:
                if not reader.readline() == self.header_line:
                    raise CsvWriterError("cannot write, wrong header")
                self.write_csv_data(reader, writer)

        self.__write_yaml()

    def write_csv_with_header(self, infile):
        with getFileHandle(self.filepath, 'wt') as writer:
            with getFileHandle(infile) as reader:
                writer.write(self.header_line)
                self.write_csv_data(reader, writer)

        self.__write_yaml()


def annotate_csv(infile, annotation_data, outfile, on="cell_id", write_header=True):
    csvinput = CsvInput(infile)

    metrics_df = csvinput.read_csv()

    merge_on = metrics_df[on]

    for cell in merge_on:
        col_data = annotation_data[cell]

        for column, value in col_data.items():
            metrics_df.loc[metrics_df[on] == cell, column] = value

    output = CsvOutput(outfile, sep=csvinput.sep, header=write_header)
    output.write_df(metrics_df)


def concatenate_csv(in_filenames, out_filename, key_column=None, write_header=True):
    if not isinstance(in_filenames, dict):
        in_filenames = dict(enumerate(in_filenames))

    data = []
    sep = None

    for key, in_filename in in_filenames.items():
        csvinput = CsvInput(in_filename)

        if not sep:
            sep = csvinput.sep
        assert sep == csvinput.sep

        df = csvinput.read_csv()

        if key_column is not None:
            df[key_column] = str(key)
        data.append(df)
    data = pd.concat(data, ignore_index=True)

    csvoutput = CsvOutput(out_filename, header=write_header, sep=sep)
    csvoutput.write_df(data)


def extrapolate_types_from_yaml_files(csv_files):
    precedence = ['str', 'float', 'int', 'bool', 'NA']

    csv_metadata = [CsvInput(csv_file) for csv_file in csv_files]

    header = set([val.header for val in csv_metadata])
    assert len(header) == 1, 'mismatched yaml files'
    header = list(header)[0]

    sep = set([val.sep for val in csv_metadata])
    assert len(sep) == 1, 'mismatched yaml files'
    sep = list(sep)[0]

    cols = set([tuple(val.columns) for val in csv_metadata])
    assert len(cols) == 1, 'mismatched yaml files'
    cols = list(cols)[0]

    dtypes = [val.dtypes for val in csv_metadata]

    merged_dtypes = {}
    for dtype_val in dtypes:
        for col, dtype in dtype_val.items():
            if col not in merged_dtypes:
                merged_dtypes[col] = dtype
            else:
                og_index = precedence.index(merged_dtypes[col])
                new_index = precedence.index(dtype)
                if new_index < og_index:
                    merged_dtypes[col] = dtype

    return header, sep, cols, merged_dtypes


def concatenate_csv_files_quick_lowmem(inputfiles, output, write_header=True):
    if isinstance(inputfiles, dict):
        inputfiles = inputfiles.values()

    header, sep, cols, dtypes = extrapolate_types_from_yaml_files(inputfiles)

    if header:
        raise CsvInputError("Attempting to concatenate files with header.")

    csvoutput = CsvOutput(
        output, header=write_header, sep=sep, columns=cols, dtypes=dtypes
    )
    csvoutput.concatenate_files(inputfiles)


def prep_csv_files(filepath, outputfile, header=False):
    """
    generate header less csv files
    :param filepath:
    :type filepath:
    :param outputfile:
    :type outputfile:
    """
    csvinput = CsvInput(filepath)

    if csvinput.header is None:
        raise CsvInputError("cannot prep csv file with no header")

    csvoutput = CsvOutput(
        outputfile, header=False, columns=csvinput.columns,
        sep=csvinput.sep, dtypes=csvinput.dtypes
    )

    if header:
        csvoutput.write_csv_with_header(filepath)
    else:
        csvoutput.write_headerless_csv(filepath)


def finalize_csv(infile, outfile):
    csvinput = CsvInput(infile)

    if csvinput.header:
        raise CsvInputError("cannot finalize file with header")

    csvoutput = CsvOutput(
        outfile, header=True, columns=csvinput.columns,
        sep=csvinput.sep, dtypes=csvinput.dtypes
    )
    csvoutput.write_csv_with_header(infile)


def merge_csv(in_filenames, out_filename, how, on, nan_val='NA', suffixes=None, write_header=True):
    if isinstance(in_filenames, dict):
        in_filenames = in_filenames.values()

    data = []
    sep = None
    for in_filename in in_filenames:
        csvinput = CsvInput(in_filename)
        indata = csvinput.read_csv()
        if not indata.empty:
            data.append(indata)

        if not sep:
            sep = csvinput.sep
        assert sep == csvinput.sep

    data = merge_frames(data, how, on, suffixes=suffixes)
    data = data.fillna(nan_val)

    csvoutput = CsvOutput(out_filename, header=write_header, sep=sep)
    csvoutput.write_df(data)


def merge_frames(frames, how, on, suffixes=None):
    '''
    annotates input_df using ref_df
    '''

    suff = ['', '']

    if ',' in on:
        on = on.split(',')

    if len(frames) == 1:
        return frames[0]
    else:
        left = frames[0]
        right = frames[1]

        cols_to_use = list(right.columns.difference(left.columns))
        cols_to_use += on
        cols_to_use = list(set(cols_to_use))

        if suffixes:
            suff = (suffixes[0], suffixes[1])

        merged_frame = pd.merge(left, right[cols_to_use],
                                how=how,
                                on=on,
                                suffixes=suff)
        for i, frame in enumerate(frames[2:]):

            if suffixes:
                suff = (suffixes[i + 2], suffixes[i + 2])

            merged_frame = pd.merge(merged_frame, frame,
                                    how=how,
                                    on=on,
                                    suffixes=suff)
        return merged_frame


def read_csv_and_yaml(infile, chunksize=None):
    return CsvInput(infile).read_csv(chunksize=chunksize)


def write_dataframe_to_csv_and_yaml(df, outfile, write_header=False, sep=','):
    csvoutput = CsvOutput(outfile, header=write_header, sep=sep)
    csvoutput.write_df(df)

def get_metadata(infile):
    csvinput = CsvInput(infile)
    return csvinput.header, csvinput.dtypes, csvinput.columns
