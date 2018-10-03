import abc
import pandas as pd

from enum import Enum
from StringIO import StringIO



class _FNLOOutputParserBase(object):

    __metaclass__ = abc.ABCMeta

    def __init__(self, filename):
        self._fname = filename
        self._data = {}

    # -- helper methods

    @staticmethod
    def _parse_table_with_pandas(columns, table_lines, add_columns_dict=None):
        # construct pandas table from preprocessed table lines
        _table_string = ' '.join(columns) + '\n' + '\n'.join(table_lines)
        _table_sio = StringIO(_table_string)
        _df = pd.read_table(_table_sio, delim_whitespace=True)

        # add new columns from dict
        if add_columns_dict is not None:
            for _col_name, _col_value in add_columns_dict.iteritems():
                _df[_col_name] = _col_value
        return _df

    # -- main parser method (implement in derived classes)

    @abc.abstractmethod
    def _parse(self):
        raise NotImplementedError

    # -- public API

    def get_var(self, var_name):
        if 'table' not in self._data:
            self._parse()
        return self._data.get(var_name, None)

    def get_table(self):
        if 'table' not in self._data:
            self._parse()
        return self._data['table']



class FNLOYodaOutputParser(_FNLOOutputParserBase):

    TOKEN_BEGIN_SECTION_UNC = "fnlo-tk-yodaout: Evaluating uncertainties"

    class ParserMode(Enum):
        general = 0
        parsing_unc_section = 1
        parsing_unc_name = 2
        parsing_header = 3
        parsing_table = 4

    @staticmethod
    def _read_header(line):
        """parse the header line and return column names"""
        line = line.replace('#', ' ').strip()
        return [_colname.strip() for _colname in line.split()]

    def _parse(self):
        self._data = dict(subtables=[])
        with open(self._fname, 'r') as _f:
            _mode = self.ParserMode.general
            for line in _f:
                if _mode == self.ParserMode.general:
                    # check for beginning of table section
                    if self.TOKEN_BEGIN_SECTION_UNC in line:
                        _new_columns_dict = {}
                        _mode = self.ParserMode.parsing_unc_section
                        continue

                elif _mode == self.ParserMode.parsing_unc_section:
                    # check for beginning of table section
                    if '#=' in line:
                        _mode = self.ParserMode.parsing_unc_name
                        continue
                    elif '#-' in line:
                        _mode = self.ParserMode.parsing_header
                        continue

                elif _mode == self.ParserMode.parsing_unc_name:
                     self._data['uncertainty_label'] = line.replace('#', ' ').strip()
                     _mode = self.ParserMode.parsing_unc_section
                     continue

                elif _mode == self.ParserMode.parsing_header:
                    # read columns and switch mode
                    _columns = self._read_header(line)
                    _table_lines = []
                    _mode = self.ParserMode.parsing_table
                    continue

                elif _mode == self.ParserMode.parsing_table:
                    # skip empty lines
                    if '#' in line:
                        line = line.split('#', 1)[0].strip()

                    if not line:
                        continue

                    # accumulate table lines
                    _table_lines.append(line)

        # construct pandas table
        self._data['table'] = self._parse_table_with_pandas(_columns, _table_lines, add_columns_dict=None)



class FNLOCppReadOutputParser(_FNLOOutputParserBase):

    TOKEN_BEGIN_SECTION_XS = "[fnlo-tk-cppread] Calculate my cross sections"
    TOKEN_SCALE_FACTORS = "The scale factors xmur, xmuf chosen here are:"
    TOKEN_FIXED_SCALES = "The fixed scales mur, muf chosen here are:"
    TOKEN_PDF_MEMBER = "The PDF member chosen here is:"


    class ParserMode(Enum):
        general = 0
        parsing_xs_section = 1
        parsing_header = 2
        parsing_table = 3


    @staticmethod
    def _read_header(line):
        """parse the header line and return column names"""
        line = line.replace('#', ' ')
        _tokens = line.split()
        _column_names = []
        _flag_reading_binlimits = False
        _flag_reading_meanvalue = False
        for _t in _tokens:
            # -- special tokens
            if _t == '[':
                # '[' marks beginning of bin limit columns
                assert not _flag_reading_binlimits
                _flag_reading_binlimits = True
                continue
            elif _t == ']':
                # ']' marks end of bin limit columns
                assert _flag_reading_binlimits
                _flag_reading_binlimits = False
                continue
            elif _t == '<':
                # '<' marks beginning of mean value column
                assert not _flag_reading_meanvalue
                _flag_reading_meanvalue = True
                continue
            elif _t == '>':
                # '>' marks end of mean value column
                assert _flag_reading_meanvalue
                _flag_reading_meanvalue = False
                continue

            # -- construct column names
            if _flag_reading_binlimits:
                _column_names.extend([_t+'_lo', _t+'_hi'])
            elif _flag_reading_meanvalue:
                _column_names.append('<{}>'.format(_t))
            else:
                _column_names.append(_t)

        return _column_names

    def _parse(self):
        self._data = dict(subtables=[])
        with open(self._fname, 'r') as _f:
            _mode = self.ParserMode.general
            for line in _f:
                if _mode == self.ParserMode.general:
                    # check for beginning of table section
                    if self.TOKEN_BEGIN_SECTION_XS in line:
                        _new_columns_dict = {}
                        _mode = self.ParserMode.parsing_xs_section
                        continue

                elif _mode == self.ParserMode.parsing_xs_section:
                    # check for line specifying scale factors/fixed scales/PDF member
                    if self.TOKEN_SCALE_FACTORS in line:
                        _token_string_parameter = line.split(self.TOKEN_SCALE_FACTORS, 1)[1]
                        _cmur, _cmuf = map(float, [s.strip() for s in _token_string_parameter.split(',')])
                        _new_columns_dict['cmur'] = _cmur
                        _new_columns_dict['cmuf'] = _cmuf
                        continue
                    elif self.TOKEN_FIXED_SCALES in line:
                        _token_string_parameter = line.split(self.TOKEN_FIXED_SCALES, 1)[1]
                        _mur, _muf = map(float, [s.strip() for s in _token_string_parameter.split(',')])
                        _new_columns_dict['mur'] = _mur
                        _new_columns_dict['muf'] = _muf
                        continue
                    elif self.TOKEN_PDF_MEMBER in line:
                        _token_string_parameter = line.split(self.TOKEN_PDF_MEMBER, 1)[1]
                        _pdf_member = int(_token_string_parameter.strip())
                        _new_columns_dict['pdf_member'] = _pdf_member
                        continue
                    elif '#-' in line:
                        _mode = self.ParserMode.parsing_header
                        continue

                elif _mode == self.ParserMode.parsing_header:
                    # read columns and switch mode
                    _columns = self._read_header(line)
                    _table_lines = []
                    _mode = self.ParserMode.parsing_table
                    continue

                elif _mode == self.ParserMode.parsing_table:
                    # end of table
                    if '#=' in line:
                        self._data['subtables'].append(self._parse_table_with_pandas(_columns, _table_lines, add_columns_dict=_new_columns_dict))

                        _mode = self.ParserMode.parsing_xs_section
                        continue

                    # skip empty lines
                    if '#' in line:
                        line = line.split('#', 1)[0].strip()

                    if not line:
                        continue

                    # accumulate table lines
                    _table_lines.append(line)

        # construct final subtable
        if _mode == self.ParserMode.parsing_table:
            self._data['subtables'].append(self._parse_table_with_pandas(_columns, _table_lines, add_columns_dict=_new_columns_dict))

        # construct concatenated table
        self._data['table'] = pd.concat(self._data['subtables'], ignore_index=True)

        # clean up subtables
        del self._data['subtables']

    def get_table(self):
        if 'table' not in self._data:
            self._parse()
        return self._data['table']



if __name__ == "__main__":
    import argparse

    AVAILABLE_FORMATS = ('yoda', 'cppread')

    ap = argparse.ArgumentParser()
    ap.add_argument("FNLO_LOG_FILE", help="Log file output by fnlo-tk-cppread.")
    ap.add_argument("-f", '--format', help="Format (available: {}).".format(", ".join(AVAILABLE_FORMATS)), default='cppread')

    args = ap.parse_args()


    # initialize parser
    _format_spec = args.format.strip()
    if _format_spec == 'cppread':
        p = FNLOCppReadOutputParser(args.FNLO_LOG_FILE)
    elif _format_spec == 'yoda':
        p = FNLOYodaOutputParser(args.FNLO_LOG_FILE)
    else:
        raise ValueError("Unknown format specification '{}'. Available formats: {}".format(_format_spec, AVAILABLE_FORMATS))

    # parse logfile and retrieve data table
    t = p.get_table()

    # print whole table
    print t

    if _format_spec == 'yoda':
        print 'uncertainty_label = "{}"'.format(p.get_var("uncertainty_label"))
