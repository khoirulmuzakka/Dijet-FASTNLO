#!/usr/bin/env python3

import numpy as np
import re
import math
import sys
import os
import nnlojet_algo

import logging

_comment_prefix = '#'


class NNLOJETHistogram():

    def __init__(self, **kwargs):
        # nx = None checks for empty Histograms
        self._nx = kwargs.get('nx', None)
        self._filename = kwargs.get('filename', None)
        # default value for nx in case we have a file attached to the object
        if self.filename is not None and self.nx is None:
            self._nx = 3
        if self._nx is not None and self._nx < 0:
            raise ValueError('invaid number of x-columns:{}'.format(nx))
        self._neval  = None  # set in `__parse_comments`
        self._labels = None  # set in `__parse_comments`
        # instead of saving labels, add a list of comments: __header, __footer
        # later loop and `not startswith(_comment_prefix): _comment_prefix + line`
        self._xval = None  # set in `__load_data`
        self._yval = None  # set in `__load_data`
        self._yerr = None  # set in `__load_data`
        self._ioverflow = None  # set in `__load_data`

        if self.filename is not None:
            columns = kwargs.get('columns', None)
            rebin = kwargs.get('rebin', None)
            if rebin is not None and self.nx <= 0:
                raise ValueError('rebinning requires x-columns!')
            self._read_dat(columns=columns, rebin=rebin)

        # to keep track of weights
        self._files_wgt = None  # list of the files
        self._wgt = None
        if 'weights' in kwargs and kwargs.get('weights', False):
            self._load_wgt()
            if kwargs.get('recursive_weights', False):
                self._load_recursive_wgt()

    @property
    def nx(self):
        return self._nx


    @property
    def filename(self):
        return self._filename


    @staticmethod
    def filename_suff(filename, suff='.npz'):
        if filename is None:
            return None
        if filename.endswith('.dat'):
            return filename[:-4] + suff
        else:
            return filename + suff


    @property
    def filename_wgt(self):
        return NNLOJETHistogram.filename_suff(self.filename)


    # without the setter, we do not allow setting the value after init
    # @nx.setter
    # def nx(self, nx): self._nx = nx


    def _read_dat(self, **kwargs):

        columns = kwargs.get('columns', None)
        rebin = kwargs.get('rebin', None)

        # step 1: parse label line
        with open(self.filename, 'rt') as histfile:
            lines = histfile.readlines()

        # fh = open(self.filename, 'rt')
        for line in lines:
            if re.match(r'^\s*' + _comment_prefix + 'labels', line, re.IGNORECASE):
                self._labels = [ re.sub(r'\[[0-9]+\]', '', lab) for lab in line.split()[1:] ]  # remove "#labels" part
                # break
            if re.match(r'^\s*' + _comment_prefix + 'nx', line, re.IGNORECASE):
                self._nx = int(line.split()[1]) # remove "#nx" part
                # print("using nx = {}".format(self._nx))
        # fh.close()

        # warning if no labels line could be found
        if self._labels is None:
            logging.warning("WARNING: couldn't find a labels line {}".format(self.filename))
            return

        # default init
        self._ioverflow = None
        nrows = 0
        ncols = len(self._labels)  # raw column size

        # build index list (only y-columns; no x-val)
        if columns is not None:
            idx_cols = []
            for col in columns:
                if re.match(r'.*_Err\b', col, re.IGNORECASE):
                    raise ValueError('columns must not contain *_Err column names!')
                idx_cols.append(self._labels.index(col))
            # overwrite labels
            if self.nx > 0:
                self._labels = self._labels[:self.nx]
            else:
                self._labels = []
            for col in columns:
                self._labels.extend( [col, col+'_Err'] )
        else:
            idx_cols = range(self.nx, ncols, 2)

        # step 2: read in the data
        comment_pattern = re.compile(r'^\s*' + _comment_prefix) 

        # dynamic arrays: python lists so much faster than numpy arrays
        data_xval = [] if self.nx > 0 else None
        data_yval = []
        data_yerr = []

        # fh = open(self.filename, 'rt')
        for line in lines:
            #----- parse comments
            if re.match(comment_pattern, line):
                # neval
                if re.match(r'^\s*' + _comment_prefix + 'neval', line, re.IGNORECASE):
                    self._neval = int(line.split()[1])
                # overflow
                if re.match(r'^\s*' + _comment_prefix + 'overflow', line, re.IGNORECASE):
                    data = line.split()  # make a list
                    if len(data) != ncols:
                        raise ValueError('ncols mismatch in overflow: {} != {} (file: {})'.format(len(data), ncols, self.filename))
                    if self.nx > 0:
                        data_xval.append( [0.] * self.nx )
                    data_yval.append( [ float(data[idx]) for idx in idx_cols ] )
                    data_yerr.append( [ float(data[idx+1]) for idx in idx_cols ] )
                    self._ioverflow = nrows
                    nrows += 1
            #----- read data
            else:
                data = line.split()
                if len(data) != ncols:
                    raise ValueError('ncols mismatch in data: {} != {} (file: {})'.format(len(data), ncols, self.filename))
                if self.nx > 0:
                    data_xval.append( [ float(data[idx]) for idx in range(self.nx) ] )
                data_yval.append( [ float(data[idx]) for idx in idx_cols ] )
                data_yerr.append( [ float(data[idx+1]) for idx in idx_cols ] )
                nrows += 1
        # fh.close()

        # set neval to 1 if there was none
        if self._neval is None:
            raise ValueError("couldn't find the neval information")
            print("couldn't find the neval info: {}".format(self.filename), file=sys.stderr)
            self._neval = 1
        
        # step 3: rebin 
        if rebin is not None:
            if self.nx != 3:
                raise ValueError('rebinning can only handle nx=3 for now')
            rebin_xval = []
            rebin_yval = []
            rebin_yerr = []
            rebin_row_val, rebin_row_err = None, None
            xlow, xupp = None, None
            i_rebin = 0
            i_overflow = None
            for i_row in range(nrows):
                # overflow is added immediately
                if i_row == self._ioverflow:
                    rebin_xval.append(data_xval[i_row])
                    rebin_yval.append(data_yval[i_row])
                    rebin_yerr.append(data_yerr[i_row])
                    i_overflow = i_rebin
                    #> do *not* increment! otherwise, we lose a bin
                    #i_rebin += 1
                    continue
                # # debug
                # print("data_xval[{}][0] = ".format(i_row), type(data_xval[i_row][0]), data_xval[i_row][0])
                # print("rebin[i_rebin] = ".format(i_rebin), type(rebin[i_rebin]), rebin[i_rebin])
                # find the first matching bin
                if data_xval[i_row][0] < rebin[i_rebin] and data_xval[i_row][2] <= rebin[i_rebin]:
                    continue
                # skip trailing bins
                if i_rebin == len(rebin)-1:
                    continue
                # init
                if xlow is None:
                    rebin_row_val = [0.] * len(data_yval[i_row])
                    rebin_row_err = [0.] * len(data_yval[i_row])
                    xlow = data_xval[i_row][0]
                # accumulate data
                if data_xval[i_row][0] >= rebin[i_rebin] and data_xval[i_row][2] <= rebin[i_rebin+1]:
                    delta = data_xval[i_row][2] - data_xval[i_row][0]
                    for i_col in range(len(rebin_row_val)):
                        rebin_row_val[i_col] +=  self._neval * delta*data_yval[i_row][i_col]
                        rebin_row_err[i_col] += (self._neval * delta*data_yerr[i_row][i_col])**2 + self._neval * (delta*data_yval[i_row][i_col])**2
                # finalise
                if i_row == nrows-1 or data_xval[i_row+1][0] >= rebin[i_rebin+1]:
                    xupp = data_xval[i_row][2]
                    delta = xupp - xlow
                    # print("delta = ", xlow, xupp, delta)
                    # print("rebin = ", rebin[i_rebin], rebin[i_rebin+1], (rebin[i_rebin+1]-rebin[i_rebin]))
                    rebin_xval.append( [xlow, (xlow+xupp)/2., xupp] )
                    for i_col in range(len(rebin_row_val)):
                        rebin_row_val[i_col] /= self._neval
                        rebin_row_err[i_col] = math.sqrt( rebin_row_err[i_col] - self._neval * rebin_row_val[i_col]**2 ) / self._neval
                        # differential:
                        rebin_row_val[i_col] /= delta
                        rebin_row_err[i_col] /= delta
                    rebin_yval.append(rebin_row_val)
                    rebin_yerr.append(rebin_row_err)
                    i_rebin += 1
                    # reset accumulators
                    rebin_row_val, rebin_row_err = None, None
                    xlow, xupp = None, None
                    continue

            self._ioverflow = i_overflow
            data_xval = rebin_xval
            data_yval = rebin_yval
            data_yerr = rebin_yerr

        # step 4: store to member vars as numpy arrays
        if data_xval is not None:
            self._xval = np.array(data_xval, dtype='float') 
        self._yval = np.array(data_yval, dtype='float') 
        self._yerr = np.array(data_yerr, dtype='float') 

        # warning if all entries are zero
        if (not np.any(self._yval)) and (not np.any(self._yerr)):
            logging.warning("WARNING: all entries zero for {}".format(self.filename))

        # warning if nan or inf
        if np.any(np.logical_not(np.isfinite(self._yval))) or np.any(np.logical_not(np.isfinite(self._yerr))):
            logging.warning("WARNING: NaN or Inf encountered in {}".format(self.filename))


    def _load_wgt(self):
        if self.filename_wgt is None:
            raise ValueError("tried to load weights although no filename attached")
        if os.path.isfile(self.filename_wgt):
            npzfile = np.load(self.filename_wgt)
            self._files_wgt = npzfile['files']
            self._wgt = npzfile['weights']
        else:  # initialise weights to one
            self._files_wgt = np.array([self.filename])
            self._wgt = np.expand_dims(np.ones(self._yval.shape, dtype='float'), axis=2)


    def _load_recursive_wgt(self):
        while True:
            qnew = False
            for (i_file, datfile) in enumerate(self._files_wgt):
                wgtfile = NNLOJETHistogram.filename_suff(datfile)
                if os.path.isfile(wgtfile):
                    # save weight parent weight factors
                    fac_wgt = np.expand_dims(self._wgt[:, :, i_file], axis=2)
                    # delete `i_file` entries
                    self._wgt = np.delete(self._wgt, i_file, axis=2)
                    self._files_wgt = np.delete(self._files_wgt, i_file)
                    # read in the sub-weights
                    npzfile = np.load(wgtfile)
                    sub_files_wgt = npzfile['files']
                    sub_wgt = npzfile['weights']
                    # multipy the parent weights
                    sub_wgt = sub_wgt * fac_wgt
                    # append the sub-weights
                    self._wgt = np.append(self._wgt, sub_wgt, axis=2)
                    self._files_wgt = np.append(self._files_wgt, sub_files_wgt)
                    qnew = True
                    break
            if not qnew: break


    def __str__(self):
        # # set output format for floats
        # np.set_printoptions(formatter={'float': '{:E}'.format})
        # construct string line by line
        lines = []
        # first the comments
        lines.append('#labels: ' + ' '.join( [ lab + '[{}]'.format(idx+1) for (idx, lab) in enumerate(self._labels) ] ))
        lines.append('#neval: {}'.format(self._neval))
        # overflow?
        for irow, yarr in enumerate(self._yval):
            line = ''
            if self._xval is not None:
                if self._ioverflow is not None and irow == self._ioverflow:
                    if self.nx == 1:
                        line += '#overflow: '
                    else:
                        line += '#overflow:lower center upper '
                else:
                    for x in self._xval[irow]:
                        line += '{:.11E} '.format(x)
            for icol, y in enumerate(yarr):
                line += '{:.11E} {:.11E} '.format(y, self._yerr[irow, icol])
            lines.append(line)
        lines.append('#nx: {}'.format(self._nx))
        return '\n'.join(lines)


    def to_APPLfast_weight(self):
        # construct string line by line
        lines = []

        # x-bins
        line = '#nx={} '.format(self.nx)
        if self.nx > 0:
            for ix, xrow in enumerate(self._xval):
                if ix == self._ioverflow:  # skipping overflow
                    continue
                line += '[{:.11E},{:.11E}] '.format(xrow[0], xrow[-1])
        lines.append(line)

        # files & weights
        nfiles = len(self._files_wgt)
        for i_file in range(nfiles):
            line = self._files_wgt[i_file]
            for iwgt, wgt in enumerate(self._wgt[:, 0, i_file]):  # only first column
                if iwgt == self._ioverflow:  # skipping overflow
                    continue
                line += ' {:.11E}'.format(wgt)
            lines.append(line)

        return '\n'.join(lines)


    def write_to_file(self, filename):

        # if self.filename is not None and self.filename != filename:
        #     print("creating duplicate: {} > {}?".format(self.filename, filename))

        # write to the file
        with open(filename, 'wt') as fout:
            print(self, file=fout)
            # override internal filename
            self._filename = filename

        if self._wgt is not None:
            np.savez(self.filename_wgt, files=self._files_wgt, weights=self._wgt)


    def __add__(self, other):
        # check for the correct type
        if not isinstance(other, NNLOJETHistogram):
            raise ValueError("not adding a histogram object")

        # catch the two cases where at least one is unset
        if self.nx is not None and other.nx is None:
            return self
        if self.nx is None:  # also covers: both None
            return other

        # check for consistency
        if self.nx != other.nx:
            raise ValueError("nx mismatch")
        if self._labels != other._labels:
            raise ValueError("labels mismatch")
        if not np.array_equal(self._xval, other._xval):
            raise ValueError("xval mismatch")
        if self._ioverflow != other._ioverflow:
            raise ValueError("ioverflow mismatch")
        if self._yval.shape != other._yval.shape:
            raise ValueError("yval shape mismatch")

        result = NNLOJETHistogram(nx=self.nx)

        result._neval = self._neval + other._neval
        result._labels = self._labels
        result._ioverflow = self._ioverflow

        result._xval = self._xval
        result._yval = np.copy(self._yval)
        result._yerr = np.square(self._yerr)

        (n_rows, n_cols) = self._yval.shape
        for i_row in range(n_rows):
            for i_col in range(n_cols):
                result._yval[i_row, i_col] += other._yval[i_row, i_col]
                result._yerr[i_row, i_col] += other._yerr[i_row, i_col]**2
        result._yerr = np.sqrt(result._yerr)

        # if both have weights: determine combined weight table
        if self._wgt is not None and other._wgt is not None:
            result._files_wgt = np.append(self._files_wgt, other._files_wgt)
            result._wgt = np.append(self._wgt, other._wgt, axis=2)

        return result


    def overwrite(self, other):
        # check for the correct type
        if not isinstance(other, NNLOJETHistogram):
            raise ValueError("not adding a histogram object")

        # catch the two cases where at least one is unset
        if self.nx is not None and other.nx is None:
            return self
        if self.nx is None:  # also covers: both None
            return other

        # check for consistency
        if self.nx != other.nx:
            raise ValueError("nx mismatch")
        if self._labels != other._labels:
            raise ValueError("labels mismatch")
        if not np.array_equal(self._xval, other._xval):
            raise ValueError("xval mismatch")
        if self._ioverflow != other._ioverflow:
            raise ValueError("ioverflow mismatch")
        if self._yval.shape != other._yval.shape:
            raise ValueError("yval shape mismatch")

        result = NNLOJETHistogram(nx=self.nx)

        result._neval = self._neval + other._neval
        result._labels = self._labels
        result._ioverflow = self._ioverflow

        result._xval = self._xval
        result._yval = self._yval
        result._yerr = self._yerr

        (n_rows, n_cols) = self._yval.shape
        for i_row in range(n_rows):
            for i_col in range(n_cols):
                if other._yval[i_row, i_col] == 0. and other._yerr[i_row, i_col] == 0.:
                    continue
                result._yval[i_row, i_col] = other._yval[i_row, i_col]
                result._yerr[i_row, i_col] = other._yerr[i_row, i_col]

        # if both have weights: determine combined weight table
        if self._wgt is not None and other._wgt is not None:
            result._files_wgt = np.append(self._files_wgt, other._files_wgt)

            result._wgt = np.copy(self._wgt)
            for i_row in range(n_rows):
                for i_col in range(n_cols):
                    if other._yval[i_row, i_col] == 0. and other._yerr[i_row, i_col] == 0.:
                        continue
                    result._wgt[i_row, i_col, :] = 0.

            result._wgt = np.append(result._wgt, other._wgt, axis=2)

        return result


    def multiply_weights(self, wgts):
        if len(wgts)+1 != self._yval.shape[0]:
            raise ValueError("length of weight array does not match")

        fac = np.array(wgts, dtype='float')
        # add an empty overflow bin
        fac = np.insert(fac, self._ioverflow, 0.)
        # reshape for numpy matrix multiplication
        fac = np.expand_dims(fac, axis=1)
        self._yval = self._yval * fac
        self._yerr = self._yerr * fac


class NNLOJETContainer():


    def __init__(self, **kwargs):
        self._nx = None
        self._filename = None
        self._neval  = None
        self._labels = None
        self._xval = None
        self._yval = None
        self._yerr = None
        self._ioverflow = None

        # the mask always has the same shape as yval & yerr
        # 0:   unmasked
        # 1:   trimmed
        # 2:   invalid (nan / unused buffer entry)
        # < 0: `-(n+1)`, where `n` is where the run was merged into
        #      (useful to reconstruct merging history of k-optimisation)
        self._mask = None

        # it's much more efficient pre-allocating a sufficiently
        # big array, instead of dynamically growingit
        self._buffer_size = kwargs.get('size', None)
        self._size = 0

        # flag to switch on/off weight tables in generated histograms
        self._qwgt = kwargs.get('weights', False)


    @property
    def nx(self): return self._nx


    def append(self, hist):

        # check for the correct type
        if not isinstance(hist, NNLOJETHistogram):
            raise ValueError("not appending a histogram object")

        if self._qwgt and hist.filename is None:
            raise ValueError("only flattened weight-tables can be appended to containers")

        if hist._labels is None:
            print("skip appending {}".format(hist._filename))
            return

        # use nx to check for initialisation
        if self._nx is None:  # not initialised, copy over all entries

            self._nx = hist.nx
            self._labels = hist._labels
            self._xval = hist._xval
            self._ioverflow = hist._ioverflow

            #@todo: filename & neval are not critical, can use a
            # normal python list instead of a numpy array
            self._filename = np.array([hist._filename])
            self._neval = np.array([hist._neval])

            # yval & yerr are the critical arrays!
            # use pre-allocation for efficiency
            # index ordering is not optimal for accumulating histo's
            # but for later processing we'll be looping over the runs
            # all the time, so we choose to have the run-index the last
            if self._buffer_size is not None:
                # we have requested pre-allocation!
                (ncols, nrows) = hist._yval.shape
                # initialise zero array of correct size
                self._yval = np.zeros((ncols, nrows, self._buffer_size), dtype='float')
                self._yerr = np.zeros((ncols, nrows, self._buffer_size), dtype='float')
                self._mask = np.ones((ncols, nrows, self._buffer_size), dtype='int')
                self._mask *= 2  # means invalid
                # copy over 
                if self._size != 0:
                    raise ValueError('first append and non-zero size?!')
                self._yval[:, :, self._size] = hist._yval[:, :]
                self._yerr[:, :, self._size] = hist._yerr[:, :]
                self._mask[:, :, self._size] = np.logical_or(np.logical_not(np.isfinite(hist._yval[:, :])), np.logical_not(np.isfinite(hist._yerr[:, :]))).astype(int) * 2
                # self._mask[:, :, self._size] = np.isnan(hist._yerr[:, :]).astype(int) * 2
                self._size += 1
            else:
                # the inefficient implementation...
                self._yval = np.expand_dims(hist._yval, axis=2)
                self._yerr = np.expand_dims(hist._yerr, axis=2)

        else:  # already initialised
            # check for consistency
            if self._nx != hist.nx:
                raise ValueError("nx mismatch")
            if self._labels != hist._labels:
                raise ValueError("labels mismatch")
            if not np.array_equal(self._xval, hist._xval):
                raise ValueError("xval mismatch")
            if self._ioverflow != hist._ioverflow:
                raise ValueError("ioverflow mismatch")

            self._filename = np.append(self._filename, hist.filename)
            self._neval = np.append(self._neval, hist._neval)

            if self._buffer_size is not None:
                if self._size == self._buffer_size:
                    raise ValueError('tried to append more runs than pre-allocated')
                self._yval[:, :, self._size] = hist._yval[:, :]
                self._yerr[:, :, self._size] = hist._yerr[:, :]
                self._mask[:, :, self._size] = np.logical_or(np.logical_not(np.isfinite(hist._yval[:, :])), np.logical_not(np.isfinite(hist._yerr[:, :]))).astype(int) * 2
                # self._mask[:, :, self._size] = np.isnan(hist._yerr[:, :]).astype(int) * 2
                self._size += 1                
            else:
                # the inefficient implementation...
                self._yval = np.append(self._yval, np.expand_dims(hist._yval, axis=2), axis=2)
                self._yerr = np.append(self._yerr, np.expand_dims(hist._yerr, axis=2), axis=2)

        if self._buffer_size is None:
            # append always resets the mask
            # and makes it dynamically grow
            self._mask = np.logical_or(np.logical_not(np.isfinite(self._yval)), np.logical_not(np.isfinite(self._yerr))).astype(int) * 2
            # self._mask = np.isnan(self._yval).astype(int) * 2


    def _recursive_k_weights(self, i_row, i_col, wgts_in):
        n_files = len(wgts_in)        
        wgts_out = [0.] * n_files

        leafs = [ [] for i in range(n_files) ]
        ntot = [0.] * n_files
        fails = []
        for i_leaf in range(n_files):
            # trimmed data have no leafs
            if self._mask[i_row, i_col, i_leaf] > 0:
                continue
            # find the target for `i_leaf`
            target = i_leaf
            i_chk = 0
            while self._mask[i_row, i_col, target] < 0:
                if i_chk < n_files:
                    i_chk += 1
                else:
                    fails.append(i_leaf)
                    break
                target = - self._mask[i_row, i_col, target] - 1
            if self._mask[i_row, i_col, target] > 0:
                raise ValueError('leaf terminating on a trimmed entry?!')
            leafs[target].append(i_leaf)
            ntot[target] += self._neval[i_leaf]

        for i_fail in fails:
            print('recursive_k_weights failed!')
            print(' > file {}, [row,col] = [{},{}] '.format(self._filename[i_fail], i_row, i_col))
            target = i_fail
            i_chk = 0
            while self._mask[i_row, i_col, target] < 0:
                if i_chk < n_files:
                    i_chk += 1
                else:
                    break
                print(" > file {}[{}], mask = {}".format(self._filename[target], target, self._mask[i_row, i_col, target]))
                target = - self._mask[i_row, i_col, target] - 1

        for i_target in range(n_files):
            # take the weight in `i_target` and distribute them
            for i_leaf in leafs[i_target]:
                wgts_out[i_leaf] = wgts_in[i_target] * self._neval[i_leaf] / ntot[i_target]

        return np.array(wgts_out)


#     def _recursive_k_weights(self, i_row, i_col, wgts_in):
#         wgts_out = np.zeros(wgts_in.shape)
#         n_files = len(wgts_in)
# 
#         # recursively search where the runs have been merged into
#         # do this once at the beginning to save time (one inner loop)
#         targets = - np.ones(wgts_in.shape)
#         for i_file in range(n_files):
#             # trimmed data have target = -1
#             if self._mask[i_row, i_col, i_file] > 0:
#                 continue
#             # find the target for `i_file`
#             target = i_file
#             while self._mask[i_row, i_col, target] < 0:
#                 target = - self._mask[i_row, i_col, target] - 1
#             if self._mask[i_row, i_col, target] > 0:
#                 raise ValueError('leaf terminating on a trimmed entry?!')
#             targets[i_file] = target
# 
#         for i_file in range(n_files):
#             # find all the leafs that were merged into current `i_file`
#             i_leafs = set(np.where(targets == i_file)[0])            
#             # compute ntot needed for the unweighted partial weight
#             ntot = 0
#             for i_leaf in i_leafs:
#                 ntot += self._neval[i_leaf]
#             # take the weight in `i_file` and distribute them
#             for i_leaf in i_leafs:
#                 wgts_out[i_leaf] = wgts_in[i_file] * self._neval[i_leaf] / ntot
# 
#         return wgts_out


    def _merge_bin(self, i_row, i_col, skip_wgt=False):
        yval, yerr = 0., 0.
        ntot = 0
        n_files = self._yval.shape[2]
        if self._buffer_size is not None:
            n_files = self._size

        # init weight array
        wgts = np.zeros(n_files) if self._qwgt else None

        for i_file in range(n_files):
            # skip trimmed data points
            if self._mask[i_row, i_col, i_file] != 0:
                continue
            ntot += self._neval[i_file]
            yval += self._neval[i_file] * self._yval[i_row, i_col, i_file]
            yerr += self._neval[i_file] * self._yval[i_row, i_col, i_file]**2 \
                + (self._neval[i_file] * self._yerr[i_row, i_col, i_file])**2
            if self._qwgt: wgts[i_file] = self._neval[i_file]

        if ntot != 0:
            yval /= ntot
            if self._qwgt and not skip_wgt: 
                wgts /= ntot
                wgts = self._recursive_k_weights(i_row, i_col, wgts)
            yerr = math.sqrt( yerr - ntot * yval**2 ) / ntot

        return (yval, yerr, wgts)


    def _merge_weighted_bin(self, i_row, i_col, skip_wgt=False):
        yval, yerr = 0., 0.
        n_files = self._yval.shape[2]
        if self._buffer_size is not None:
            n_files = self._size

        # init weight array
        wgts = np.zeros(n_files) if self._qwgt else None

        for i_file in range(n_files):
            # skip trimmed data points
            if self._mask[i_row, i_col, i_file] != 0:
                # if self._mask[i_row, i_col, i_file] == 1: print("* [{}] trimmed: (yval, yerr) = ({}, {})".format(i_file, self._yval[i_row, i_col, i_file], self._yerr[i_row, i_col, i_file]))
                # if self._mask[i_row, i_col, i_file] == 2: print("* [{}] invalid: (yval, yerr) = ({}, {})".format(i_file, self._yval[i_row, i_col, i_file], self._yerr[i_row, i_col, i_file]))
                # if self._mask[i_row, i_col, i_file] <  0: print("* [{}]  -> {}:  (yval, yerr) = ({}, {})".format(i_file, -1-self._mask[i_row, i_col, i_file], self._yval[i_row, i_col, i_file], self._yerr[i_row, i_col, i_file]))
                continue
            # skip zero's to avoid nan's
            if self._yerr[i_row, i_col, i_file] == 0.:
                continue
            weight = 1. / self._yerr[i_row, i_col, i_file]**2
            # print("_merge_weighted_bin", i_row, i_col, i_file, weight)
            yval += self._yval[i_row, i_col, i_file] * weight
            yerr += weight
            if self._qwgt: wgts[i_file] = weight

        if yerr != 0.:
            yval /= yerr
            if self._qwgt and not skip_wgt:
                wgts /= yerr
                wgts = self._recursive_k_weights(i_row, i_col, wgts)
            yerr = math.sqrt(1. / yerr)

        return (yval, yerr, wgts)


    def merge(self, weighted=False):
        # init an empty histogram with correct nx
        result = NNLOJETHistogram(nx=self.nx)

        result._labels = self._labels
        result._xval = self._xval
        result._ioverflow = self._ioverflow

        if self._yval is None:
            raise ValueError('tried to merge a container with no data (zero files?)')

        (n_rows, n_cols, n_files) = self._yval.shape
        if self._buffer_size is not None:
            n_files = self._size

        # init result histogram
        result._neval = np.sum(self._neval)
        result._yval = np.zeros((n_rows, n_cols), dtype='float')
        result._yerr = np.zeros((n_rows, n_cols), dtype='float')

        if self._qwgt:
            result._files_wgt = self._filename
            result._wgt = np.zeros((n_rows, n_cols, n_files), dtype='float')

        for i_row in range(n_rows):
            for i_col in range(n_cols):

                if weighted:
                    ( result._yval[i_row, i_col], 
                      result._yerr[i_row, i_col], 
                      wgts ) = self._merge_weighted_bin(i_row, i_col)

                else:
                    ( result._yval[i_row, i_col], 
                      result._yerr[i_row, i_col], 
                      wgts ) = self._merge_bin(i_row, i_col)

                if self._qwgt:
                    result._wgt[i_row, i_col, :] = wgts[:]
        
        return result


    def unmask(self):
        # only reset trimming flags but not invalid (=2)
        # negative values come from k-merging and that is irreversible
        self._mask[self._mask == 1] = 0
        # if self._buffer_size is not None:
        #     # simplify our life a bit and introduce mask=2 for
        #     # unset buffer entries and only reset self._mask == 1?
        #     self._buffer_size[:, :, self._size:] = 1


    def mask_outliers(self, thresh=4., maxfrac=0.05):
        if thresh <= 0.:
            raise ValueError("non-positive threshold value: {}".format(thresh))
        if maxfrac <= 0. or maxfrac >= 1.:
            raise ValueError("invalid maxfrac range: {}".format(maxfrac))

        self.unmask()

        # what about running this *after* optimise_k ?! (mask[i,j,k] < 0)

        (n_rows, n_cols, n_files) = self._yval.shape
        if self._buffer_size is not None:
            n_files = self._size

        for i_row in range(n_rows):
            for i_col in range(n_cols):
                dyn_thresh = thresh
                # mask = nnlojet_algo.is_outlier_doubleMAD(self._yval[i_row, i_col, :], dyn_thresh).astype(int)
                mask = nnlojet_algo.is_outlier_IQR(self._yval[i_row, i_col, :], dyn_thresh).astype(int)
                cutfrac = np.count_nonzero(mask) / n_files
                # print("#", i_row, i_col, dyn_thresh, cutfrac, maxfrac, file=sys.stderr)
                while cutfrac > maxfrac:
                    dyn_thresh *= 1.1
                    if dyn_thresh > 100.:
                        # unmask this bin for safety?
                        break
                    # mask = nnlojet_algo.is_outlier_doubleMAD(self._yval[i_row, i_col, :], dyn_thresh).astype(int)
                    mask = nnlojet_algo.is_outlier_IQR(self._yval[i_row, i_col, :], dyn_thresh).astype(int)
                    cutfrac = np.count_nonzero(mask) / n_files
                    # print(">", i_row, i_col, dyn_thresh, cutfrac, maxfrac, file=sys.stderr)
                # print("#", i_row, i_col, dyn_thresh, cutfrac, maxfrac, file=sys.stderr)
                self._mask[i_row, i_col, mask == True] = 1


    def optimise_k(self, **kwargs):

        # parse arguments for the merge settings
        maxdev_unwgt = kwargs.get('maxdev_unwgt', None)
        maxdev_steps = kwargs.get('maxdev_steps', None)
        nsteps = kwargs.get('nsteps', None)

        # set up params
        if nsteps is None: maxdev_steps = None
        if maxdev_unwgt is None and maxdev_steps is None: return
            # raise ValueError("optimise_k termination condition missing")

        # reference for termination is the unweighted combination
        if maxdev_unwgt is not None:
            ref_hist = self.merge()

        (n_rows, n_cols, n_files) = self._yval.shape
        if self._buffer_size is not None:
            n_files = self._size

        for i_row in range(n_rows):
            for i_col in range(n_cols):
                # to do this for each bin, need to copy neval
                # we can also have different maskings for the bins
                neval = np.copy(self._neval)

                # we keep track of the history to use it as a termination condition
                history = []

                # print("------", i_row, i_col)

                n_unmasked = neval.size - np.count_nonzero(self._mask[i_row, i_col, :])

                while n_unmasked > 1:

                    # do a weighted average over all runs
                    (yval_wgt, yerr_wgt, wgts) = self._merge_weighted_bin(i_row, i_col, skip_wgt=True)
                    # save history
                    history.append((yval_wgt, yerr_wgt))

                    term_cond1 = False
                    term_cond2 = False

                    # --- termination condition 1 ---
                    # check for `maxdev_unwgt` * sigma compatibility to unweighted reference
                    if maxdev_unwgt is not None:
                        delta = abs(yval_wgt - ref_hist._yval[i_row, i_col])
                        sigma = math.sqrt(yerr_wgt**2 + ref_hist._yerr[i_row, i_col]**2)
                        term_cond1 = delta < maxdev_unwgt * sigma
                        # print(">>>", yval_wgt, yerr_wgt, ref_hist._yval[i_row, i_col], ref_hist._yerr[i_row, i_col], delta, sigma, delta/sigma)

                    # --- termination condition 2 ---
                    # check for `maxdev_steps` * sigma compatibility the last `nsteps` steps
                    if maxdev_steps is not None and len(history) >= nsteps:
                        term_cond2 = True
                        for istep in range(-1,-nsteps-1,-1):
                            for jstep in range(istep-1,-nsteps-1,-1):
                                delta = abs(history[istep][0] - history[jstep][0])
                                sigma = math.sqrt(history[istep][1]**2 + history[jstep][1]**2)
                                term_cond2 = term_cond2 and delta < maxdev_steps * sigma

                    # --- termination of merges ---
                    if term_cond1 or term_cond2:
                        break

                    # print("# unmasked:", n_unmasked)
                    sort_index = np.argsort(neval)

                    # init indices
                    ilow = 0
                    iupp = neval.size - 1
                    low = sort_index[ilow]
                    upp = sort_index[iupp]

                    # print("neval:\n", neval)
                    # print("mask:\n", self._mask[i_row, i_col, :])
                    # print("sort_index:\n", sort_index)
                    # input("Press Enter to continue...")

                    # loop over pairs
                    while ilow < iupp:

                        # skip invalid lower
                        while ilow < neval.size - 1:
                            low = sort_index[ilow]
                            if (self._mask[i_row, i_col, low] == 0 and
                                 neval[low] != 0):
                                break
                            ilow += 1

                        # skip invalid upper
                        while iupp > 0:
                            upp = sort_index[iupp]
                            if (self._mask[i_row, i_col, upp] == 0 and 
                                neval[upp] != 0):
                                break
                            iupp -= 1

                        if ilow >= iupp:
                            break

                        # combine two runs (unweighted)
                        ntot_pair = neval[low] + neval[upp]
                        yval_pair = (
                                + neval[low] * self._yval[i_row, i_col, low]
                                + neval[upp] * self._yval[i_row, i_col, upp]
                            ) / ntot_pair
                        yerr_pair = math.sqrt(
                                + neval[low] * self._yval[i_row, i_col, low]**2
                                + neval[upp] * self._yval[i_row, i_col, upp]**2
                                + (neval[low] * self._yerr[i_row, i_col, low])**2
                                + (neval[upp] * self._yerr[i_row, i_col, upp])**2
                                - ntot_pair * yval_pair**2
                            ) / ntot_pair

                        # print("pair", low, upp, neval[low], neval[upp], ">>>", ntot_pair, yval_pair, yerr_pair)

                        # merge low into upp
                        self._mask[i_row, i_col, low] = -(upp + 1)  # upp can be zero: shift!
                        self._yval[i_row, i_col, low] = np.nan
                        self._yval[i_row, i_col, upp] = yval_pair
                        self._yerr[i_row, i_col, low] = np.nan
                        self._yerr[i_row, i_col, upp] = yerr_pair
                        neval[low] = 0
                        neval[upp] = ntot_pair

                        # move to next pair
                        ilow += 1
                        iupp -= 1

                    n_unmasked = neval.size - np.count_nonzero(self._mask[i_row, i_col, :])






def main():

    pass


if __name__ == "__main__":
    main()
