#!/usr/bin/env python3

import sys
import os
import argparse
import configparser
import glob
import re
import ast
import multiprocessing as mp

from nnlojet_util import NNLOJETHistogram, NNLOJETContainer


class Task_part():
    # merging options and default values
    _trim_threshold     = 4.0
    _trim_max_frac      = 0.1
    _kscan_maxdev_unwgt = None  # needs another safety feature (inactive by default)
    _kscan_nsteps       = 2
    _kscan_maxdev_steps = 0.5
    _merge_weighted     = True

    def __init__(self, **kwargs):
        self._nx = kwargs.get('nx', 3)
        self._files = kwargs.get('files', None)
        self._outfile = kwargs.get('outfile', None)
        self._qweights = kwargs.get('weights', False)
        self._columns = kwargs.get('columns', None)
        self._rebin = kwargs.get('rebin', None)

    def __str__(self):
        return '{} [{} file(s)]'.format(self._outfile, len(self._files))

    def __call__(self):
        # container = NNLOJETContainer()
        # pre-allocate the internal arrays for performance
        container = NNLOJETContainer(size=len(self._files), weights=self._qweights)
        for file in self._files:
            try:
                container.append(NNLOJETHistogram(nx=self._nx, filename=file, columns=self._columns, rebin=self._rebin))
            except ValueError as e:
                print(e)
                print("error reading file:", file)

        # start the merge
        
        outfile_unwgt = self._outfile[:-4] + '.UNWGT.dat'
        hist = container.merge()
        hist.write_to_file(outfile_unwgt)

        outfile_wgt = self._outfile[:-4] + '.WGT.dat'
        hist = container.merge(weighted=True)
        hist.write_to_file(outfile_wgt)

        # trim outliers
        if self._trim_threshold is not None:
            container.mask_outliers(self._trim_threshold, self._trim_max_frac)

            outfile_trim = self._outfile[:-4] + '.TRIM.dat'
            hist = container.merge()
            hist.write_to_file(outfile_trim)

        # k-optimisation
        container.optimise_k(maxdev_unwgt=self._kscan_maxdev_unwgt, nsteps=self._kscan_nsteps, maxdev_steps=self._kscan_maxdev_steps)
        hist = container.merge(weighted=self._merge_weighted)
        hist.write_to_file(self._outfile)


def process_parts_queue(parts_queue):
    while True:
        task = parts_queue.get()
        if task is None:  # the poison pill
            print(' > terminating [worker-id: {}]...'.format(os.getpid()))
            parts_queue.task_done()
            break
        print(' > processing {} [worker-id: {}]...'.format(task, os.getpid()))
        task()
        parts_queue.task_done()


def read_APPLfast(wgt_file):
    with open(wgt_file, 'rt') as file:
        hist_acc = NNLOJETHistogram()
        nx = 0
        for line in file:
            if line.startswith("#"):  # skip comment lines
                # parse nx
                m = re.search(r'#nx=([0-9]+).*', line)
                if m: nx = int(m.groups()[0])
                continue
            array = line.split()
            datfile = array.pop(0)  # first entry is the filename
            hist = NNLOJETHistogram(nx=nx, filename=datfile)
            hist.multiply_weights(array)
            hist_acc = hist_acc + hist
        print(hist_acc)


def obs_generator(obs_list):
    cross_pattern = re.compile(r'.*cross.*')
    for obs_line in obs_list:
        obs_parse = [ it.strip() for it in obs_line.split('>') ]
        if len(obs_parse) == 1:
            obs_in = obs_parse[0]
            obs_out = obs_parse[0]
        elif len(obs_parse) == 2:
            obs_in = obs_parse[0]
            obs_out = obs_parse[1]
        else:
            raise ValueError('invalid observables specification: {}'.format(obs_line))
        nx = 0 if re.match(cross_pattern, obs_in) else 3
        yield (obs_line, obs_in, obs_out, nx)


def main():
    # Give help description
    parser = argparse.ArgumentParser(description='Merge histogram files')
    # config file
    parser.add_argument('-C', '--config', action='store', default='combine.ini',
                        help='configuration file for the combine')
    # multi-process
    parser.add_argument('-j', '--jobs', type=int, nargs='?', action='store', default='1',
                        help='Specifies the number of jobs to run simultaneously.')
    # read in APPLfast weight table
    parser.add_argument('--APPLfast', action='store', default=None,
                        help='Read in an APPLfast weight table and combine.')
    # Parse the input arguments!
    args = parser.parse_args()

    config = configparser.ConfigParser(allow_no_value=True, delimiters=('=', ':'), comment_prefixes=('#',), inline_comment_prefixes=('#',), empty_lines_in_values=False)
    config.optionxform = lambda option: option  # do not convert to lower case
    config.read(args.config)

    # read in APPLfast weight table?
    wgt_file = args.APPLfast
    if wgt_file is not None:
        read_APPLfast(wgt_file)
        return

    cross_pattern = re.compile(r'.*cross.*')
    # if we have dots in the observable this is too naive
    # obs_pattern = re.compile(r'.*\.(.+)\.s[0-9]+\.dat')
    obs_pattern = re.compile(r'.*?/?([^./]+\.){3}([^/]+)\.s[0-9]+\.dat')

    outdir_parts = config.get('Paths', 'out_dir') + '/Parts'
    os.makedirs(outdir_parts, exist_ok=True)

    outdir_final = config.get('Paths', 'out_dir') + '/Final'
    os.makedirs(outdir_final, exist_ok=True)

    # recursive search for histogram files?
    qrecursive = config.getboolean('Options', 'recursive', fallback=False)
    print("recursive search: {}".format('ON' if qrecursive else 'OFF'))

    # output weight data for the combined files?
    qweights = config.getboolean('Options', 'weights', fallback=False)
    print("weights: {}".format('ON' if qweights else 'OFF'))

    # restrict the first `columns` columns?
    columns = config.get('Options', 'columns', fallback=None)
    if columns is not None:
        columns = ast.literal_eval(columns)
        print("only columns: {}".format(columns))


    # read optional trim settings (default: on)
    trim = config.get('Options', 'trim', fallback=None)
    if trim is not None:
        trim = ast.literal_eval(trim)
        if type(trim) is bool:
            if not trim:
                Task_part._trim_threshold = None
                Task_part._trim_max_frac  = None
        elif type(trim) is int or type(trim) is float:
            Task_part._trim_threshold = trim
            Task_part._trim_max_frac  = None
        elif type(trim) is list or type(trim) is tuple:
            if len(trim) == 2:
                Task_part._trim_threshold = trim[0]
                Task_part._trim_max_frac  = trim[1]
            else:
                print("trim list invalid: {}".format(trim))
                raise
        else:
            print("trim option invalid: {}".format(trim))
            raise
    # digest settings
    print("trim: ({},{})".format(Task_part._trim_threshold, Task_part._trim_max_frac))

    # read optional flag to control weighted vs. unweighted average
    Task_part._merge_weighted = config.get('Options', 'weighted', fallback=True)
    print("weighted: {}".format(Task_part._merge_weighted))

    # read optional k-optimisation settings
    kscan = config.get('Options', 'k-scan', fallback=None)
    if kscan is not None:
        kscan = ast.literal_eval(kscan)
        if type(kscan) is bool:
            if not kscan:
                Task_part._kscan_maxdev_unwgt = None
                Task_part._kscan_nsteps       = None
                Task_part._kscan_maxdev_steps = None
        elif type(kscan) is int or type(kscan) is float:
            Task_part._kscan_maxdev_unwgt = kscan
            Task_part._kscan_nsteps       = None
            Task_part._kscan_maxdev_steps = None
        elif type(kscan) is list or type(kscan) is tuple:
            if len(kscan) == 2:
                Task_part._kscan_maxdev_unwgt = None
                Task_part._kscan_nsteps       = kscan[0]
                Task_part._kscan_maxdev_steps = kscan[1]
            elif len(kscan) == 3:
                Task_part._kscan_maxdev_unwgt = kscan[0]
                Task_part._kscan_nsteps       = kscan[1]
                Task_part._kscan_maxdev_steps = kscan[2]
            else:
                print("k-scan list invalid: {}".format(kscan))
                raise
        else:
            print("k-scan option invalid: {}".format(kscan))
            raise
    # digest settings
    print("k-scan: ({},{},{})".format(Task_part._kscan_maxdev_unwgt, Task_part._kscan_nsteps, Task_part._kscan_maxdev_steps))

    # set up for parallel merge
    nworkers = args.jobs
    if nworkers is None:
        nworkers = mp.cpu_count()
    print("running the merge on {} processes".format(nworkers))
    parts_queue = mp.JoinableQueue()
    for i in range(nworkers):
        worker = mp.Process(target=process_parts_queue, args=(parts_queue,))
        worker.start()


    #@todo: checks
    # do all Part directories exist?
    # are there name clashes in the aliases?
    # are there name clashes between Parts & Merge?

    # check if we should auto-scan for observables
    obs_list = config.options('Observables')
    if 'ALL' in obs_list:
        print("automatically scanning for observables...")
        # obs_list = []  # start with an empty list
        obs_list.remove('ALL')
        for pt in config.options('Parts'):
            if qrecursive:
                files = glob.glob(config.get('Paths', 'raw_dir') + '/' + pt + '/**/*.dat', recursive=True)
            else:
                files = glob.glob(config.get('Paths', 'raw_dir') + '/' + pt + '/*.dat')
            for file in files:
                m = re.search(obs_pattern, file)
                if m:
                    obs = m.groups()[-1]
                    if not obs in obs_list:  # unique entries
                        obs_list.append(obs)
                else:
                    print("couldn't extract observable name from file: {}".format(file))
        print("found observables: {}".format(obs_list))


    print("""
        
///////////
// Parts //
///////////
""")

    for (obs_line, obs_in, obs_out, nx) in obs_generator(obs_list):
        print("processing observable {}... (nx = {} bins)".format(obs_out, nx))

        # combine the different parts
        for pt in config.options('Parts'):
            if qrecursive:
                files = glob.glob(config.get('Paths', 'raw_dir') + '/' + pt + '/**/*.' + obs_in + '.*.dat', recursive=True)
            else:
                files = glob.glob(config.get('Paths', 'raw_dir') + '/' + pt + '/*.' + obs_in + '.*.dat')
            mappt = config.get('Parts', pt)
            if mappt is None:
                mappt = pt
            outfile = outdir_parts + '/' + mappt + '.' + obs_out + '.dat'

            rebin = config.get('Observables', obs_line, fallback=None)
            if rebin is not None:
                rebin = ast.literal_eval(rebin)
                #rebin = [ float(i) for i in rebin ]
                #print('rebin = ', rebin)

            task = Task_part(nx=nx, files=files, outfile=outfile, weights=qweights, columns=columns, rebin=rebin)
            print(" > submitting {}".format(outfile))
            parts_queue.put(task)

    # add poison pills for each worker
    for i in range(nworkers):
        parts_queue.put(None)
    parts_queue.close()
    # wait until all Parts have been processed
    parts_queue.join()


    # merge 
    if config.has_section('Merge'):
        print("""

///////////
// Merge //
///////////
""")
        for (obs_line, obs_in, obs_out, nx) in obs_generator(obs_list):
            print("processing observable {}... (nx = {} bins)".format(obs_out, nx))
            for mrg in config.options('Merge'):
                parts = config.get('Merge', mrg)
                if int("|" in parts) + int("&" in parts) + int("+" in parts) > 1:
                    raise ValueError("not allowed to mix '|', '&' or '+' in a merge:", parts)

                elif "+" in parts:
                    parts = list(map(str.strip, parts.split('+')))
                    print(" > Merge[+]: {} = {}...".format(mrg, parts))
                    outfile = outdir_parts + '/' + mrg + '.' + obs_out + '.dat'
                    hist = NNLOJETHistogram()
                    for mappt in parts:
                        ptfile = outdir_parts + '/' + mappt + '.' + obs_out + '.dat'
                        hist = hist + NNLOJETHistogram(nx=nx, filename=ptfile, weights=qweights)
                    hist.write_to_file(outfile)

                elif "|" in parts:
                    parts = list(map(str.strip, parts.split('|')))
                    print(" > Merge[|]: {} = {}...".format(mrg, parts))
                    outfile = outdir_parts + '/' + mrg + '.' + obs_out + '.dat'
                    hist = NNLOJETHistogram()
                    for mappt in parts:
                        ptfile = outdir_parts + '/' + mappt + '.' + obs_out + '.dat'
                        hist = hist.overwrite(NNLOJETHistogram(nx=nx, filename=ptfile, weights=qweights))
                    hist.write_to_file(outfile)
    
                elif "&" in parts:
                    parts = list(map(str.strip, parts.split('&')))
                    print(" > Merge[&]: {} = {}...".format(mrg, parts))
                    outfile = outdir_parts + '/' + mrg + '.' + obs_out + '.dat'
    
                    container = NNLOJETContainer(size=len(parts), weights=qweights)
                    for pt in parts:
                        ptfile = outdir_parts + '/' + pt + '.' + obs_out + '.dat'
                        try:
                            container.append(NNLOJETHistogram(nx=nx, filename=ptfile))
                        except ValueError as e:
                            print(e)
                            print("error reading file:", ptfile)
                    # weighted average
                    hist = container.merge(weighted=True)
                    hist.write_to_file(outfile)
    
                else:
                    raise ValueError("couldn't find '|' nor '&' in the merge:", parts)
 

    print("""

///////////
// Final //
///////////
""")
    for (obs_line, obs_in, obs_out, nx) in obs_generator(obs_list):
        print("processing observable {}... (nx = {} bins)".format(obs_out, nx))
        # assemble final histograms
        for fin in config.options('Final'):
            parts = list(map(str.strip, config.get('Final', fin).split('+')))
            print(" > Final: {} = {}...".format(fin, parts))
            outfile = outdir_final + '/' + fin + '.' + obs_out + '.dat'
            hist = NNLOJETHistogram()
            for mappt in parts:
                ptfile = outdir_parts + '/' + mappt + '.' + obs_out + '.dat'
                hist = hist + NNLOJETHistogram(nx=nx, filename=ptfile, weights=qweights, recursive_weights=qweights)
            hist.write_to_file(outfile)
            # write to the file
            if qweights:
                with open(outdir_final + '/' + fin + '.' + obs_out + '.APPLfast.txt', 'wt') as fout:
                    print(hist.to_APPLfast_weight(), file=fout)
            

if __name__ == "__main__":
    main()
