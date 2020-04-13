#!/usr/bin/env python3
#-*- coding:utf-8 -*-
#
########################################################################
#
# Plot the closure of APPLfast versus NNLOJET
# from ensemble of one-to-one comparisons
#
# Created by K. Rabbertz, 13.07.2017
# Modified by B. Schillinger, 04.09.2018
# Prepared for python3 by K. Rabbertz, 07.03.2020
#
########################################################################
#
# python2 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3
import argparse
import glob
import os
import re
import string
import sys
import timeit
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.ticker import (FormatStrFormatter, LogFormatter, NullFormatter, ScalarFormatter, AutoMinorLocator, MultipleLocator)
from matplotlib import cm
# We do not want any interactive plotting! Figures are saved to files instead.
# This also avoids the ANNOYANCE of frequently missing Tkinter/tkinter (python2/3) GUI backends!
# To produce scalable graphics for publication use eps, pdf, or svg as file format.
# For this to work we try the Cairo backend, which can do all of these plus the raster format png.
# If this is not usable, we fall back to the Agg backend capable only of png for nice web plots.
#ngbackends = mpl.rcsetup.non_interactive_bk
#print('[fastnnlo_pdfunc]: Non GUI backends are: ', ngbackends)
# 1st try cairo
backend = 'cairo'
usecairo = True
try:
    import cairocffi as cairo
except ImportError:
    try:
        import cairo
    except ImportError:
        usecairo = False
#        print('[fastnnlo_pdfunc]: Can not use cairo backend :-(')
#        print('                   cairocffi or pycairo are required to be installed')
    else:
        if cairo.version_info < (1, 11, 0):
            # Introduced create_for_data for Py3.
            usecairo = False
#            print('[fastnnlo_pdfunc]: Can not use cairo backend :-(')
#            print('                   cairo {} is installed; cairo>=1.11.0 is required'.format(cairo.version))
if usecairo:
    mpl.use('cairo')
else:
    backend = 'agg'
    useagg = True
    try:
        mpl.use(backend, force=True)
        print('[fastnnlo_pdfunc]: Warning! Could not import cairo backend :-( Using agg instead for raster plots only!')
    except:
        useagg = False
        print('[fastnnlo_pdfunc]: Can not use agg backend :-(')
        raise ImportError('[fastnnlo_pdfunc]: Neither cairo nor agg backend found :-( Cannot produce any plots. Good bye!')
    mpl.use('agg')
import matplotlib.pyplot as plt
# numpy
import numpy as np


# Action class to allow comma-separated (or empty) list in options
class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values:
            setattr(namespace, self.dest, values[0].split(','))
        else:
            setattr(namespace, self.dest, [''])


# Some global definitions
_formats = {'eps': 0, 'pdf': 1, 'png': 2, 'svg': 3}
_debug = False
_nmax = 99999 # Default maximum number of files to evaluate

########################################################################

def main():
    ################################################################################
    # Style settings
    # Location of matplotlibrc
    # print mpl.matplotlib_fname()
    # List active style settings
    # print mpl.rcParams
    # Needs matplotlib > 1.3
    # plt.style.use('ggplot')
    # plt.style.use('presentation')
    params = {'legend.fontsize': 'x-large',
              'figure.figsize': (16, 12),
              'mathtext.fontset': "stix",
              'axes.labelsize':  'x-large',
              'axes.titlesize':  'x-large',
              #'axes.linewidth':  2, #increase default linewidth
              'xtick.labelsize': 'x-large',
              'ytick.labelsize': 'x-large'}
    mpl.rcParams.update(params)

    # Start timer
    # just for measuring wall clock time - not necessary
    start_time = timeit.default_timer()

    # Define arguments & options
    parser = argparse.ArgumentParser(epilog='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Positional arguments
    parser.add_argument('table', nargs=1, type=argparse.FileType('r'),
                        help='Example filename, i.e. with seed index, of fastNLO tables to be evaluated. By default the same basename is assumed to match NNLOJET dat and weight files as well as fastNLO log files; figures are stored as tabfilebasename.sclx.ratio|asymm_nfiles.png')
    # Optional arguments
    parser.add_argument('-d', '--datfile', required=False, nargs=1, type=argparse.FileType('r'),
                        help='Example NNLOJET dat file for evaluation.')
    parser.add_argument('-l', '--logfile', required=False, nargs=1, type=argparse.FileType('r'),
                        help='Example fastNLO log file for evaluation. A direct evaluation of table files within this code is possible, but too time consuming running over hundreds of tables.')
    parser.add_argument('-w', '--weightfile', required=False, nargs=1, type=argparse.FileType('r'),
                        help='NNLOJET .txt file containing merging weights for each datfile and observable bin.')
    parser.add_argument('-n', '--nfiles', required=False, nargs=1, type=int, default=[0],
                        help='Number of result pairs to compare. Zero means all available.')
    parser.add_argument('-f', '--filename', default=None, type=str,
                        help='Replace grid file basename by string in output figure name (optional).')
    parser.add_argument('--format', required=False, nargs=1, type=str, action=SplitArgs,
                        help='Comma-separated list of plot formats to use: eps, pdf, png, svg. If nothing is chosen, png is used.')
    parser.add_argument('-s', '--scale', default=1, required=False, nargs='?', type=int, choices=[*range(-12, 0), *range(1, 8)], metavar='[-12,..,-1,1,..7]',
                        help='Scale combination for mu_r, mu_f to be used in comparison. scale=1: Only central scale; scale in [2,7]: Central & add. scale factor combinations up to given index; scale in [-6,-1]: Central & low fixed-scale variations down to given index; scale in [-12,-7]: Central & high fixed-scale variations down to given index.')
    parser.add_argument('--single', action="store_true",
                        help='If --single is chosen, only the scale index given via -s is used.')
    parser.add_argument('--title', default=None, type=str,
                        help='Replace table name as default title by given string.')
    parser.add_argument('-v', '--verbose', action="store_true",
                        help="Increase output verbosity.")
    parser.add_argument('--xlabel', default=None, type=str,
                        help='Replace x axis default label by given string.')
    parser.add_argument('--ylabel', default=None, type=str,
                        help='Replace y axis default label by given string.')

    #
    # Print header
    #
    print("\n###########################################################################################")
    print("# fastnnlo_multigridclosure:")
    print("# Plot statistical evaluation of APPLfast closure vs. NNLOJET")
    print("###########################################################################################\n")

    #
    # Parse arguments
    #
    args = vars(parser.parse_args())

    # Table name
    tabfile = args['table'][0].name
    print('\n')
    print('[fastnnlo_multigridclosure]: Analysing tables: {}'.format(tabfile))
    # Eliminate extensions
    tabfile = re.sub('.gz$', '', tabfile)
    tabfile = re.sub('.tab$', '', tabfile)

    # Take care of the input files
    datfile = None
    logfile = None
    wgtfile = None
    if args['datfile']:
        datfile = args['datfile'][0].name
    if args['logfile']:
        logfile = args['logfile'][0].name
    if args['weightfile']:
        wgtfile = args['weightfile'][0].name
    if not datfile:
        datfile = tabfile + '.dat'
    print('[fastnnlo_multigridclosure]: NNLOJET dat file: {:s}'.format(datfile))
    if not logfile:
        logfile = tabfile + '.log'
    print('[fastnnlo_multigridclosure]: fastNLO log file: {:s}'.format(logfile))
    if wgtfile:
        print('[fastnnlo_multigridclosure]: NNLOJET weight file: {:s}'.format(wgtfile))

    # Number of result pairs to compare
    nfil = args['nfiles'][0]

    # Scale settings
    scaleindex = args['scale']
    sclvar = not args['single'] and scaleindex != 1

    # Plot formats to use
    formats = args['format']
    if formats is None:
        formats = ['png']
    for fmt in formats:
        if fmt not in _formats:
            print('[fastnnlo_multigridclosure]: Illegal format specified, aborted!')
            print('[fastnnlo_multigridclosure]: Format list:', args['format'])
            sys.exit(1)
        elif fmt != 'png' and not usecairo:
            print('[fastnnlo_multigridclosure]: Vector format plots not possible without cairo backend, aborted!')
            print('[fastnnlo_multigridclosure]: Format list:', args['format'])
            sys.exit(2)

    # Plot labelling
    nice_title = args['title']
    nice_xlabel = args['xlabel']
    nice_ylabel = args['ylabel']

    # Verbosity
    verb = args['verbose']
    if verb:
        print('[fastnnlo_multigridclosure]: Using matplotlib version {:s}'.format(mpl.__version__))
        print('                       from location {:s}'.format(mpl.__file__))

    # Decrypt arguments from table name
    tabbase = os.path.basename(tabfile)
    tabargs = tabbase.split(".")

    # Decrypt table name structure (tabfile extensions have been cut here!)
    nlargs = len(tabargs)
    if nlargs == 3:
        proc, jobn, obsv = tabargs
    elif nlargs == 4:
        proc, jobn, kinn, obsv = tabargs
    elif nlargs == 5:
        proc, jobn, kinn, obsv, seed = tabargs
    else:
        print('[fastnnlo_multigridclosure]: ERROR! Unknown table name structure. Aborted!')
        print('[fastnnlo_multigridclosure]: Found {:d} instead of 4-6 point-separated substrings in {:s}.'.format(nlargs,tabbase))
        sys.exit(3)

    print('[fastnnlo_multigridclosure]: NNLOJET process: {:s}'.format(proc))
    print('[fastnnlo_multigridclosure]: NNLOJET job name: {:s}'.format(jobn))
    if nlargs>3:
        print('[fastnnlo_multigridclosure]: NNLOJET kinematics acronym: {:s}'.format(kinn))
    print('[fastnnlo_multigridclosure]: Observable: {:s}'.format(obsv))
    if nlargs>4:
        print('[fastnnlo_multigridclosure]: NNLOJET seed index: {:s}'.format(seed))

    if nlargs != 5:
        print('[fastnnlo_multigridclosure]: ERROR! Individual seed index results are required, but seed index could not be identified! Aborted.')
        sys.exit(4)

    # Extract order/contribution from job type (substring before first '-')
    order = jobn.split('-')[0]

    # Prepare result arrays which are not used by Read_XS() function, but later
    xm = []      # bin "center"
    ntab = 0
    xs_fnlt = []  # fastNLO results
    nlin = 0
    weights = []  # Weight factors

    # Read binning and cross sections from NNLOJET dat files
    # Insist on constant seed index length; otherwise log files might not be present when adding grids from external production with other seed ranges
    seedglob = (len(seed)-1) * '?'
    datglob = datfile[:-(len(seed)+len('dat'))] + seedglob + '.dat'
    print('[fastnnlo_multigridclosure]: dat file glob is: {}'.format(datglob))
    datfiles = glob.glob(datglob)
    if not datfiles:
        print('[fastnnlo_multigridclosure]: ERROR! No NNLOJET dat files matching {} found! Aborted.'.format(datglob))
        sys.exit(5)
    datfiles.sort()

    xs_nnlo = None  # will be updated in case only one scale is being looked at
    # otherwise xs_nnlo_list will be used (and all scales investigated)
    xs_nnlo_list = None  # see comment above (just the other way round)
    xs_fnll = None  # same logic as above
    xs_fnll_list = None  # same logic as above

    # Analyse chosen scale indices
    if scaleindex < -6:
        sclind = -5 - scaleindex
    elif scaleindex < 0:
        sclind = 1 - scaleindex
    else:
        sclind = scaleindex
    if not sclvar:
        if sclind in np.arange(1, 8, 1):  # if single valid fscl is chosen
            xl, xu, ndat, xs_nnlo, nobs, seeds = Read_XS(datfiles, sclind, nfil, verb)
    else: # option --single not used
        xs_nnlo_list = []  # will become a list of arrays
        for fscl in np.arange(1, sclind+1, 1):
            # take bin bounds from first iteration (fscl=1), they stay the same
            if fscl == 1:
                xl, xu, ndat, xs_nnlo_item, nobs, seeds = Read_XS(datfiles, fscl, nfil, verb)
                xs_nnlo_list.append(xs_nnlo_item)
            else:
                xs_nnlo_item = Read_XS(datfiles, fscl, nfil, verb)[3]
                xs_nnlo_list.append(xs_nnlo_item)
                # xs_nnlo_list contains 7 arrays now

#    if xs_nnlo is not None and verb:
#        print('[fastnnlo_multigridclosure]: shape of xs_nnlo: ', np.shape(xs_nnlo))
#    if xs_nnlo_list is not None and verb:
#        print('[fastnnlo_multigridclosure]: shape of xs_nnlo_list', np.shape(xs_nnlo_list))

    # Read weights per file per bin from NNLOJET
    if wgtfile:
        wgtnams = np.genfromtxt(wgtfile, dtype=str, usecols=0)
        wgttmps = np.loadtxt(wgtfile, usecols=(list(range(1, nobs+1))))
        if verb:
            print('[fastnnlo_multigridclosure]: wgttmps', wgttmps)
        ntmp = len(wgtnams)
        # Combine to key-value tuple ( name, weights[] ) and sort according to key=name
        wgttup = list(zip(wgtnams, wgttmps))
        wgttup.sort(key=lambda row: (row[0]))
        # Unzip again
        allnames, allweights = list(zip(*wgttup))
        if verb:
            print('[fastnnlo_multigridclosure]: datfile names in weightfile: \n', np.array(allnames), '\n')

        #print 'datfiles array: \n', np.array(datfiles), '\n'
        #print "weights..", weights
        if verb:
            print ('[fastnnlo_multigridclosure]: allnames', allnames)
        for dfile in datfiles:  # does not have to be changed, as choice of sclind makes no difference here
            newna = './'+order+'/'+os.path.basename(dfile)
            if newna in allnames:
                indexlin = allnames.index(newna)
                if verb:
                    print('[fastnnlo_multigridclosure]: Weight file line no. ', indexlin, ' is for ', allnames[indexlin])
                nlin += 1
                if (nfil == 0) or (nfil > 0 and nlin <= nfil):
                    weights.append(allweights[indexlin])
                else:
                    nlin -= 1
                    break
            else:
                print('[fastnnlo_multigridclosure]: ERROR: Did not find matching weight file line! Aborted.')
                sys.exit(6)
        print('[fastnnlo_multigridclosure]: Using %s weight lines.' % nlin)

    # Empty list gets also transformed into np array!
    weights = np.array(weights)

    # Evaluate cross sections from fastNLO tables
    # INFO=0, WARNING=1
    # SetGlobalVerbosity(1)
    #fnlotabs = glob.glob(pord+'/'+scen+'.'+proc+'.'+pord+'-'+ecms+'.???.'+obsv+'.*.tab.gz')
    # fnlotabs.sort()
    # for fnlotab in fnlotabs:
    #    print 'fastNLO table no. ', ntab, ' is ', fnlotab
    #    fnlo = fastNLOLHAPDF(fnlotab)
    #    fnlo.SetLHAPDFFilename('CT14nnlo')
    #    fnlo.SetLHAPDFMember(0)
    #    fnlo.SetMuFFunctionalForm(1); # kScale1=0
    #    fnlo.SetMuRFunctionalForm(1); # kScale2=1
    #    fnlo.CalcCrossSection()
    #    xs_fnla.append(fnlo.GetCrossSection())
    #    ntab += 1
    #    if ntab==_nmax: break
    #print 'Using ', ntab, 'table files.'
    #xs_fnla = np.array(xs_fnla)

    # Evaluate cross sections from pre-evaluated fastNLO tables
    # assume that corresponding logfiles are located in the same directory as the datfiles
    log_list = []

    nd = 0
    for dfile in datfiles:
        nd += 1
        if (nfil == 0) or (nfil > 0 and nd <= nfil):
            log_list.append(dfile[:-(len('dat')+1)]+'.log')
        else:
            break
    fnlologs = np.array(log_list)
    if verb:
        print('[fastnnlo_multigridclosure]: fnlologs: \n', fnlologs, '\n')

    fnlologs.sort()
    if (xs_nnlo is not None):  # Only one specific scale choice sclind is considered
        xs_fnll, nlog = Read_logfile(fnlologs, sclind, seeds, nobs, verb)
    elif (xs_nnlo_list is not None):  # Have list with values for all sclind variations
        xs_fnll_list = []
        for fscl in range(1, sclind+1, 1):
            xs_fnll_item, nlog = Read_logfile(fnlologs, fscl, seeds, nobs, verb)
            xs_fnll_list.append(xs_fnll_item)

    # Check on identical file numbers, either for table or log files, dat files, and weight lines if given
    if ndat == nlog == nlin and ndat*nlog*nlin != 0 and wgtfile:
        print('[fastnnlo_multigridclosure]: OK: Have equal number of dat files, fastNLO results, and weight lines.')
        print('[fastnnlo_multigridclosure]: ndat = ', ndat, ', nlog = ', nlog, ', nlin = ', nlin)
    elif ndat == nlog and ndat*nlog != 0 and not wgtfile:
        print('[fastnnlo_multigridclosure]: OK: Have equal number of dat files and fastNLO results.')
        print('[fastnnlo_multigridclosure]: ndat = ', ndat, ', nlog = ', nlog)
    else:
        print('[fastnnlo_multigridclosure]: ERROR: No matching file found or file number mismatch! ndat = ' +
              str(ndat)+', nlog = '+str(nlog)+', nlin = '+str(nlin))
        sys.exit(7)

    # prepare info_values tuple for plotting function
    info_values = (weights, nobs, ndat, nlin, proc, jobn, order, obsv, args)

    # Use plotting function:
    if xs_fnll is not None:  # should be calculated already, in case only one sclind is investigated
        Statistics_And_Plotting(sclind, xs_fnll, xs_nnlo, info_values)
        # plot here and save the two plots
        print("[fastnnlo_multigridclosure]: Plotting done. (fscl=%s)" % sclind)
    elif xs_fnll_list is not None:  # that list exists if all fscl from 1 to sclind are considered, otherwise it is None
        r_f2nl_list = []
        a_f2nl_list = []
        for fscl in range(1, sclind+1, 1):
            xs_fnll = xs_fnll_list[fscl-1]
            xs_nnlo = xs_nnlo_list[fscl-1]

            # plot here and save plots
            Statistics_And_Plotting(fscl, xs_fnll, xs_nnlo, info_values)
            print("[fastnnlo_multigridclosure]: Plotting done for fscl= %s" % fscl)

    stop_time = timeit.default_timer()
    timediff = stop_time-start_time
    print('[fastnnlo_multigridclosure]: Elapsed time: %s sec = %s min' %
          (timediff, round(timediff/60., 2)))

    exit(0)
    # what are these plots like?
    #    for i_bin, (nnlo, fnlo, ratios, asyms, wgts) in enumerate(zip(xs_nnlo.T, xs_fnll.T, r_f2nl.T, a_f2nl.T, weights.T)):
    #
    #        if i_bin > 0:
    #            exit(2222)
    #        print('NNLOJET cross sections')
    #        bin_dists(nnlo, 'NNLOJET', r'$\sigma$ [pb]', 'blue', 'violet')
    #        print('fastNLO cross sections')
    #        bin_dists(fnlo, 'fastNLO', r'$\sigma$ [pb]', 'blue', 'violet')
    #        print('Ratios')
    #        bin_dists(ratios, 'Ratio', 'fastNLO/NNLOJET', 'blue', 'violet', True,
    #                  rmed_f2nl[i_bin]-10.*riqd_f2nl[i_bin], rmed_f2nl[i_bin]+10.*riqd_f2nl[i_bin])
    #        print('Asymmetries')
    #        bin_dists(asyms, 'Asymmetry', '(fastNLO-NNLOJET)/(fastNLO+NNLOJET)', 'blue', 'violet',
    #                  True, amed_f2nl[i_bin]-10.*aiqd_f2nl[i_bin], amed_f2nl[i_bin]+10.*aiqd_f2nl[i_bin])

################################################################################################
##################################### END OF MAIN ##############################################
################################################################################################

### FUNCTIONS to be used in main() ###

################################################################################
# Define method for histogram statistics
################################################################################
def hist_stats(x, weights, axis=None):

    # In case weight array is empty, create array filled with unity weights
    if not weights.size:
        weights = np.ones_like(x)

    # Unweighted mean and sample variance (not population variance)
    _ave = np.mean(x, axis=axis)
    if x.shape[axis] == 1:
        _var = np.zeros_like(_ave)
    else:
        _var = np.var(x, ddof=1, axis=axis)
    _std = np.sqrt(_var)

    # Useful sums for weighted quantities
    sumw = np.sum(weights, axis=axis)
    sumw2 = np.sum(weights * weights, axis=axis)
    sumwx = np.sum(weights * x, axis=axis)
    sumwx2 = np.sum(weights * x * x, axis=axis)

    # Weighted mean and reliability weighted sample variance
    _wave = []
    _wvar = []
    for i in range(len(sumw)):
        if sumw[i] == 0:
            _wave.append(np.nan)
        else:
            _wave.append(sumwx[i]/sumw[i])

    if x.shape[axis] == 1:
        _wvar = np.zeros_like(_ave)
    else:
        for i in range(len(sumw)):
            if sumw[i] == 0:
                _wvar.append(np.nan)
            else:
                _wvar.append(
                    (sumwx2[i] - sumwx[i]*sumwx[i]/sumw[i]) / (sumw[i] - sumw2[i]/sumw[i]))

    _wstd = np.sqrt(_wvar)

# Median and half interquartile distance
    _med = np.median(x, axis=axis)
    _med_err = np.subtract(
        *np.percentile(x, [75, 25], axis=axis, interpolation='linear'))/2.

    return dict(mean=_ave,
                stddev=_std,
                weighted_mean=_wave,
                weighted_stddev=_wstd,
                median=_med,
                iqd2=_med_err)
################################################################################
# End of method for histogram statistics
################################################################################

################################################################################
# Function to read bin bounds and xs_nnlo from datfiles
################################################################################
def Read_XS(datfiles, fscl, nfil, verb):  # takes list of datfiles (different seeds) and fscl
    # Prepare result arrays
    xl = []      # left bin border
    xu = []      # right bin border
    ndat = 0
    xs_nnlo = []  # NNLOJET results
    seeds = []   # Seed numbers for matching

    print('[fastnnlo_multigridclosure]: Reading datfiles for fscl=%s' % fscl)
    ixscol = 3 + 2 * (fscl-1)
    for datfile in datfiles:
        if verb:
            print('[fastnnlo_multigridclosure]: Datfile no. ', ndat, ' is ', datfile)
        if ndat == 0:
            xl.append(np.loadtxt(datfile, usecols=(0,)))
            xu.append(np.loadtxt(datfile, usecols=(2,)))
        ndat += 1
        if ((nfil == 0) or (nfil > 0 and ndat <= nfil)) and not ndat == _nmax:
            xs_nnlo.append(np.loadtxt(datfile, usecols=(ixscol,)))
            parts = datfile.split(".")
            seed = parts[len(parts)-2]
            seeds.append(seed)
        else:
            ndat -= 1
            break
    print('[fastnnlo_multigridclosure]: Using ', ndat, 'dat files.')
    xl = np.array(xl)
    xu = np.array(xu)
    xs_nnlo = np.array(xs_nnlo)/1000.  # Conversion of fb to pb

    # Determine no. of observable bins
    nobs = xl.size
    print('[fastnnlo_multigridclosure]: Number of observable bins: ', nobs)
    return xl, xu, ndat, xs_nnlo, nobs, seeds
###############################################################################
# End of function to read bin bounds and xs_nnlo from datfiles
###############################################################################

###############################################################################
# Function to read XS from logfile
###############################################################################
def Read_logfile(fnlologs, fscl, seeds, nobs, verb):  # takes list of logfiles and fscl,
    # as well as list of seeds and nobs (number of observable bins)
    nlog = 0
    xs_fnll = []
    for fnlolog in fnlologs:
        if verb:
            print('[fastnnlo_multigridclosure]: fastNLO file no. ', nlog, ' is ', fnlolog)
        if not seeds[nlog] in fnlolog:
            print('[fastnnlo_multigridclosure]: Mismatch in result sort order between NNLOJET and fastNLO. Aborted!')
            print('[fastnnlo_multigridclosure]:seeds[', nlog, '] = ', seeds[nlog], ', fastNLO file = ', fnlolog)
            sys.exit(8)
        else:
            if verb:
                print('[fastnnlo_multigridclosure]: NNLOJET and fastNLO result correctly matched. seed is ', seeds[nlog])

        xs_tmp = np.loadtxt(fnlolog, usecols=(
            6,), comments=['#', ' #', 'C', 'L', 'N'])
        #print "xs_tmp \n", xs_tmp

        indi = (fscl-1)*nobs  # skip lines of lower fscl
        indf = indi + nobs  # last line in logfile which belongs to certain fscl
        xs_sub = xs_tmp[indi:indf]
        xs_fnll.append(xs_sub)
        nlog += 1
        if nlog == _nmax:
            break
    print('[fastnnlo_multigridclosure]: Using ', nlog, 'pre-evaluated table files.')
    xs_fnll = np.array(xs_fnll)
    return xs_fnll, nlog
###############################################################################
# End of function to read XS from logfile
###############################################################################

###############################################################################
# Function for Statistics and Plotting
###############################################################################
# takes pre-calculated values for specific fscl (process variables as tuple)
def Statistics_And_Plotting(fscl, xs_fnll, xs_nnlo, info_values):
    # get required variables from tuple
    weights, nobs, ndat, nlin, proc, jobn, order, obsv, args = info_values

    r_f2nl = xs_fnll/xs_nnlo
    a_f2nl = (xs_fnll - xs_nnlo)/(xs_fnll + xs_nnlo)

    # Ratio statistics
    _ratio_stats = hist_stats(r_f2nl, weights=weights, axis=0)
    rave_f2nl = _ratio_stats['mean']
    rstd_f2nl = _ratio_stats['stddev']
    rwav_f2nl = _ratio_stats['weighted_mean']
    rwst_f2nl = _ratio_stats['weighted_stddev']
    rmed_f2nl = _ratio_stats['median']
    riqd_f2nl = _ratio_stats['iqd2']

    # Asymmetry statistics
    _asymm_stats = hist_stats(a_f2nl, weights=weights, axis=0)
    aave_f2nl = _asymm_stats['mean']
    astd_f2nl = _asymm_stats['stddev']
    awav_f2nl = _asymm_stats['weighted_mean']
    awst_f2nl = _asymm_stats['weighted_stddev']
    amed_f2nl = _asymm_stats['median']
    aiqd_f2nl = _asymm_stats['iqd2']

    # Settings for x-axis of plots
    # For setting small xticks for each bin with only certain labels
    # Set minimum > unity in MultipleLocator
    majorLocator = MultipleLocator(max((nobs//6), 1.1))
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(1)  # set minor tick for each observable bin

    # Ratio plot
    limfs = 'x-large'
    # renamed fig to fig1, so function can return it
    fig1 = plt.figure(figsize=(16, 12))
    ax1 = fig1.gca()
    plt.title(r'Ratio: {} {} for {} using scale {}'.format(
        proc, order, obsv, fscl), fontsize='x-large', fontweight='bold', loc='left')
    plt.xlabel('Observable bin index', horizontalalignment='right',
               x=1.0, verticalalignment='top', y=1.0)
    plt.ylabel('fastNLO/NNLOJET', horizontalalignment='right',
               x=1.0, verticalalignment='bottom', y=1.0)
    plt.axhline(y=1.001, linestyle='--', linewidth=1.0, color='black')
    plt.axhline(y=0.999, linestyle='--', linewidth=1.0, color='black')
    plt.fill_between([0.0, nobs+1.1], 0.999, 1.001, color='black', alpha=0.1)
    # plt.text(33.6, 1.00085, u'+1‰')
    # plt.text(33.7, 0.99885, u'–1‰')
    plt.text(nobs+1.1, +1.00085, u'+1‰')
    # location flexible adjusted to plotting range
    plt.text(nobs+1.1, +0.99885, u'–1‰')

    if nlin > 0:
        dx = 1./3.
    else:
        dx = 1./4.
    xa = np.arange(1-dx, nobs-dx+1.e-6)
    xb = np.arange(1, nobs + 1.e-6)
    xc = np.arange(1+dx, nobs+dx+1.e-6)

    rmean = plt.errorbar(xa, rave_f2nl, yerr=0.*rstd_f2nl, marker='v',
                         linestyle='none', label=r'$\mu$ $\pm$ $\Delta\mu$', color='blue')
    plt.errorbar(xa, rave_f2nl, yerr=rstd_f2nl/np.sqrt(ndat),
                 marker='.', linestyle='none', color='blue')
    if nlin > 0:
        rmwgt = plt.errorbar(xb, rwav_f2nl, yerr=0.*rwst_f2nl, marker='^',
                             linestyle='none', label=r'$\mu_w$ $\pm$ $\Delta\mu_w$', color='violet')
        plt.errorbar(xb, rwav_f2nl, yerr=rwst_f2nl/np.sqrt(nlin),
                     marker='.', linestyle='none', color='violet')
    rmedian = plt.errorbar(xc, rmed_f2nl, yerr=riqd_f2nl, marker='s',
                           linestyle='none', label=r'median $\pm$ IQD/2', color='red')

    # plt.xlim(0.0,34.0)
    plt.xlim(0, nobs+1)
    #plt.xticks(range(0, nobs+2, nobs//6))
    ax1.xaxis.set_major_locator(majorLocator)
    ax1.xaxis.set_major_formatter(majorFormatter)
    ax1.xaxis.set_minor_locator(minorLocator)
    ax1.tick_params(which='both', width=1)
    ax1.tick_params(which='major', length=7)
    ax1.tick_params(which='minor', length=5)
    plt.ylim(0.99, 1.01)
    # plt.ylim(-50, 12000) #to display points for fscl=0 (log6file), was just a test
    if nlin > 0:
        handles = [rmean, rmwgt, rmedian]
    else:
        handles = [rmean, rmedian]
    labels = [h.get_label() for h in handles]
    legend = ax1.legend(handles, labels, title=r'No. of entries/bin = {}'.format(
        ndat), loc='upper left', numpoints=1, handlelength=0)
    legend.get_title().set_fontsize(limfs)

    if args['filename'] is None:
        fignam = proc+'.'+jobn+'.'+obsv+'.'+'scl' + \
            str(fscl)+'.'+'ratio_'+str(ndat)+'.png'
    else:
        fignam = args['filename']+'.scl' + \
            str(fscl)+'.ratio_'+str(ndat)+'.png'
    plt.savefig(fignam)
    #print "Ratio plot saved as: %s" %fignam
    # plt.show()

    # Asymmetry plot
    limfs = 'x-large'
    fig2 = plt.figure(figsize=(16, 12))  # renamed fig to fig2
    ax2 = fig2.gca()
    plt.title(r'Asymmetry: {} {} for {} using scale {}'.format(
        proc, order, obsv, fscl), fontsize='x-large', fontweight='bold', loc='left')
    plt.xlabel('Observable bin index', horizontalalignment='right',
               x=1.0, verticalalignment='top', y=1.0)
    plt.ylabel('(fastNLO-NNLOJET)/(fastNLO+NNLOJET)',
               horizontalalignment='right', x=1.0, verticalalignment='bottom', y=1.0)
    plt.axhline(y=-0.001, linestyle='--', linewidth=1.0, color='black')
    plt.axhline(y=+0.001, linestyle='--', linewidth=1.0, color='black')
    plt.fill_between([0.0, nobs+1.1], -0.001, 0.001, color='black', alpha=0.1)
    # plt.text(33.6, +0.00085, u'+1‰')
    # plt.text(33.7, -0.00115, u'–1‰')
    plt.text(nobs+1.1, +0.00085, u'+1‰')
    # location flexible adjusted to plotting range
    plt.text(nobs+1.1, -0.00115, u'–1‰')

    if nlin > 0:
        dx = 1./3.
    else:
        dx = 1./4.
    xa = np.arange(1-dx, nobs-dx+1.e-6)
    xb = np.arange(1, nobs + 1.e-6)
    xc = np.arange(1+dx, nobs+dx+1.e-6)

    amean = plt.errorbar(xa, aave_f2nl, yerr=0.*astd_f2nl, marker='<',
                         linestyle='none', label='$\mu$ $\pm$ $\Delta\mu$', color='orange')
    plt.errorbar(xa, aave_f2nl, yerr=astd_f2nl/np.sqrt(ndat),
                 marker='.', linestyle='none', color='orange')
    if nlin > 0:
        amwgt = plt.errorbar(xb, awav_f2nl, yerr=0.*awst_f2nl, marker='>',
                             linestyle='none', label=r'$\mu_w$ $\pm$ $\Delta\mu_w$', color='brown')
        plt.errorbar(xb, awav_f2nl, yerr=awst_f2nl/np.sqrt(ndat),
                     marker='.', linestyle='none', color='brown')
    amedian = plt.errorbar(xc, amed_f2nl, yerr=aiqd_f2nl, marker='8',
                           linestyle='none', label=r'median $\pm$ IQD/2', color='green')

    # plt.xlim(0.0,34.0)
    plt.xlim(0, nobs+1)
    ax2.xaxis.set_major_locator(majorLocator)
    ax2.xaxis.set_major_formatter(majorFormatter)
    ax2.xaxis.set_minor_locator(minorLocator)
    ax2.tick_params(which='both', width=1)
    ax2.tick_params(which='major', length=7)
    ax2.tick_params(which='minor', length=5)
    plt.ylim(-0.01, 0.01)
    # plt.ylim(-1.02, 1.02) ## to display points for fscl=0 (log6file), just a test
    if nlin > 0:
        handles = [amean, amwgt, amedian]
    else:
        handles = [amean, amedian]
    labels = [h.get_label() for h in handles]
    legend = ax2.legend(handles, labels, title=r'No. of entries/bin = {}'.format(
        ndat), loc='upper left', numpoints=1, handlelength=0)
    legend.get_title().set_fontsize(limfs)

    if args['filename'] is None:
        fignam = proc+'.'+jobn+'.'+obsv+'.'+'scl' + \
            str(fscl)+'.'+'asymm_'+str(ndat)+'.png'
    else:
        # should fscl necessarily be added, too?
        fignam = args['filename']+'.scl' + \
            str(fscl)+'.asymm_'+str(ndat)+'.png'
    plt.savefig(fignam)
    #print "Asymmetry plot saved as: %s" %fignam
    # plt.show()

    # produces 1 plot for asymmetry and one for ratio (so it is unnecessary to return all intern variables)
    return fig1, fig2
##############################################################################
# End of function for Statistics and Plotting
##############################################################################

# The following method was not used in default settings of original script

################################################################################
# Define method for binwise distribution plots
################################################################################
def bin_dists(q, label, xlabel, col1, col2, indicators=False, lval=0., rval=0.):
    _min, _max = np.min(q), np.max(q)
    _range = _max - _min
    _nbins = min(int(ndat/5), 50)
    _bins = np.linspace(_min, _max, _nbins + 1)

    fig_dist = plt.figure(figsize=(12, 12))

    _h1, _c1, _ = plt.hist(q, bins=_bins, color=col1)
    _h2, _c2, _ = plt.hist(q, weights=wgts/np.sum(wgts)
                           * ndat, bins=_bins, color=col2)

    plt.gcf().clf()
    plt.bar(_bins[:-1], _h1, color=col1, width=0.5*_range/_nbins)
    plt.bar(_bins[:-1]+0.5*_range/_nbins, _h2,
            color=col2, width=0.5*_range/_nbins)

    plt.title(r'{}: {} {} for {}'.format(label, proc, order, obsv),
              fontsize='x-large', fontweight='bold', loc='left')
    plt.xlabel(xlabel, horizontalalignment='right',
               x=1.0, verticalalignment='top', y=1.0)
    plt.ylabel('Frequency', horizontalalignment='right',
               x=1.0, verticalalignment='bottom', y=1.0)
    raw = mpatches.Patch(color=col1, label='raw')
    wgt = mpatches.Patch(color=col2, label='weighted')
    if indicators:
        vl1 = plt.axvline(x=lval, linestyle='--', linewidth=2.0,
                          color='red', label=r'$\pm10$ IQD/2')
        vl2 = plt.axvline(x=rval, linestyle='--', linewidth=2.0, color='red')
        handles = [vl1, raw, wgt]
    else:
        handles = [raw, wgt]
    labels = [h.get_label() for h in handles]
#    plt.fill_between([-0.5,33.5],0.999,1.001, color='black', alpha=0.1)
#    plt.text(33.6, 1.00085, u'+1‰')
#    plt.text(33.7, 0.99885, u'–1‰')

    plt.legend(handles, labels, title=r'Observable bin no. {}'.format(
        i_bin), loc='upper left', numpoints=1)

    plt.yscale('log')
    plt.ylim((0.1, None))

    fig_dist_name = "{}.{}.{}.approx_{:d}_{}_bin{:d}.png".format(
        proc, jobn, obsv, ndat, label, i_bin)
    fig_dist.savefig(fig_dist_name)
    plt.close(fig_dist)

###################################################################################
# End of method for binwise distribution plots
###################################################################################

if __name__ == "__main__":
    main()
