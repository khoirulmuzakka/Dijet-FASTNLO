#!/usr/bin/env python3
#-*- coding:utf-8 -*-
#
########################################################################
#
# Plot the closure of APPLfast versus NNLOJET
#
# Created by K. Rabbertz, 22.03.2018
# Modified by B. Schillinger, 10.07.2018
# Modified by K. Rabbertz, 28.02.2020
# Prepared for python3 by K. Rabbertz, 28.02.2020
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
# fastNLO for direct evaluation of interpolation grids
# Not mandatory here since using evalution in log files
# TODO: Currently installed only for Python 2!
#import fastnlo
#from fastnlo import fastNLOLHAPDF
#from fastnlo import SetGlobalVerbosity

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
          'ytick.labelsize': 'x-large',
          'lines.linewidth': 2,
          #          'lines.markeredgewidth': 2,
          'lines.markersize': 10}
mpl.rcParams.update(params)

# Action class to allow comma-separated (or empty) list in options
class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values:
            setattr(namespace, self.dest, values[0].split(','))
        else:
            setattr(namespace, self.dest, [''])


# Some global definitions
_formats = {'eps': 0, 'pdf': 1, 'png': 2, 'svg': 3}

# Define arguments & options
parser = argparse.ArgumentParser(epilog='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# Positional arguments
parser.add_argument('table', nargs='?', type=argparse.FileType('r'),
                    help='Filename of fastNLO table to be evaluated. By default the same basename is assumed for the NNLOJET dat file and the fastNLO log file; figures are stored as tabfilebasename.sclx.scaleclosure-x.xx.png')
# Optional arguments
parser.add_argument('-d', '--datfile', required=False, nargs=1, type=argparse.FileType('r'),
                    help='NNLOJET dat file for evaluation.')
parser.add_argument('-l', '--logfile', required=False, nargs=1, type=argparse.FileType('r'),
                    help='fastNLO log file for evaluation. A direct evaluation of the table file within this code is possible, but not yet implemented.')
parser.add_argument('-f', '--filename', default=None, type=str,
                    help='Replace grid file basename by string in output figure name (optional).')
parser.add_argument('--format', required=False, nargs=1, type=str, action=SplitArgs,
                    help='Comma-separated list of plot formats to use: eps, pdf, png, svg. If nothing is chosen, png is used.')
parser.add_argument('--scalename', default=None, type=str,
                    help='Replace default scale name by given string.')
parser.add_argument('--title', default=None, type=str,
                    help='Replace table name as default title by given string.')
parser.add_argument('-v', '--verbose', action="store_true",
                    help="Increase output verbosity.")
parser.add_argument('--xlabel', default=None, type=str,
                    help='Replace x axis default label by given string.')
parser.add_argument('--ylabel', default=None, type=str,
                    help='Replace y axis default label by given string.')
parser.add_argument('--ylim2', default=0.01, type=float,
                    help='y axis limits for lower (ratio) panel.')

# Print header
print("\n###########################################################################################")
print("# fastnnlo_scaleclosure:")
print("# Plot APPLfast closure versus NNLOJET for all scale choices")
print("###########################################################################################\n")

# Parse arguments
args = vars(parser.parse_args())

# Table name if given (nargs='?' --> not a list object)
tabfile = None
if args['table']:
    tabfile = args['table'].name
if tabfile:
    print('[fastnnlo_scaleclosure]: fastNLO table: {:s}'.format(tabfile))
    # Eliminate extensions
    tabfile = re.sub('.gz$','',tabfile)
    tabfile = re.sub('.tab$','',tabfile)

# Take care of the input files
datfile = None
logfile = None
if args['datfile']:
    datfile = args['datfile'][0].name
if args['logfile']:
    logfile = args['logfile'][0].name
if not datfile:
    if tabfile:
        datfile = tabfile + '.dat'
    else:
        print('[fastnnlo_scaleclosure]: ERROR! No dat file given, neither directly nor via tablefile. Aborted!')
        sys.exit(1)
if not logfile:
    if tabfile:
        logfile = tabfile + '.log'
    else:
        print('[fastnnlo_scaleclosure]: ERROR! No log file given, neither directly nor via tablefile. Aborted!')
        sys.exit(2)
print('[fastnnlo_scaleclosure]: NNLOJET dat file: {:s}'.format(datfile))
print('[fastnnlo_scaleclosure]: fastNLO log file: {:s}'.format(logfile))

# Scale name
nice_scalename = args['scalename']

# Plot formats to use
formats = args['format']
if formats is None:
    formats = ['png']
for fmt in formats:
    if fmt not in _formats:
        print('[fastnnlo_scaleclosure]: Illegal format specified, aborted!')
        print('[fastnnlo_scaleclosure]: Format list:', args['format'])
        exit(1)
    elif fmt != 'png' and not usecairo:
        print('[fastnnlo_scaleclosure]: Vector format plots not possible without cairo backend, aborted!')
        print('[fastnnlo_scaleclosure]: Format list:', args['format'])
        exit(1)

# Plot labelling
nice_title = args['title']
nice_xlabel = args['xlabel']
nice_ylabel = args['ylabel']
ylim2 = args['ylim2']

# Verbosity
verb = args['verbose']
if verb:
    print('[fastnnlo_scaleclosure]: Using matplotlib version {:s}'.format(mpl.__version__))
    print('                       from location {:s}'.format(mpl.__file__))

# Decrypt arguments from log filename
log0base = os.path.basename(logfile)
log0args = log0base.split(".")
# Some default values: No seed index (merged grids); show reasonably large stat. uncertainty
seed  = ''
dstat = 1
# No. of scale settings to investigate (central + scale factor or fixed-scale variations)
nscl = 1

# Decrypt log filename structure
nlargs = len(log0args)
if nlargs == 4:
    proc, jobn, obsv, ext = log0args
elif nlargs == 5:
    proc, jobn, kinn, obsv, ext = log0args
elif nlargs == 6:
    proc, jobn, kinn, obsv, seed, ext = log0args
    # Do not show huge stat. uncertainty
    dstat = 0
else:
    print('[fastnnlo_scaleclosure]: ERROR! Unknown log filename structure. Aborted!')
    print('[fastnnlo_scaleclosure]: Found {:d} instead of 4-6 point-separated substrings in {:s}.'.format(nlargs,log0base))
    sys.exit(4)

print('[fastnnlo_scaleclosure]: NNLOJET process: {:s}'.format(proc))
print('[fastnnlo_scaleclosure]: NNLOJET job name: {:s}'.format(jobn))
print('[fastnnlo_scaleclosure]: Observable: {:s}'.format(obsv))
if nlargs>4:
    print('[fastnnlo_scaleclosure]: NNLOJET kinematics acronym: {:s}'.format(kinn))
if nlargs>5:
    print('[fastnnlo_scaleclosure]: NNLOJET seed index: {:s}'.format(seed))

# Extract order/contribution from job name (substring before first '-')
# LO or single contributions like R, V, RRa, RRb, RV, and VV always appear in column 6
order = jobn.split('-')[0]
ordcol = 6
if order == 'NLO':
    ordcol = 7
elif order == 'NNLO':
    ordcol = 8
print('[fastnnlo_scaleclosure]: Contribution|order: {:s}'.format(order))

# Prepare result arrays
xl = []      # left bin border
xm = []      # bin "center"
xu = []      # right bin border
xs_nnlo = []  # NNLOJET results
xs_fnlt = []  # fastNLO results
xs_fnll = []  # Pre-evaluated fastNLO results (log file)

# Read binning and cross sections from NNLOJET dat file
# Merge multiple occurences of '.' in filename by one '.'
#datfile = re.sub('\.+','.', datfile)
if verb: print('[fastnnlo_scaleclosure]: Reading from NNLOJET datfile: {:s}'.format(datfile))
xs_all = np.loadtxt(datfile, usecols=list(range(0, 17)))
xl = xs_all[:, 0]
xm = xs_all[:, 1]
xu = xs_all[:, 2]
xs_nnlo = []
dxs_nnlo = []

# Determine no. of observable bins
nobs = xl.size
if verb: print('[fastnnlo_scaleclosure]: Number of observable bins: {:d}'.format(nobs))
buf = 0
bof = nobs+1
xb = np.arange(1, nobs+1.e-6)

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
#    if ntab==nmax: break
#print 'Using ', ntab, 'table files.'
#xs_fnla = np.array(xs_fnla)

# Evaluate cross sections from pre-evaluated fastNLO tables
scalename = r'$\bf \mu$\_GeV'
for file in [logfile]:
    if verb: print('[fastnnlo_scaleclosure]: Reading from fastNLO logfile: {:s}'.format(file))
    # Skip all lines starting with "#", "C", or "L" as first non-whitespace character
    with open(file, 'r') as f:
        # C for CT, L for LHAPDF, N for NNPDF
        # TODO: Find a better solution!
        # Empty string
        data = ''
        # Read everything line-by-line
        lines = f.readlines()
        # Counter for scale variations
        nvar = 0
        sclfac = []
        murvar = []
        mufvar = []
        for line in lines:
            # Find header line(s) matching '#IObs'
            if re.search(r'#IObs', line):
                nvar = nvar+1
                # Extract scale name from within <> in this line
                scl = re.search('<(.*)>',line)
                if scl:
                    scalename = scl.group(1)
                    # Eliminate possible [] around units to avoid problems with filenames
                    scalename = re.sub(r'[\[\]]','',scalename)
                    if verb: print('[fastnnlo_scaleclosure]: Detected scale definition: {:s}'.format(scalename))
            # Scale factors
            elif re.search(r'xmur, xmuf', line):
                lend = re.sub('# The scale factors xmur, xmuf chosen here are:','',line)
                mus = lend.split(',')
                sclfac.append(True)
                murvar.append(float(mus[0]))
                mufvar.append(float(mus[1]))
            # Fixed scale settings
            elif re.search(r'mur, muf', line):
                lend = re.sub('# The fixed scales mur, muf chosen here are:','',line)
                mus = lend.split(',')
                sclfac.append(False)
                murvar.append(float(mus[0]))
                mufvar.append(float(mus[1]))
            # Skip all other lines except the ones with cross sections
            elif not re.search(r'\s*[#CLN].*', line):
                # Accumulate all numbers into one string
                data = data+line
        all_contrib = np.genfromtxt(StringIO(data), usecols=(ordcol,))
        # Calculate number of scale variations (default nscl=1, see above)
        nscl = len(all_contrib)//nobs
        # This must be identical to nvar!
        if nscl != nvar:
            print('[fastnnlo_scaleclosure]: ERROR! Found inconsistent number for scale variations in log file. Aborted!')
            print('                       nvar = {:d} while nscl = {:d}'.format(nvar,nscl))
            sys.exit(5)
        if verb: print('[fastnnlo_scaleclosure]: Number of scale variations: {:d}'.format(nscl))
        # Read the nscl*nobs values into nscl arrays of nobs entries
        for i in range(nscl):
            a = []
            for j in range(nobs):
                ind = i*nobs+j
                a.append(all_contrib[ind])
            xs_fnll.append(a)

for i in range(nscl):
    xs_nnlo.append(xs_all[:, 2*i+3]/1000)  # Conversion of fb to pb
    dxs_nnlo.append(xs_all[:, 2*i+4]/1000)

xs_fnll = np.array(xs_fnll)

# Ratio and asymmetry+1
r_fl2nn = np.divide(xs_fnll, xs_nnlo, out=np.ones_like(
    xs_fnll), where=xs_nnlo != 0)
a_fl2nn = np.divide(xs_fnll-xs_nnlo, xs_fnll+xs_nnlo,
                    out=np.zeros_like(xs_fnll), where=xs_nnlo != 0) + 1.

# Rel. stat. uncertainty from NNLOJET
dst_nnlo = np.divide(dxs_nnlo, xs_nnlo, out=np.ones_like(dxs_nnlo), where=xs_nnlo != 0)
dr_fl2nn = dstat*np.multiply(r_fl2nn, dst_nnlo)
da_fl2nn = dstat*np.multiply(a_fl2nn, dst_nnlo)

# Prepare plotting
titwgt = 'bold'
limfs = 'x-large'
if nice_scalename: scalename = nice_scalename
sclnam = [scalename]
for i in range(nscl):
    if sclfac[i]:
        sclnam.append(r'$\bf(\mu_r/\mu_0,\mu_f/\mu_0) = $ ({:4.1f},{:4.1f})'.format(murvar[i],mufvar[i]))
    else:
        sclnam.append(r'$\bf(\mu_r,\mu_f) = $ ({:4.1f},{:4.1f})'.format(murvar[i],mufvar[i]))
xmin = xl[0]
xmax = xu[nobs-1]


for i in range(nscl):

    # Closure plot
    fig = plt.figure()
    ax = fig.gca()

    if nice_title:
        title = nice_title
    elif seed == '' or seed == '_':
        title = 'Merged grid:'
    else:
        title = 'Single grid:'
    plt.title(r'{} {} {} {} for scale choice {}'.format(
        title, proc, jobn, obsv, sclnam[i]), fontweight=titwgt, y=1.02)
    if nice_xlabel:
        xlabel = nice_xlabel
    else:
        xlabel = 'Observable bin index'
    plt.xlabel(xlabel, horizontalalignment='right',
               x=1.0, verticalalignment='top', y=1.0)
    if nice_ylabel:
        ylabel = nice_ylabel
    else:
        ylabel = 'Closure quality'
    plt.ylabel(ylabel, horizontalalignment='right',
               x=1.0, verticalalignment='top', y=1.0, labelpad=30)
    plt.axhline(y=1.001, linestyle='--', linewidth=1.0, color='black')
    plt.axhline(y=0.999, linestyle='--', linewidth=1.0, color='black')
    plt.fill_between([buf, bof], 0.999, 1.001, color='black', alpha=0.1)
    if ylim2 > 0.1:
        plt.text(bof+0.7, 0.99900, u'±1‰', fontsize=limfs)
    else:
        plt.text(bof+0.6, 1.00085, u'+1‰', fontsize=limfs)
        plt.text(bof+0.7, 0.99885, u'–1‰', fontsize=limfs)

    dx = 1./8.
    x = np.arange(1, nobs+1.e-6)
    xa = np.arange(1-dx, nobs-dx+1.e-6)
    xb = np.arange(1+dx, nobs+dx+1.e-6)

#    asymm = plt.errorbar(x, a_fl2nn[i,], yerr=0.*x, marker='s', linestyle='none', label=r'Asymmetry (fastNLO-NNLOJET)/(fastNLO+NNLOJET) + 1', color='blue')
    asymm = plt.errorbar(xa, a_fl2nn[i, ], yerr=da_fl2nn[i, ], marker='s', linestyle='none',
                         label=r'Asymmetry (APPLfast-NNLOJET)/(APPLfast+NNLOJET) + 1', color='blue')
#    ratio = plt.errorbar(x, r_fl2nn[i,], yerr=0.*x, marker='o', linestyle='none', label=r'Ratio fastNLO/NNLOJET', color='orange')
    ratio = plt.errorbar(xb, r_fl2nn[i, ], yerr=dr_fl2nn[i, ], marker='o',
                         linestyle='none', label=r'Ratio APPLfast/NNLOJET', color='orange')

    plt.xlim(0.0, nobs+1)
    plt.ylim(1.-ylim2, 1.+ylim2)

    handles = [asymm, ratio]
#    handles = [ratio]
    labels = [h.get_label() for h in handles]

    legend = ax.legend(handles, labels, title=r'Closure APPLfast vs. NNLOJET',
                       loc='upper left', numpoints=1, frameon=True)
    legend.get_title().set_fontsize(limfs)
    legend._legend_box.align = 'left'

    if args['filename'] is None:
        fignam = proc+'.'+jobn+'.'+obsv+'.scl'+str(i+1)+'.scaleclosure-'+str(ylim2)+'.png'
    else:
        fignam = args['filename']+'.scl'+str(i+1)+'.scaleclosure-'+str(ylim2)+'.png'
    plt.savefig(fignam)

exit(0)
