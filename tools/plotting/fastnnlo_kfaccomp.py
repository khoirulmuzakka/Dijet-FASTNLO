#!/usr/bin/env python3
#-*- coding:utf-8 -*-
#
########################################################################
#
# Compare k factors of APPLfast versus NNLOJET incl. stat. uncertainties
#
# Created by K. Rabbertz, 18.04.2018
# Modified by B. Schillinger, 10.07.2018
# Modified by K. Rabbertz, 29.02.2020
# Prepared for python3 by K. Rabbertz, 29.02.2020
#
########################################################################
#
# python2 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# TODO: Adapt to change of genfromtxt in Python 3?
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
#print('[fastnnlo_kfaccomp]: Non GUI backends are: ', ngbackends)
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
#        print('[fastnnlo_kfaccomp]: Can not use cairo backend :-(')
#        print('                   cairocffi or pycairo are required to be installed')
    else:
        if cairo.version_info < (1, 11, 0):
            # Introduced create_for_data for Py3.
            usecairo = False
#            print('[fastnnlo_kfaccomp]: Can not use cairo backend :-(')
#            print('                   cairo {} is installed; cairo>=1.11.0 is required'.format(cairo.version))
if usecairo:
    mpl.use('cairo')
else:
    backend = 'agg'
    useagg = True
    try:
        mpl.use(backend, force=True)
        print('[fastnnlo_kfaccomp]: Warning! Could not import cairo backend :-( Using agg instead for raster plots only!')
    except:
        useagg = False
        print('[fastnnlo_kfaccomp]: Can not use agg backend :-(')
        raise ImportError('[fastnnlo_kfaccomp]: Neither cairo nor agg backend found :-( Cannot produce any plots. Good bye!')
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

# Optionally set font to Computer Modern to avoid common missing font errors
#mpl.rc('font', family='serif', serif='cm10')
#mpl.rc('text', usetex=True)
#mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

params = {'font.size': 16,
          'legend.fontsize': 'x-large',
          'figure.figsize': (16, 12),
          'mathtext.fontset': "stix",
          'axes.labelsize':  'x-large',
          'axes.titlesize':  'x-large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large',
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
#_order_color = {'LO': 'g', 'NLO': 'b', 'NNLO': 'r'}
_order_color = {'LO': 'olivedrab', 'NLO': 'dodgerblue', 'NNLO': 'orangered'}
_order_symbol = {'LO': ',', 'NLO': 's', 'NNLO': 'o'}

# Define arguments & options
parser = argparse.ArgumentParser(epilog='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# Positional arguments
parser.add_argument('datfile', nargs=1, type=argparse.FileType('r'),
                    help='Filename of highest available order NNLOJET result (.dat extension) to be evaluated. By default figures are stored as datfilebasename.kfac-index.png')
# Optional arguments
parser.add_argument('-f', '--filename', default=None, type=str,
                    help='Replace datfile basename by string in output figure name (optional).')
parser.add_argument('--format', required=False, nargs=1, type=str, action=SplitArgs,
                    help='Comma-separated list of plot formats to use: eps, pdf, png, svg. If nothing is chosen, png is used.')
parser.add_argument('-l', '--logfile', required=False, nargs='?', type=str, default='',
                    help='By default only NNLOJET results are shown. If logfile==on, the same basename is assumed for the log file (.log extension) to also show fastNLO k factors. Alternatively, an explicit filename can be specified for the log file to compare.')
parser.add_argument('--title', default=None, type=str,
                    help='Replace basefilename as default title by given string.')
parser.add_argument('-v', '--verbose', action="store_true",
                    help="Increase output verbosity.")
parser.add_argument('--xlabel', default=None, type=str,
                    help='Replace x axis default label by given string.')
parser.add_argument('--ylabel', default=None, type=str,
                    help='Replace y axis default label by given string.')

# Print header
print("\n###########################################################################################")
print("# fastnnlo_kfaccomp:")
print("# Plot k factors of APPLfast versus NNLOJET incl. stat. uncertainties")
print("###########################################################################################\n")

# Parse arguments
args = vars(parser.parse_args())

# Dat file
datfile = args['datfile'][0].name # nargs=1 gives list with exactly one entry
datbase = os.path.basename(datfile)
print('[fastnnlo_kfaccomp]: Analysing datfile: {:s}'.format(datbase))
datfiles = {}
datbases = {}
max_order = datbase.split('.')[1] # Must follow naming convention: proc.order.rest.dat
if max_order == 'LO':
    print('[fastnnlo_kfaccomp]: ERROR! k factors always need > LO to be known. Aborted!')
    print('[fastnnlo_kfaccomp]: Found {:s} instead as highest order.'.format(max_order))
    sys.exit(2)
elif max_order != 'NLO' and max_order != 'NNLO':
    print('[fastnnlo_kfaccomp]: ERROR! Unknown higher order {:s}. Aborted!'.format(max_order))
    sys.exit(3)
else:
    if max_order == 'NLO':
        datfiles['NLO'] = datfile
        datbases['NLO'] = datbase
        lofile = re.sub('.NLO.','.LO.',datfile)
    else:
        datfiles['NNLO'] = datfile
        datbases['NNLO'] = datbase
        nlofile = re.sub('.NNLO.','.NLO.',datfile)
        if not os.path.isfile(nlofile):
            print('[fastnnlo_kfaccomp]: ERROR! Could not find NLO NNLOJET dat file. Aborted!')
            print('[fastnnlo_kfaccomp]: Filename is: {:s}'.format(nlofile))
            sys.exit(4)
        print('[fastnnlo_kfaccomp]: Found NLO file: {}'.format(nlofile))
        datfiles['NLO'] = nlofile
        datbases['NLO'] = os.path.basename(nlofile)
        lofile = re.sub('.NNLO.','.LO.',datfile)
    if not os.path.isfile(lofile):
        print('[fastnnlo_kfaccomp]: ERROR! Could not find LO NNLOJET dat file. Aborted!')
        print('[fastnnlo_kfaccomp]: Filename is: {:s}'.format(lofile))
        sys.exit(5)
    print('[fastnnlo_kfaccomp]: Found LO file: {}'.format(lofile))
    datfiles['LO'] = lofile
    datbases['LO'] = os.path.basename(lofile)

# Given filename
given_filename = args['filename']

# Log file with all necessary fastNLO cross sections
logfile = args['logfile']
print (logfile)
if logfile == 'on':
    logfile = re.sub('dat$', 'log', datfiles[max_order])
if logfile:
    if not os.path.isfile(logfile):
        print('[fastnnlo_kfaccomp]: ERROR! Could not find fastNLO log file. Aborted!')
        print('[fastnnlo_kfaccomp]: Filename is: {:s}'.format(logfile))
        sys.exit(6)
    logbase = os.path.basename(logfile)
    print('[fastnnlo_kfaccomp]: Found logfile: {}'.format(logfile))
else:
    logbase = ''
    print('[fastnnlo_kfaccomp]: Only NNLOJET results are shown.')

# Plot formats to use
formats = args['format']
if formats is None:
    formats = ['png']
for fmt in formats:
    if fmt not in _formats:
        print('[fastnnlo_kfaccomp]: Illegal format specified, aborted!')
        print('[fastnnlo_kfaccomp]: Format list:', args['format'])
        exit(1)
    elif fmt != 'png' and not usecairo:
        print('[fastnnlo_kfaccomp]: Vector format plots not possible without cairo backend, aborted!')
        print('[fastnnlo_kfaccomp]: Format list:', args['format'])
        exit(1)

# Plot labelling
nice_title = args['title']
nice_xlabel = args['xlabel']
nice_ylabel = args['ylabel']

# Verbosity
verb = args['verbose']
if verb:
    print('[fastnnlo_kfaccomp]: Using matplotlib version {:s}'.format(mpl.__version__))
    print('                     from location {:s}'.format(mpl.__file__))

# Decrypt arguments from dat filename
datargs = datbase.split(".") # Array containing filename-parts
# Decrypt dat filename structure
ndargs = len(datargs)
if ndargs == 4:
    proc, jobn, obsv, ext = datargs
elif ndargs == 5:
    proc, jobn, kinn, obsv, ext = datargs
elif ndargs == 6:
    proc, jobn, kinn, obsv, seed, ext = datargs
    # Do not show huge stat. uncertainty
    dstat = 0
else:
    print('[fastnnlo_kfaccomp]: ERROR! Unknown filename structure. Aborted!')
    print('[fastnnlo_kfaccomp]: Found {:d} instead of 4-6 point-separated substrings in {:s}.'.format(ndargs,datbase))
    sys.exit(2)

print('[fastnnlo_kfaccomp]: NNLOJET process: {:s}'.format(proc))
print('[fastnnlo_kfaccomp]: NNLOJET order name: {:s}'.format(jobn))
if ndargs>4:
    print('[fastnnlo_kfaccomp]: NNLOJET kinematics acronym: {:s}'.format(kinn))
    print('[fastnnlo_kfaccomp]: ERROR! Only full cross sections (no kinematics acronym) allowed. Aborted!')
    sys.exit(3)
if ndargs>5:
    print('[fastnnlo_kfaccomp]: NNLOJET seed index: {:s}'.format(seed))
    print('[fastnnlo_kfaccomp]: ERROR! Only full cross sections (no seed index) allowed. Aborted!')
    sys.exit(4)
print('[fastnnlo_kfaccomp]: Observable: {:s}'.format(obsv))

# some default values
xaxe = 'bins'  # x axis with bin numbers ('bins') or physics observable

# Prepare result arrays
xl = []  # left bin border
xm = []  # bin "center"
xu = []  # right bin border

# Read binning and cross sections from NNLOJET dat files
# Loop over orders
for order in datfiles.keys():
    datfile = datfiles[order]
    print('[fastnnlo_kfaccomp]: Reading from NNLOJET dat file {:s}'.format(datfile))
    xs_all = np.loadtxt(datfile, usecols=range(0, 17))
    if order == 'LO':
        xl = xs_all[:, 0]  # Binning must be identical for all orders
        xm = xs_all[:, 1]
        xu = xs_all[:, 2]
        xs_lo = xs_all[:, 3]/1000.  # Conversion of fb to pb
        dxs_lo = xs_all[:, 4]/1000.
    elif order == 'NLO':
        xs_nlo = xs_all[:, 3]/1000.
        dxs_nlo = xs_all[:, 4]/1000.
    else:
        xs_nnlo = xs_all[:, 3]/1000.
        dxs_nnlo = xs_all[:, 4]/1000.

xs_lo = np.array(xs_lo)
dxs_lo = np.array(dxs_lo)
xs_nlo = np.array(xs_nlo)
dxs_nlo = np.array(dxs_nlo)
if max_order == 'NNLO':
    xs_nnlo = np.array(xs_nnlo)
    dxs_nnlo = np.array(dxs_nnlo)

# Determine no. of observable bins
nobs = xl.size
print('[fastnnlo_kfaccomp]: Number of observable bins: {:d}'.format(nobs))
dx  = 1./8.
xb  = np.arange(1, nobs+1.e-6)
xsa = np.arange(1-dx, nobs-dx+1.e-6)
xsb = np.arange(1+dx, nobs+dx+1.e-6)

# Read cross sections from pre-evaluated fastNLO tables
ordcol = {'LO': 6, 'NLO': 7, 'NNLO': 8}
if logfile:
    # Skip all lines starting with "#", "C", or "L" as first non-whitespace character
    # why don't we start these lines (C, L...) also with a "#" ???
    # Read first nobs values for central scale result
    print('[fastnnlo_kfaccomp]: Reading from fastNLO log file {:s}'.format(logfile))
    with open(logfile, 'r') as f:
        data = re.sub(r'\s*[#CLN].*', '', f.read())
        lo = np.genfromtxt(StringIO(data), usecols=ordcol['LO'],)
        nlo = np.genfromtxt(StringIO(data), usecols=ordcol['NLO'],)
        if max_order == 'NNLO':
            nnlo = np.genfromtxt(StringIO(data), usecols=ordcol['NNLO'],)
        xs_flo = []
        xs_fnlo = []
        if max_order == 'NNLO':
            xs_fnnlo = []
        for irow in range(nobs):
            xs_flo.append(lo[irow])
            xs_fnlo.append(nlo[irow])
            if max_order == 'NNLO':
                xs_fnnlo.append(nnlo[irow])

    xs_flo = np.array(xs_flo)
    xs_fnlo = np.array(xs_fnlo)
    if max_order == 'NNLO':
        xs_fnnlo = np.array(xs_fnnlo)

# Successive k factors
kn_nlo  = np.divide(xs_nlo, xs_lo, out=np.ones_like(xs_nlo), where=xs_lo != 0)
if logfile: kf_nlo  = np.divide(xs_fnlo, xs_flo, out=np.ones_like(xs_fnlo), where=xs_flo!=0)
if max_order == 'NNLO':
    kn_nnlo = np.divide(xs_nnlo, xs_nlo, out=np.ones_like(xs_nnlo), where=xs_nlo != 0)
    if logfile: kf_nnlo = np.divide(xs_fnnlo, xs_fnlo, out=np.ones_like(xs_fnnlo), where=xs_fnlo!=0)

# Rel. stat. uncertainty of K factors from NNLOJET
dst_lo = np.divide(dxs_lo, xs_lo, out=np.ones_like(dxs_lo), where=xs_lo != 0)
dst_nlo = np.divide(dxs_nlo, xs_nlo, out=np.ones_like(dxs_nlo), where=xs_nlo != 0)
if max_order == 'NNLO':
    dst_nnlo = np.divide(dxs_nnlo, xs_nnlo, out=np.ones_like(dxs_nnlo), where=xs_nnlo != 0)
    dst_nnlo = np.sqrt(np.multiply(dst_nnlo, dst_nnlo) + np.multiply(dst_nlo, dst_nlo))
    dk_nnlo = np.multiply(kn_nnlo, dst_nnlo)
dst_nlo = np.sqrt(np.multiply(dst_nlo, dst_nlo) + np.multiply(dst_lo, dst_lo))
dk_nlo = np.multiply(kn_nlo, dst_nlo)

# Prepare plotting
titwgt = 'bold'
limfs = 'x-large'

if xaxe == 'bins':
    x = xb
    dx = 0.5*np.ones_like(x)
    xscl = 'linear'
    xmin = 0
    xmax = nobs+1
else:
    x = xm
    dx = (xu-xl)/2.
#    xscl = 'linear'
    xscl = 'log'
    xmin = xl[0]
    xmax = xu[nobs-1]

fig = plt.figure()
gs  = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0)
ax  = plt.subplot(gs[0])
ax.set_autoscalex_on(False)
plt.setp(ax.get_xticklabels(), visible=False)

# Upper subplot setup
#
# Set limits and scaling of x axis through coupled axis from lower plot
ax.set_yscale('log')
if nice_ylabel:
    ylabel = nice_ylabel
else:
    ylabel = '$\sigma \pm \Delta\sigma_\mathrm{stat}$'
ax.set_ylabel(r'%s' % ylabel, horizontalalignment='right',
              x=1.0, verticalalignment='top', y=1.0, labelpad=20)
if nice_title:
    title = nice_title
else:
    title = re.sub('.dat$', '', datbases[max_order])
ax.set_title(r'{:s}'.format(title), loc='left', fontweight=titwgt, y=1.05)
#ax.set_title('%s' % title)

if nice_xlabel:
    xlabel = nice_xlabel
else:
    xlabel = 'Observable bin index'

# Cross sections
axhandles = []
abslo = ax.errorbar(x, xs_lo, xerr=dx, capsize=0, yerr=dxs_lo, linestyle='none',
                    label=r'{:s}'.format('LO'), marker=_order_symbol['LO'], color=_order_color['LO'])
axhandles.append(abslo)
absnlo = ax.errorbar(x, xs_nlo, xerr=dx, capsize=0, yerr=dxs_nlo, linestyle='none',
                     label=r'{:s}'.format('NLO'), marker=_order_symbol['NLO'], color=_order_color['NLO'])
axhandles.append(absnlo)
if max_order == 'NNLO':
    absnnlo = ax.errorbar(x, xs_nnlo, xerr=dx, capsize=0, yerr=dxs_nnlo, linestyle='none',
                          label=r'{:s}'.format('NNLO'), marker=_order_symbol['NNLO'], color=_order_color['NNLO'],
                          mfc='none', mec=_order_color['NNLO'], mew=2)
    axhandles.append(absnnlo)

# Lower subplot setup
#
axr = plt.subplot(gs[1], sharex=ax)
# Set common x axis for both subplots
axr.set_xlim(xmin, xmax)
axr.set_xscale(xscl)
axr.set_yscale('linear')
axr.set_xlabel(r'{:s}'.format(xlabel), horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
axr.set_ylabel(r'K factors', horizontalalignment='center', verticalalignment='center', labelpad=24)
#axr.fill_between(x, 1.-abs(dst_nnlo[i]), 1.+abs(dst_nnlo[i]), edgecolor=col1, facecolor=col1b, alpha=0.5)
axr.axhline(1.0, color='black', ls='--')

if logfile: x = xsa
rkn = axr.errorbar(x, kn_nlo, yerr=dk_nlo, marker='s', linestyle='none', label='', color=_order_color['NLO'])
if logfile: rkf = axr.errorbar(xsb, kf_nlo, yerr=dk_nlo, marker='s', linestyle='none', label='APPLfast NLO', color='cyan', mfc='none', mec='cyan', mew=2)
axrhandles = []
if logfile:
    axrhandles.append(rkf)
if max_order == 'NNLO':
    rknn = axr.errorbar(x, kn_nnlo, yerr=dk_nnlo, marker='o', linestyle='none', label='', color=_order_color['NNLO'])
    if logfile: rknf = axr.errorbar(xsb, kf_nnlo, yerr=dk_nnlo, marker='o', linestyle='none', label='APPLfast NNLO', color='orange', mfc='none', mec='orange', mew=2)
    if logfile:
        axrhandles.append(rknf)

axlabels = [h.get_label() for h in axhandles]
axrlabels = [h.get_label() for h in axrhandles]

legend = ax.legend(numpoints=1)
if logfile: legend = axr.legend(axrhandles, axrlabels, fontsize=12, numpoints=1)

if xaxe == 'bins':
    if args['filename'] is None:
        fignam = proc+'.'+jobn+'.'+obsv+'.'+'kfac-index'+'.png'
    else:
        fignam = args['filename']+'.kfac-index.png'
else:
    if args['filename'] is None:
        fignam = proc+'.'+jobn+'.'+obsv+'.'+'kfac-obs'+'.png'
    else:
        fignam = args['filename']+'.kfac-obs.png'

print('[fastnnlo_kfaccomp]: Writing figure {:s}'.format(fignam))

plt.savefig(fignam)

exit(0)
