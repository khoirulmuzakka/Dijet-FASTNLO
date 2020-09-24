#!/usr/bin/env python3
#-*- coding:utf-8 -*-
#
########################################################################
#
# Plot absolute cross sections of APPLfast versus NNLOJET incl. stat. uncertainties
#
# Created by K. Rabbertz, 18.04.2018
# Modified by B. Schillinger, 10.07.2018
# Modified by K. Rabbertz, 31.10.2019
# Prepared for python3 by K. Rabbertz, 09.03.2020
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
#print('[fastnnlo_absolute]: Non GUI backends are: ', ngbackends)
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
#        print('[fastnnlo_absolute]: Can not use cairo backend :-(')
#        print('                   cairocffi or pycairo are required to be installed')
    else:
        if cairo.version_info < (1, 11, 0):
            # Introduced create_for_data for Py3.
            usecairo = False
#            print('[fastnnlo_absolute]: Can not use cairo backend :-(')
#            print('                   cairo {} is installed; cairo>=1.11.0 is required'.format(cairo.version))
if usecairo:
    mpl.use('cairo')
else:
    backend = 'agg'
    useagg = True
    try:
        mpl.use(backend, force=True)
        print('[fastnnlo_absolute]: Warning! Could not import cairo backend :-( Using agg instead for raster plots only!')
    except:
        useagg = False
        print('[fastnnlo_absolute]: Can not use agg backend :-(')
        raise ImportError('[fastnnlo_absolute]: Neither cairo nor agg backend found :-( Cannot produce any plots. Good bye!')
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
_fntrans = str.maketrans({'[': '', ']': '', '(': '', ')': '', ',': ''}) # Filename translation table
_formats = {'eps': 0, 'pdf': 1, 'png': 2, 'svg': 3}
#_order_color = {'LO': 'g', 'NLO': 'b', 'NNLO': 'r'}
_order_color = {'LO': 'olivedrab', 'NLO': 'dodgerblue', 'NNLO': 'orangered'}
_order_symbol = {'LO': ',', 'NLO': 's', 'NNLO': 'o'}

# Define arguments & options
parser = argparse.ArgumentParser(epilog='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# Positional arguments
parser.add_argument('datfile', nargs=1, type=argparse.FileType('r'),
                    help='Filename of NNLOJET result (.dat extension) to be evaluated. By default the same basename is assumed for the fastNLO log file; figures are stored as datfilebasename.sclx.absolute-index-x.xx.png')
# Optional arguments
parser.add_argument('-f', '--filename', default=None, type=str,
                    help='Replace datfile basename by string in output figure name (optional).')
parser.add_argument('--format', required=False, nargs=1, type=str, action=SplitArgs,
                    help='Comma-separated list of plot formats to use: eps, pdf, png, svg. If nothing is chosen, png is used.')
parser.add_argument('-l', '--logfile', required=False, nargs=1, type=argparse.FileType('r'),
                    help='fastNLO log file for evaluation. A direct evaluation of the table file within this code is possible, but not yet implemented.')
parser.add_argument('--scalename', default=None, type=str,
                    help='Replace default scale name by given string.')
parser.add_argument('--title', default=None, type=str,
                    help='Replace basefilename as default title by given string.')
parser.add_argument('-v', '--verbose', action="store_true",
                    help="Increase output verbosity.")
parser.add_argument('--xlabel', default=None, type=str,
                    help='Replace x axis default label by given string.')
parser.add_argument('--ylabel', default=None, type=str,
                    help='Replace y axis default label by given string.')
parser.add_argument('--ylim2', default=None, type=float,
                    help='y axis limits for lower (ratio) panel. If not set, automatic choice by matplotlib.')

# Print header
print("\n###########################################################################################")
print("# fastnnlo_absolute:")
print("# Plot absolute cross sections of APPLfast versus NNLOJET incl. stat. uncertainties")
print("###########################################################################################\n")

# Parse arguments
args = vars(parser.parse_args())

# Dat file
datfile = args['datfile'][0].name # nargs=1 gives list with exactly one entry
datbase = os.path.basename(datfile)
print('[fastnnlo_absolute]: NNLOJET dat file: {:s}'.format(datfile))

# Take care of the input files
logfile = None
if args['logfile']:
    logfile = args['logfile'][0].name
if not logfile:
    if datfile:
        logfile = re.sub('.dat$', '.log', datfile)
    else:
        print('[fastnnlo_absolute]: ERROR! No log file given, neither directly nor via datfile. Aborted!')
        sys.exit(1)
print('[fastnnlo_absolute]: fastNLO log file: {:s}'.format(logfile))

# Scale name
nice_scalename = args['scalename']

# Plot formats to use
formats = args['format']
if formats is None:
    formats = ['png']
for fmt in formats:
    if fmt not in _formats:
        print('[fastnnlo_absolute]: Illegal format specified, aborted!')
        print('[fastnnlo_absolute]: Format list:', args['format'])
        exit(1)
    elif fmt != 'png' and not usecairo:
        print('[fastnnlo_absolute]: Vector format plots not possible without cairo backend, aborted!')
        print('[fastnnlo_absolute]: Format list:', args['format'])
        exit(1)

# Plot labelling
nice_title = args['title']
nice_xlabel = args['xlabel']
nice_ylabel = args['ylabel']
ylim2 = args['ylim2']

# Verbosity
verb = args['verbose']
if verb:
    print('[fastnnlo_absolute]: Using matplotlib version {:s}'.format(mpl.__version__))
    print('                       from location {:s}'.format(mpl.__file__))

# Decrypt arguments from log filename
log0base = os.path.basename(logfile)
log0args = log0base.split(".")
# Some default values: No seed index (merged grids); show reasonably large stat. uncertainty
seed  = ''
dstat = 1
# No. of scale settings to investigate (central + scale factor or fixed-scale variations)
nscl = 1
ylim = 0.01
xaxe = 'bins'  # x axis with bin numbers ('bins') or physics observable
# TODO shall xaxe be optional?

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
    print('[fastnnlo_absolute]: ERROR! Unknown log filename structure. Aborted!')
    print('[fastnnlo_absolute]: Found {:d} instead of 4-6 point-separated substrings in {:s}.'.format(nlargs,log0base))
    sys.exit(4)

print('[fastnnlo_absolute]: NNLOJET process: {:s}'.format(proc))
print('[fastnnlo_absolute]: NNLOJET job name: {:s}'.format(jobn))
print('[fastnnlo_absolute]: Observable: {:s}'.format(obsv))
if nlargs>4:
    print('[fastnnlo_absolute]: NNLOJET kinematics acronym: {:s}'.format(kinn))
if nlargs>5:
    print('[fastnnlo_absolute]: NNLOJET seed index: {:s}'.format(seed))

# Extract order/contribution from job name (substring before first '-')
# LO or single contributions like R, V, RRa, RRb, RV, and VV always appear in column 6
order = jobn.split('-')[0]
ordcol = 6
if order == 'NLO':
    ordcol = 7
elif order == 'NNLO':
    ordcol = 8
print('[fastnnlo_absolute]: Contribution|order: {:s}'.format(order))

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
if verb: print('[fastnnlo_absolute]: Reading from NNLOJET dat file: {:s}'.format(datfile))
xs_all = np.loadtxt(datfile, usecols=range(0, 17))
xl = xs_all[:, 0]
xm = xs_all[:, 1]
xu = xs_all[:, 2]
xs_nnlo = []
dxs_nnlo = []

# Determine no. of observable bins
nobs = xl.size
if verb: print('[fastnnlo_absolute]: Number of observable bins: {:d}'.format(nobs))
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
#for file in [logfile]:
#    if verb: print('[fastnnlo_absolute]: Reading from fastNLO logfile: {:s}'.format(file))
#    # Skip all lines starting with "#", "C", or "L" as first non-whitespace character
#    with open(file, 'r') as f:
#        data = re.sub(r'\s*[#CLN].*', '', f.read())
#        all_contrib = np.genfromtxt(StringIO(data), usecols=(ordcol,))
#        # Calculate number of scale variations (default nscl=1, see above)
#        nscl = len(all_contrib)/nobs
#        print "nscl", nscl
#        # Read the nscl*nobs values into nscl arrays of nobs entries
#        for i in range(nscl):  # consider all scale variations given in logfile
#            a = []
#            for j in range(nobs):
#                ind = i*nobs+j
#                a.append(all_contrib[ind])
#            xs_fnll.append(a)

# Evaluate cross sections from pre-evaluated fastNLO tables
scalename = r'$\bf \mu$\_GeV'
for file in [logfile]:
    if verb: print('[fastnnlo_absolute]: Reading from fastNLO logfile: {:s}'.format(file))
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
                    # Do not use characters defined in _fntrans for filenames
                    scalename = scalename.translate(_fntrans)
                    if verb: print('[fastnnlo_absolute]: Detected scale definition: {:s}'.format(scalename))
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
            print('[fastnnlo_absolute]: ERROR! Found inconsistent number for scale variations in log file. Aborted!')
            print('                       nvar = {:d} while nscl = {:d}'.format(nvar,nscl))
            sys.exit(5)
        if verb: print('[fastnnlo_absolute]: Number of scale variations: {:d}'.format(nscl))
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
r_nn2nn = np.ones_like(xs_nnlo)
r_fl2nn = np.divide(xs_fnll, xs_nnlo, out=np.ones_like(
    xs_fnll), where=xs_nnlo != 0)
a_fl2nn = np.divide(xs_fnll-xs_nnlo, xs_fnll+xs_nnlo,
                    out=np.zeros_like(xs_fnll), where=xs_nnlo != 0) + 1.

# Rel. stat. uncertainty from NNLOJET
dst_nnlo = np.divide(dxs_nnlo, xs_nnlo, out=np.ones_like(dxs_nnlo), where=xs_nnlo != 0)
dr_fl2nn = np.multiply(r_fl2nn, dst_nnlo)
da_fl2nn = np.multiply(a_fl2nn, dst_nnlo)

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


# TODO: Needed?
if xaxe == 'bins':
    x = xb
    dx = 0.5*np.ones_like(x)
    xscl = 'linear'
    xlab = 'bin index'
    xmin = 0
    xmax = nobs+1
else:
    x = xm
    dx = (xu-xl)/2.
#    xscl = 'linear'
    xscl = 'log'
    xlab = xaxe
    xmin = xl[0]
    xmax = xu[nobs-1]


for i in range(nscl):

    # Absolute predictions
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0)
    ax = plt.subplot(gs[0])
    axr = plt.subplot(gs[1], sharex=ax)

    if nice_title:
        title = nice_title
    elif seed == '' or seed == '_':
        title = 'Merged grid:'
    else:
        title = 'Single grid:'
    ax.set_title(r'{} {} {} {} for scale choice {}'.format(title, proc, jobn, obsv, sclnam[i]), fontweight=titwgt, y=1.02)
    if nice_xlabel:
        xlabel = nice_xlabel
    else:
        xlabel = 'Observable bin index'
    plt.setp(ax.get_xticklabels(), visible=False)
    axr.set_xlabel(xlabel, horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
    if nice_ylabel:
        ylabel = nice_ylabel
    else:
        ylabel = r'$\bf\left|\mathrm{d}\sigma/\mathrm{d}O\right|$ [pb/GeV]'
    ax.set_ylabel(ylabel, horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, labelpad=20)
    axr.set_ylabel('Ratio', horizontalalignment='center', verticalalignment='center', labelpad=20)
    #    ax.set_xscale(xscl)
    axr.set_xscale(xscl)
    ax.set_yscale('log')
    axr.set_yscale('linear')

    # Switch of idiotic change to exponential tick labels; TODO: Do this in general in other plots as well
    axr.get_yaxis().get_major_formatter().set_useOffset(False)

    axhandles = []
    # Booleans for + and - x sections
    xmp = xs_nnlo[i] > 0
    xmn = xs_nnlo[i] < 0
    if len(x[xmp]):
        abs1p = ax.errorbar(x[xmp], +xs_nnlo[i][xmp], xerr=dx[xmp], capsize=0, yerr=dxs_nnlo[i]
                            [xmp], marker=',', linestyle='none', label=r'NNLOJET $\bf\pm\Delta_{stat}$', color=_order_color['LO'])
        axhandles.append(abs1p)
    if len(x[xmn]):
        abs1n = ax.errorbar(x[xmn], -xs_nnlo[i][xmn], xerr=dx[xmn], capsize=0, yerr=dxs_nnlo[i]
                            [xmn], marker=',', linestyle='none', label=r'NNLOJET', color=_order_color['LO'])
        #        axhandles.append(abs1n)

    xmp = xs_fnll[i] > 0
    xmn = xs_fnll[i] < 0
    if len(x[xmp]):
        abs2p = ax.errorbar(x[xmp], +xs_fnll[i][xmp], marker='d', linestyle='none',
                            label=r'APPLfast grid ($\bf\sigma>0$)', color=_order_color['NNLO'])
        axhandles.append(abs2p)
    if len(x[xmn]):
        abs2n = ax.errorbar(x[xmn], -xs_fnll[i][xmn], marker='d', linestyle='none',
                            label=r'APPLfast grid ($\bf\sigma<0$)', color=_order_color['NLO'], markerfacecolor=_order_color['NLO'])
        axhandles.append(abs2n)

    ax.set_xlim(xmin, xmax)
    axr.set_xlim(xmin, xmax)
    if ylim2:
        axr.set_ylim(1.-ylim2, 1.+ylim2)
    axr.fill_between(x, 1.-abs(dst_nnlo[i]), 1.+abs(dst_nnlo[i]),
                     edgecolor=_order_color['LO'], facecolor=_order_color['LO'], alpha=0.5)
    axr.axhline(1.0, color=_order_color['LO'])

    rnn = axr.errorbar(x, r_nn2nn[i], marker=',',
                       linestyle='none', label='', color=_order_color['LO'])
    rfp = axr.errorbar(x[xmp], r_fl2nn[i][xmp], marker='d',
                       linestyle='none', label='', color=_order_color['NNLO'])
    rfn = axr.errorbar(x[xmn], r_fl2nn[i][xmn], marker='d',
                       linestyle='none', label='', color=_order_color['NLO'], markerfacecolor=_order_color['NLO'])
    axrhandles = [rnn, rfp, rfn]

    axlabels = [h.get_label() for h in axhandles]
    axrlabels = [h.get_label() for h in axrhandles]

    legend = ax.legend(axhandles, axlabels, # title=r'Stat. uncertainty & grid closure',
                       loc='upper right', numpoints=1, frameon=True)
    legend.get_title().set_fontsize(limfs)

    if xaxe == 'bins':
        if args['filename'] is None:
            fignam = proc+'.'+jobn+'.'+obsv+'.scl'+str(i+1)+'.absolute-index-'+str(ylim2)+'.png'
        else:
            fignam = args['filename']+'.scl'+str(i+1)+'.absolute-index-'+str(ylim2)+'.png'
    else:
        if args['filename'] is None:
            fignam = proc+'.'+jobn+'.'+obsv+'.scl'+str(i+1)+'.absolute-obs-'+str(ylim2)+'.png'
        else:
            fignam = args['filename']+'.scl'+str(i+1)+'.absolute-index-'+str(ylim2)+'.png'
    if verb: print('[fastnnlo_absolute]: Writing figure: {:s}'.format(fignam))
    plt.savefig(fignam, bbox_inches='tight')

exit(0)
