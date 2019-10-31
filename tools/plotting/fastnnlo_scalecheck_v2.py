#!/usr/bin/env python2
#-*- coding:utf-8 -*-
import argparse
import matplotlib as mpl
mpl.use('Agg')  # otherwise some matplotlib error occurs
import glob
import os
import sys
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
# from copy import deepcopy
from matplotlib import cm
from fastnlo import fastNLOLHAPDF
from fastnlo import SetGlobalVerbosity
import re
from StringIO import StringIO
#import warnings
# warnings.filterwarnings("error")

# Style settings
# Location of matplotlibrc
# print mpl.matplotlib_fname()
# List active style settings
# print mpl.rcParams
# Needs matplotlib > 1.3
# plt.style.use('ggplot')
# plt.style.use('presentation')


parser = argparse.ArgumentParser()

# add arguments (one .dat and two .log files)
parser.add_argument('-d', '--datfile', default='file.dat', required=True,
                    help='.dat file for evaluation.')
parser.add_argument('-l', '--logfile', default='file.log', required=True,
                    help='.log file, need one *.log file')

parser.add_argument('-o', '--outputfilename', required=False, nargs='?', type=str,
                    help='Customise the first part of the output filename.'
                    'Default: Same structure as datfile name.')

# parse arguments
args = vars(parser.parse_args())
namesp = parser.parse_args()

# take care of the input files
datfile = args['datfile']
logfile = args['logfile']
print logfile
print type(logfile)

params = {'legend.fontsize': 'x-large',
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

# arguments
log0base = os.path.basename(logfile)
#proc = os.path.splitext(logfile)[0]
log0args = log0base.split(".")
print log0args

# some default values
seed = ''  # default
# no. of scale settings to investigate (central + scale factor or fixed scale variations)
nscl = 1
ylim = 0.01

# arguments from logfile
if len(log0args) == 5:
    proc, jobn, kinn, obsv, ext = log0args
    print 'proc: ', proc
    print 'jobn: ', jobn
    print 'kinn: ', kinn
    print 'obsv: ', obsv
    print '\n'
elif len(log0args) == 6:
    proc, jobn, kinn, obsv, seed, ext = log0args
    print 'proc: ', proc
    print 'jobn: ', jobn
    print 'kinn: ', kinn
    print 'obsv: ', obsv
    print 'seed: ', seed
    print '\n'
elif len(log0args) == 4:
    proc, jobn, obsv, ext = log0args
    print 'proc: ', proc
    print 'jobn: ', jobn
    print 'obsv: ', obsv
    print '\n'


# Extract order/contribution from job type (substring before first '-')
order = jobn.split('-')[0]
ordcol = 6
if order == 'NLO':
    ordcol = 7
elif order == 'NNLO':
    ordcol = 8

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
print 'Reading from NNLOJET dat file ', datfile
xs_all = np.loadtxt(datfile, usecols=range(0, 17))
xl = xs_all[:, 0]
xm = xs_all[:, 1]
xu = xs_all[:, 2]
xs_nnlo = []
dxs_nnlo = []

# Determine no. of observable bins
nobs = xl.size
print 'Number of observable bins: ', nobs
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
for file in [logfile]:
    print 'Reading from fastNLO log file ', file
    # Skip all lines starting with "#", "C", or "L" as first non-whitespace character
    with open(file, 'r') as f:
        # C for CT, L for LHAPDF, N for NNPDF
        # TODO: Find a better solution!
        data = re.sub(r'\s*[#CLN].*', '', f.read())
        all_contrib = np.genfromtxt(StringIO(data), usecols=(ordcol,))
        # Calculate number of scale variations (default nscl=1, see above)
        nscl = len(all_contrib)/nobs
        print "nscl", nscl
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

r_fl2nn = np.divide(xs_fnll, xs_nnlo, out=np.ones_like(
    xs_fnll), where=xs_nnlo != 0)
a_fl2nn = np.divide(xs_fnll-xs_nnlo, xs_fnll+xs_nnlo,
                    out=np.zeros_like(xs_fnll), where=xs_nnlo != 0) + 1.

# Rel. stat. uncertainty from NNLOJET
dst_nnlo = np.divide(dxs_nnlo, xs_nnlo, out=np.ones_like(
    dxs_nnlo), where=xs_nnlo != 0)
dr_fl2nn = np.multiply(r_fl2nn, dst_nnlo)
da_fl2nn = np.multiply(a_fl2nn, dst_nnlo)

# Prepare plotting
titwgt = 'bold'
limfs = 'x-large'
sclnam = [r'$\bf p_{T,max}$', r'$\bf\mu_r = \mu_f$ small', r'$\bf\mu_r = \mu_f$ large',
          r'$\bf\mu_r < \mu_f$', r'$\bf\mu_r > \mu_f$', r'$\bf\mu_r \ll \mu_f$', r'$\bf\mu_r \gg \mu_f$']
xmin = xl[0]
xmax = xu[nobs-1]


for i in range(nscl):

    # Closure plot
    fig = plt.figure()
    ax = fig.gca()

    if seed == '' or seed == '_':
        plt.title(r'Merged grid: {} {} {} for scale choice {}'.format(
            proc, jobn, obsv, sclnam[i]), fontweight=titwgt, y=1.05)
    else:
        plt.title(r'Single grid: {} {} {} for scale choice {}'.format(
            proc, jobn, obsv, sclnam[i]), fontweight=titwgt, y=1.05)
    plt.xlabel('Observable bin index', horizontalalignment='right',
               x=1.0, verticalalignment='top', y=1.0)
#    plt.ylabel('Closure fastNLO vs. NNLOJET', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
    plt.ylabel('Closure quality', horizontalalignment='right',
               x=1.0, verticalalignment='top', y=1.0, labelpad=30)
    plt.axhline(y=1.001, linestyle='--', linewidth=1.0, color='black')
    plt.axhline(y=0.999, linestyle='--', linewidth=1.0, color='black')
    plt.fill_between([0.0, 34.0], 0.999, 1.001, color='black', alpha=0.1)
    if ylim > 0.1:
        plt.text(34.7, 0.99900, u'±1‰', fontsize=limfs)
    else:
        plt.text(34.6, 1.00085, u'+1‰', fontsize=limfs)
        plt.text(34.7, 0.99885, u'–1‰', fontsize=limfs)

    dx = 1./4.
    x = np.arange(1, nobs+1.e-6)
    xa = np.arange(1-dx, nobs-dx+1.e-6)
    xb = np.arange(1+dx, nobs+dx+1.e-6)

#    asymm = plt.errorbar(x, a_fl2nn[i,], yerr=0.*x, marker='s', linestyle='none', label=r'Asymmetry (fastNLO-NNLOJET)/(fastNLO+NNLOJET) + 1', color='blue')
    asymm = plt.errorbar(x, a_fl2nn[i, ], yerr=da_fl2nn[i, ], marker='s', linestyle='none',
                         label=r'Asymmetry (APPLfast-NNLOJET)/(APPLfast+NNLOJET) + 1', color='blue')
#    ratio = plt.errorbar(x, r_fl2nn[i,], yerr=0.*x, marker='o', linestyle='none', label=r'Ratio fastNLO/NNLOJET', color='orange')
    ratio = plt.errorbar(x, r_fl2nn[i, ], yerr=dr_fl2nn[i, ], marker='o',
                         linestyle='none', label=r'Ratio APPLfast/NNLOJET', color='orange')

    plt.xlim(0.0, nobs+1)
    plt.ylim(1.-ylim, 1.+ylim)

    handles = [ratio, asymm]
#    handles = [ratio]
    labels = [h.get_label() for h in handles]

    legend = ax.legend(handles, labels, title=r'Closure APPLfast vs. NNLOJET',
                       loc='upper right', numpoints=1, frameon=False)
    legend.get_title().set_fontsize(limfs)

    if args['outputfilename'] is None:
        fignam = proc+'.'+jobn+'.'+obsv+'.'+'scale-no-' + \
            str(i+1)+'.scalecheck-'+str(ylim)+'.png'
    else:
        fignam = args['outputfilename']+'.scale-no-' + \
            str(i+1)+'.scalecheck-'+str(ylim)+'.png'
    plt.savefig(fignam)

#    plt.show()

exit(0)
