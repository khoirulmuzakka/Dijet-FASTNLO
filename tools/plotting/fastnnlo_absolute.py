#!/usr/bin/env python2
#-*- coding:utf-8 -*-
import glob, os, pylab, sys
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
# from copy import deepcopy
from matplotlib import cm
from fastnlo import fastNLOLHAPDF
from fastnlo import SetGlobalVerbosity
import re
from StringIO import StringIO
#import warnings
#warnings.filterwarnings("error")

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
pylab.rcParams.update(params)

# Default arguments
proc = '1jet'
jobn = 'LO-CMS7'
kinn = 'vBa'
obsv = 'fnl2332d_xptji_y1'
seed = ''
nscl = 6  # central + nscl fixed scales
xaxe = 'bins' # x axis with bin numbers ('bins') or physics observable
# Get from cmdline argument
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
if len(sys.argv) > 1:
    proc = sys.argv[1]
if len(sys.argv) > 2:
    jobn = sys.argv[2]
if len(sys.argv) > 3:
    kinn = sys.argv[3]
if len(sys.argv) > 4:
    obsv = sys.argv[4]
if len(sys.argv) > 5:
    seed = sys.argv[5]
if len(sys.argv) > 6:
    if not sys.argv[6]=='_': 
        nscl = int(sys.argv[6])
if len(sys.argv) > 7:
    xaxe = sys.argv[7]

# Extract order/contribution from job type (substring before first '-')
order  = jobn.split('-')[0]
ordcol = 6
if order == 'NLO':
    ordcol = 7
elif order == 'NNLO':
    ordcol = 8

# Prepare result arrays
xl = []      # left bin border
xm = []      # bin "center"
xu = []      # right bin border
xs_nnlo = [] # NNLOJET results
xs_fnlt = [] # fastNLO results
xs_fnll = [] # Pre-evaluated fastNLO results (log file)

# Read binning and cross sections from NNLOJET dat file
if kinn == '_':
    datfile = order+'.'+obsv+'.dat'
elif seed == '' or seed == '_':
    datfile = order+'.'+obsv+'.dat'
else:
    datfile = proc+'.'+jobn+'.'+kinn+'.'+obsv+'.'+seed+'.dat'
# Merge multiple occurences of '.' in filename by one '.'
#datfile = re.sub('\.+','.', datfile)
print 'Reading from NNLOJET dat file ', datfile
xs_all = np.loadtxt(datfile,usecols=range(0,17))
xl = xs_all[:,0]
xm = xs_all[:,1]
xu = xs_all[:,2]
xs_nnlo  = []
dxs_nnlo = []
for i in range(nscl+1):
    xs_nnlo.append(xs_all[:,2*i+3]/1000) # Conversion of fb to pb
    dxs_nnlo.append(xs_all[:,2*i+4]/1000)

# Determine no. of observable bins
nobs = xl.size
print 'Number of observable bins: ', nobs
xb = np.arange(1, nobs+1.e-6)

# Evaluate cross sections from fastNLO tables
# INFO=0, WARNING=1
#SetGlobalVerbosity(1)
#fnlotabs = glob.glob(pord+'/'+scen+'.'+proc+'.'+pord+'-'+ecms+'.???.'+obsv+'.*.tab.gz')
#fnlotabs.sort()
#for fnlotab in fnlotabs:
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
if kinn == '_':
    log0file = proc+'.'+jobn+'.'+obsv+'_0.log'
    log6file = proc+'.'+jobn+'.'+obsv+'_6.log'
elif seed == '' or seed == '_':
    log0file = proc+'.'+jobn+'.'+kinn+'.'+obsv+'_0.log'
    log6file = proc+'.'+jobn+'.'+kinn+'.'+obsv+'_6.log'
else:
    log0file = proc+'.'+jobn+'.'+kinn+'.'+obsv+'.'+seed+'_0.log'
    log6file = proc+'.'+jobn+'.'+kinn+'.'+obsv+'.'+seed+'_6.log'
for file in [log0file, log6file]:
    print 'Reading from fastNLO log file ', file
    # Skip all lines starting with "#", "C", or "L" as first non-whitespace character
    with open(file, 'r') as f:
        data = re.sub(r'\s*[#CL].*', '', f.read())
        all = np.genfromtxt(StringIO(data),usecols=(ordcol,))
        # Read the (nscl+1)*nobs values into nscl arrays of nobs entries
        ns = nscl
        if file==log0file:
            ns = 1
        for i in range(ns):
            a = []
            for j in range(nobs):
                ind = i*nobs+j
                a.append(all[ind])
            xs_fnll.append(a)

xs_fnll = np.array(xs_fnll)

# Ratio and asymmetry+1
r_nn2nn = np.ones_like(xs_nnlo)
r_fl2nn = np.divide(xs_fnll, xs_nnlo, out=np.ones_like(xs_fnll), where=xs_nnlo!=0)
a_fl2nn = np.divide(xs_fnll-xs_nnlo, xs_fnll+xs_nnlo, out=np.zeros_like(xs_fnll), where=xs_nnlo!=0) + 1.

# Rel. stat. uncertainty from NNLOJET
dst_nnlo = np.divide(dxs_nnlo, xs_nnlo, out=np.ones_like(dxs_nnlo), where=xs_nnlo!=0)
dr_fl2nn = np.multiply(r_fl2nn,dst_nnlo)
da_fl2nn = np.multiply(a_fl2nn,dst_nnlo)

# Prepare plotting
titwgt = 'bold'
limfs = 'x-large'
sclnam = [r'$\bf p_{T,max}$', r'$\bf\mu_r = \mu_f$ small', r'$\bf\mu_r = \mu_f$ large', r'$\bf\mu_r < \mu_f$', r'$\bf\mu_r > \mu_f$', r'$\bf\mu_r \ll \mu_f$', r'$\bf\mu_r \gg \mu_f$']
col1  = 'olivedrab'
col1b = 'yellowgreen'
col2 = 'orangered'
col3 = 'dodgerblue'

if xaxe=='bins':
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


for i in range(nscl+1):

# Absolute predictions
    fig = plt.figure()
    gs  = gridspec.GridSpec(2,1,height_ratios=[3, 1],hspace=0.08)
    ax  = plt.subplot(gs[0])
    axr = plt.subplot(gs[1])

    ax.set_xscale(xscl)
    axr.set_xscale(xscl)
    ax.set_yscale('log')
    axr.set_yscale('linear')

    axr.set_xlabel(xlab, horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
    ax.set_ylabel(r'$\bf\left|d\sigma/dp_T\right|$ [pb/GeV]', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, labelpad=20)
    axr.set_ylabel('Ratio', horizontalalignment='center', verticalalignment='center', labelpad=20)
    ax.set_xticklabels([])

    if seed == '' or seed == '_':
        ax.set_title(r'Merged grid: {} {} {} for scale choice {}'.format(proc, jobn, obsv, sclnam[i]), fontweight=titwgt, y=1.05)
    else:
        ax.set_title(r'Single grid: {} {} {} for scale choice {}'.format(proc, jobn, obsv, sclnam[i]), fontweight=titwgt, y=1.05)

    axhandles = []
# Booleans for + and - x sections
    xmp = xs_nnlo[i]>0
    xmn = xs_nnlo[i]<0
    if len(x[xmp]):
        abs1p = ax.errorbar(x[xmp], +xs_nnlo[i][xmp], xerr=dx[xmp], capsize=0, yerr=dxs_nnlo[i][xmp], marker=',', linestyle='none', label=r'NNLOJET $\bf\pm\Delta_{stat}$', color=col1)
        axhandles.append(abs1p)
    if len(x[xmn]):
        abs1n = ax.errorbar(x[xmn], -xs_nnlo[i][xmn], xerr=dx[xmn], capsize=0, yerr=dxs_nnlo[i][xmn], marker=',', linestyle='none', label=r'NNLOJET', color=col1)
#        axhandles.append(abs1n)

    xmp = xs_fnll[i]>0
    xmn = xs_fnll[i]<0
    if len(x[xmp]):
        abs2p = ax.errorbar(x[xmp], +xs_fnll[i][xmp], marker='d', linestyle='none', label=r'APPLfast grid ($\bf\sigma>0$)', color=col2)
        axhandles.append(abs2p)
    if len(x[xmn]):
        abs2n = ax.errorbar(x[xmn], -xs_fnll[i][xmn], marker='d', linestyle='none', label=r'APPLfast grid ($\bf\sigma<0$)', color=col3, markerfacecolor=col3)
        axhandles.append(abs2n)

    ax.set_xlim(xmin,xmax)
    axr.set_xlim(xmin,xmax)
    axr.set_ylim(0.99,1.01)
    axr.fill_between(x, 1.-abs(dst_nnlo[i]), 1.+abs(dst_nnlo[i]), edgecolor=col1, facecolor=col1b, alpha=0.5)
    axr.axhline(1.0,color=col1)

    rnn = axr.errorbar(x, r_nn2nn[i], marker=',', linestyle='none', label='', color=col1)
    rfp = axr.errorbar(x[xmp], r_fl2nn[i][xmp], marker='d', linestyle='none', label='', color=col2)
    rfn = axr.errorbar(x[xmn], r_fl2nn[i][xmn], marker='d', linestyle='none', label='', color=col3, markerfacecolor=col3)
    axrhandles = [ rnn, rfp, rfn ]

    axlabels  = [h.get_label() for h in axhandles]
    axrlabels  = [h.get_label() for h in axrhandles]

    legend = ax.legend(axhandles, axlabels, title=r'Stat. uncertainty & grid closure', loc='upper right', numpoints=1, frameon=False)
    legend.get_title().set_fontsize(limfs)

    if xaxe=='bins':
        fignam = proc+'.'+jobn+'.'+obsv+'.'+'scale-no-'+str(i+1)+'.absolute-index'+'.png'
    else:
        fignam = proc+'.'+jobn+'.'+obsv+'.'+'scale-no-'+str(i+1)+'.absolute-obs'+'.png'
    print 'Writing figure', fignam
    plt.savefig(fignam)

#    plt.show()

exit(0)