#!/usr/bin/env python2
#-*- coding:utf-8 -*-
import argparse
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

parser = argparse.ArgumentParser()

#add arguments (three .dat and one .log files)
parser.add_argument('-d','--datfiles', default='file.dat', required=True, nargs=3,
                    help='.dat files for evaluation.')
parser.add_argument('-l', '--logfile', default='file.log', required=True, nargs='?',
                    help='.log file, need one.')
                    
parser.add_argument('-o', '--outputfilename', required=False, nargs='?', type=str,
                    help='Customise the first part of the output filename.'
                            'Default: Same structure as datfile name.')

#parse arguments
args = vars(parser.parse_args())
namesp = parser.parse_args()


#take care of the input files
log0file = args['logfile']
#create dictionary for order and correspoinding datfiles
pathdatfiles = {}
datfiles = {} ###### use basename!!!

for arg in args['datfiles']:
    arg0 = os.path.basename(arg)
    pre = arg0.split('.')[1]    #new datfile name-format with <proc>.<ord>....dat
    if pre in ['LO', 'NLO', 'NNLO']:
        pathdatfiles[pre]=arg
        datfiles[pre]=arg0
    else:
        print "Problem with datfiles."
        sys.exit("Exit: Input ERROR.")

print '\n', "datfiles: ", datfiles, '\n'
if len(datfiles)!=3:
    print "Check whether all required datfiles are given."
    sys.exit("ERROR: Wrong amount of .dat files per order.")
    #for instance if user gives two times the NLO .dat file and no NNLO


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


#arguments
print "log0file: ", log0file, '\n'
log0base = os.path.basename(log0file)
log0args = log0base.split(".") #array containing filename-parts

#some default values
seed = '' #default
nscl = 6 # central + nscl fixed scales
xaxe = 'bins' # x axis with bin numbers ('bins') or physics observable

#arguments from logfile
print 'arguments in logfile name: ', len(log0args)
print log0args
if len(log0args)==4:
    proc, jobn, obsv, ext = log0args
elif len(log0args)==5:
    proc, jobn, kinn, obsv, ext = log0args
elif len(log0args)==6:
    proc, jobn, kinn, obsv, seed, ext = log0args


print 'proc: ', proc
print 'jobn: ', jobn
print 'obsv: ', obsv
print '\n'


# Extract data set from job type (substring after first '-')
print jobn.split('-')
pord = jobn.split('-')[0]
#dset = jobn.split('-')[1]  #not existent in some(!) newer tabels (ZJtriple)
if pord != 'NNLO':
    print 'Error! Need NNLO results for this script, aborted'
    print 'jobn is ', jobn
    exit(1)

# Prepare result arrays
xl = [] # left bin border
xm = [] # bin "center"
xu = [] # right bin border

# Read binning and cross sections from NNLOJET dat file
# Loop over orders
orders = [ 'LO', 'NLO', 'NNLO' ]
for order in orders:
    #datfile = datfiles[order] (this is only basename)
    datfile = pathdatfiles[order] #(total path)    
    print 'Reading from NNLOJET dat file ', datfile
    xs_all = np.loadtxt(datfile,usecols=range(0,17))
    if order=='LO':
        xl = xs_all[:,0] # Binning must be identical for all orders
        xm = xs_all[:,1]
        xu = xs_all[:,2]
        xs_lo    = xs_all[:,3]/1000. # Conversion of fb to pb
        dxs_lo   = xs_all[:,4]/1000.
    elif order=='NLO':
        xs_nlo   = xs_all[:,3]/1000.
        dxs_nlo  = xs_all[:,4]/1000.
    else:
        xs_nnlo  = xs_all[:,3]/1000.
        dxs_nnlo = xs_all[:,4]/1000.

xs_lo    = np.array(xs_lo)
dxs_lo   = np.array(dxs_lo)
xs_nlo   = np.array(xs_nlo)
dxs_nlo  = np.array(dxs_nlo)
xs_nnlo  = np.array(xs_nnlo)
dxs_nnlo = np.array(dxs_nnlo)

# Determine no. of observable bins
nobs = xl.size
print 'Number of observable bins: ', nobs
xb = np.arange(1, nobs+1.e-6)

#print 'LO', xs_lo
#print 'NLO', xs_nlo
#print 'NNLO', xs_nnlo

# Read cross sections from pre-evaluated fastNLO tables
ordcol = { 'LO': 6, 'NLO': 7, 'NNLO': 8 }
#logfile = proc+'.'+jobn+'.'+obsv+'_0.log'
##print 'Reading from fastNLO log file ', log0file

# Skip all lines starting with "#", "C", or "L" as first non-whitespace character
### why don't we start these lines (C, L...) also with a "#" ???
# Read first nobs values for central scale result
'''
with open(log0file, 'r') as f:
    data = re.sub(r'\s*[#CL].*', '', f.read())
    lo   = np.genfromtxt(StringIO(data),usecols=ordcol['LO'],)
    nlo  = np.genfromtxt(StringIO(data),usecols=ordcol['NLO'],)
    nnlo = np.genfromtxt(StringIO(data),usecols=ordcol['NNLO'],)
    xs_flo   = []
    xs_fnlo  = []
    xs_fnnlo = []
    for irow in range(nobs):
        xs_flo.append(lo[irow])
        xs_fnlo.append(nlo[irow])
        xs_fnnlo.append(nnlo[irow])

xs_flo    = np.array(xs_flo)
xs_fnlo   = np.array(xs_fnlo)
xs_fnnlo  = np.array(xs_fnnlo)

print "xs_flo: \n", xs_flo
print "xs_fnlo: \n", xs_fnlo
print "xs_fnnlo: \n", xs_fnnlo
'''

# K factors
kn_nlo  = np.divide(xs_nlo, xs_lo, out=np.ones_like(xs_nlo), where=xs_lo!=0)
##kf_nlo  = np.divide(xs_fnlo, xs_flo, out=np.ones_like(xs_fnlo), where=xs_flo!=0)
kn_nnlo = np.divide(xs_nnlo, xs_nlo, out=np.ones_like(xs_nnlo), where=xs_nlo!=0)
##kf_nnlo = np.divide(xs_fnnlo, xs_fnlo, out=np.ones_like(xs_fnnlo), where=xs_fnlo!=0)

#print 'knn', kn_nlo
#print 'knf', kf_nlo
#print 'knnn', kn_nnlo
#print 'knnf', kf_nnlo
print "-------------"

# Rel. stat. uncertainty of K factors from NNLOJET, assume dominant uncertainty from numerator!
dst_nlo  = np.divide(dxs_nlo, xs_nlo, out=np.ones_like(dxs_nlo), where=xs_nlo!=0)
dst_nnlo = np.divide(dxs_nnlo, xs_nnlo, out=np.ones_like(dxs_nnlo), where=xs_nnlo!=0)
dk_nlo   = np.multiply(kn_nlo,dst_nlo)
dk_nnlo  = np.multiply(kn_nnlo,dst_nnlo)

##print dst_nlo
##print dst_nnlo

##print dk_nlo
##print dk_nnlo


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
ax.set_ylabel(r'$\bf d\sigma/dp_T$ [pb/GeV]', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, labelpad=20)
axr.set_ylabel('K factor', horizontalalignment='center', verticalalignment='center', labelpad=20)
ax.set_xticklabels([])

ax.set_title(r'Cross sections and K factors: {} {} {}'.format(proc, jobn, obsv), fontweight=titwgt, y=1.05)

axhandles = []

xlo   = xs_lo
xnlo  = xs_nlo
xnlo  = xs_nnlo
knlo  = kn_nlo
knnlo = kn_nnlo

abslo   = ax.errorbar(x, xs_lo, xerr=dx, capsize=0, yerr=dxs_lo, marker=',', linestyle='none', label=r'LO $\bf\pm\Delta_{stat}$', color=col1)
axhandles.append(abslo)
absnlo  = ax.errorbar(x, xs_nlo, xerr=dx, capsize=0, yerr=dxs_nlo, marker='s', linestyle='none', label=r'NLO $\bf\pm\Delta_{stat}$', color=col3)
axhandles.append(absnlo)
absnnlo = ax.errorbar(x, xs_nnlo, xerr=dx, capsize=0, yerr=dxs_nnlo, marker='o', linestyle='none', label=r'NNLO $\bf\pm\Delta_{stat}$', color=col2, mfc='none', mec=col2, mew=3)
axhandles.append(absnnlo)

ax.set_xlim(xmin,xmax)
axr.set_xlim(xmin,xmax)
axr.set_ylim(0.8,1.6)
#axr.fill_between(x, 1.-abs(dst_nnlo[i]), 1.+abs(dst_nnlo[i]), edgecolor=col1, facecolor=col1b, alpha=0.5)
axr.axhline(1.0,color='black')

rkn  = axr.errorbar(x, knlo, yerr=dk_nlo, marker='s', linestyle='none', label='', color=col3)
rknn = axr.errorbar(x, knnlo, yerr=dk_nnlo, marker='o', linestyle='none', label='', color=col2)
axrhandles = [ rkn, rknn ]

axlabels  = [h.get_label() for h in axhandles]
axrlabels = [h.get_label() for h in axrhandles]

legend = ax.legend(axhandles, axlabels, title=r'NNLOJET inclusive jet $\bf p_T$', loc='upper right', numpoints=1, frameon=False)
legend.get_title().set_fontsize(limfs)

if xaxe=='bins':
    if args['outputfilename'] is None:
        fignam = proc+'.'+jobn+'.'+obsv+'.'+'absolute-kfac-index'+'.png'
    else:
        fignam = args['outputfilename']+'.absolute-kfac-index.png'
else:
    if args['outputfilename'] is None:
        fignam = proc+'.'+jobn+'.'+obsv+'.'+'absolute-kfac-obs'+'.png'
    else:
        fignam = args['outputfilename']+'.absolute-kfac-obs.png'

print 'Writing figure', fignam
plt.savefig(fignam)

exit(0)
