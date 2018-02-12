#!/usr/bin/env python2
#-*- coding:utf-8 -*-
import glob, os, pylab, sys
import matplotlib as mpl
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
#warnings.filterwarnings("error")

# Style settings
# Location of matplotlibrc
# print mpl.matplotlib_fname()
# List active style settings
# print mpl.rcParams
# Needs matplotlib > 1.3
# plt.style.use('ggplot')
# plt.style.use('presentation')

# Default arguments
proc = '1jet'
jobn = 'LO-7TeV'
kinn = 'vBa'
obsv = 'xptj0_y1'
seed = 's1234'
fnlo = 'scls'
nscl = 6
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
    fnlo = sys.argv[6]
if len(sys.argv) > 7:
    nscl = sys.argv[7]

# Prepare result arrays
xl = []      # left bin border
xm = []      # bin "center"
xu = []      # right bin border
xs_nnlo = [] # NNLOJET results
xs_fnlt = [] # fastNLO results
xs_fnll = [] # Pre-evaluated fastNLO results (log file)

# Read binning and cross sections from NNLOJET dat file
datfile = proc+'.'+jobn+'.'+kinn+'.'+obsv+'.'+seed+'.dat'
xs_scls = []
print 'Reading from NNLOJET dat file ', datfile
xl.append(np.loadtxt(datfile,usecols=(0,)))
xu.append(np.loadtxt(datfile,usecols=(2,)))
xl = np.array(xl)
xu = np.array(xu)
for i in range(nscl):
    xs_scls.append(np.loadtxt(datfile,usecols=(2*i+5,)))
xs_nnlo = np.array(xs_scls)/1000. # Conversion of fb to pb
#print "xs_nnlo", xs_nnlo

# Determine no. of observable bins
nobs = xl.size
print 'Number of observable bins: ', nobs

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
logfile = proc+'.'+jobn+'.'+kinn+'.'+obsv+'.'+seed+'.log'
print 'Reading from fastNLO log file ', logfile
# Skip all lines starting with "#", "C", or "L" as first non-whitespace character
with open(logfile, 'r') as f:
    data = re.sub(r'\s*[#CL].*', '', f.read())
all = np.genfromtxt(StringIO(data),usecols=(6,))
# Read the nscl*nobs values into nscl arrays of nobs entries
for i in range(nscl):
    a = []
    for j in range(nobs):
        ind = i*nobs+j
        a.append(all[ind])
    xs_fnll.append(a)

xs_fnll = np.array(xs_fnll)
#print "xs_fnll", xs_fnll

r_fl2nn = np.divide(xs_fnll, xs_nnlo, out=np.ones_like(xs_fnll), where=xs_nnlo!=0)
a_fl2nn = np.divide(xs_fnll-xs_nnlo, xs_fnll+xs_nnlo, out=np.zeros_like(xs_fnll), where=xs_nnlo!=0) + 1.

#print 'ratio', r_fl2nn
#print 'asymmetry + 1', a_fl2nn

#exit(0)

for i in range(nscl):

    # Closure plot
    fig = plt.figure(figsize=(16,12))
    ax  = fig.gca()
    plt.title(r'Single grid: {} {} for {} at scale var {}'.format(proc, jobn, obsv, i+1), fontsize='xx-large', fontweight='bold')
    plt.xlabel('Observable bin index', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
    plt.ylabel('Closure fastNLO vs. NNLOJET', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
    plt.axhline(y=1.001, linestyle='--', linewidth=1.0, color='black')
    plt.axhline(y=0.999, linestyle='--', linewidth=1.0, color='black')
    plt.fill_between([0.0,34.0],0.999,1.001, color='black', alpha=0.1)
    plt.text(34.6,1.00085,u'+1‰')
    plt.text(34.7,0.99885,u'–1‰')

    dx = 1./4.
    x  = np.arange(1   , nobs+1.e-6)
    xa = np.arange(1-dx, nobs-dx+1.e-6)
    xb = np.arange(1+dx, nobs+dx+1.e-6)

    asymm = plt.errorbar(x, a_fl2nn[i,], yerr=0.*x, marker='s', linestyle='none', label=r'Asymmetry (fastNLO-NNLOJET)/(fastNLO-NNLOJET) + 1', color='blue')
    ratio = plt.errorbar(x, r_fl2nn[i,], yerr=0.*x, marker='o', linestyle='none', label=r'Ratio fastNLO/NNLOJET', color='orange')

    plt.xlim(0.0,34.0)
    plt.ylim(0.99,1.01)

    handles = [ratio,asymm]
    labels  = [h.get_label() for h in handles]

    legend = ax.legend(handles, labels, title=r'Closure of fastNLO vs. NNLOJET', loc='upper left', numpoints=1, handlelength=0)
    legend.get_title().set_fontsize(20)

    fignam = proc+'.'+jobn+'.'+obsv+'.'+'scale-no-'+str(i+1)+'.scalecheck'+'.png'
    plt.savefig(fignam)

#    plt.show()

exit(0)
