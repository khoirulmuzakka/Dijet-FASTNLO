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

# Get process order/type from cmdline argument
scen = 'fnl2332d'
proc = '1jet'
ecms = '7TeV'
obsv = 'xptj0_y1'
scle = 0
comb = ''
labl = 'comb'
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
if len(sys.argv) > 1:
    scen = sys.argv[1]
if len(sys.argv) > 2:
    proc = sys.argv[2]
if len(sys.argv) > 3:
    ecms = sys.argv[3]
if len(sys.argv) > 4:
    obsv = sys.argv[4]
if len(sys.argv) > 5:
    scle = int(sys.argv[5])
if len(sys.argv) > 6:
    comb = '.'+sys.argv[6]
    labl = sys.argv[6]

# Prepare result arrays
xl = []       # left bin border
xm = []       # bin "center"
xu = []       # right bin border
xs_nnlo  = [] # NNLOJET results
dxs_nnlo = [] # NNLOJET uncertainty
xs_fnlo  = [] # fastNLO results
sclecol  = 2 * (scle+1) + 1 # scale columns are 3, 5, 7, ... when starting count with zero
sclenam  = ['scl2', 'fscl1', 'fscl2', 'fscl3', 'fscl4', 'fscl5', 'fscl6']
print 'Selected scale column is:', sclecol
print 'Selected scale name is:', sclenam[scle]

# Read summed-up cross sections from NNLOJET dat files
sums = ['LO', 'NLO', 'NNLO', 'R', 'V', 'RRa', 'RRb', 'RV', 'VV']
#sums = ['LO', 'R', 'V', 'RRa', 'RRb', 'RV', 'VV']
#sums = ['LO', 'NLO', 'NNLO']

sumfiles = []
for i in range(len(sums)):
    print 'Reading', sums[i], 'NNLOJET cross sections'
    if labl=='single':
        sumfiles.append(sums[i]+'/'+scen+'.'+proc+'.'+sums[i]+'-'+ecms+'.'+obsv+comb+'.dat')
    else:
        sumfiles.append('Combined/Final/'+sums[i]+'.'+obsv+'_scl0.dat')
    print 'The filename is:', sumfiles[i]
    if i == 0:
        xl.append(np.loadtxt(sumfiles[i],usecols=(0,)))
        xu.append(np.loadtxt(sumfiles[i],usecols=(2,)))
    xs_nnlo.append(np.loadtxt(sumfiles[i],usecols=(sclecol,)))
    dxs_nnlo.append(np.loadtxt(sumfiles[i],usecols=(sclecol+1,)))

xl = np.array(xl)
xu = np.array(xu)
xs_nnlo  = np.array(xs_nnlo)/1000.  # Conversion of fb to pb
dxs_nnlo = np.array(dxs_nnlo)/1000. # Conversion of fb to pb

# Determine no. of observable bins
nobs = xl.size
print 'Number of observable bins: ', nobs

# Read cross sections from pre-evaluated fastNLO tables
logfiles = []
for i in range(len(sums)):
    print 'Reading', sums[i], 'fastNLO cross sections'
    if 'NLO' in sums[i]:
        logfiles.append(sums[i]+'_Combined/'+scen+'.'+proc+'.'+sums[i]+'-'+ecms+'.'+obsv+'_'+sclenam[scle]+'.log')
    else:
        logfiles.append(sums[i]+'/'+scen+'.'+proc+'.'+sums[i]+'-'+ecms+'.'+obsv+comb+'_'+sclenam[scle]+'.log')
    print 'The filename is:', logfiles[i]
    if sums[i]=='NLO':
        xs_fnlo.append(np.loadtxt(logfiles[i],usecols=(7,),comments=['#',' #','C','L']))
    elif sums[i]=='NNLO':
        xs_fnlo.append(np.loadtxt(logfiles[i],usecols=(8,),comments=['#',' #','C','L']))
    else:
        xs_fnlo.append(np.loadtxt(logfiles[i],usecols=(6,),comments=['#',' #','C','L']))

xs_fnlo = np.array(xs_fnlo)

# Ratio and asymmetry
r_f2nb  = xs_fnlo/xs_nnlo
dr_f2nb = dxs_nnlo/xs_nnlo*(xs_fnlo/xs_nnlo)
a_f2nb  = (xs_fnlo - xs_nnlo)/(xs_fnlo + xs_nnlo)

print r_f2nb

# Ratio plots
for i in range(len(sums)):
    fig = plt.figure(figsize=(16,12))
    ax  = fig.gca()
    plt.title(r'Ratio: {} {} scale {} {} for {} at {} ({})'.format(proc, sums[i], scle, labl, obsv, ecms, scen), fontsize='x-large', fontweight='bold', loc='left')
    plt.xlabel('Observable bin index', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
    plt.ylabel('fastNLO/NNLOJET', horizontalalignment='right', x=1.0, verticalalignment='bottom', y=1.0)
    plt.axhline(y=1.001, linestyle='--', linewidth=1.0, color='black')
    plt.axhline(y=0.999, linestyle='--', linewidth=1.0, color='black')
    plt.fill_between([0.,nobs+1],0.999,1.001, color='black', alpha=0.1)

    xc = np.arange(nobs)+1

    rwgt = plt.errorbar(xc, r_f2nb[i], yerr=dr_f2nb[i], marker='v', linestyle='none', label=r'r $\pm$ $\Delta{}$r(stat)', color='blue')
    plt.errorbar(xc, r_f2nb[i], yerr=dr_f2nb[i], marker='.', linestyle='none', color='blue')
    plt.xlim(0.,nobs+1)
    llim = 0.99
    ulim = 1.01
    if ulim-llim>0.1:
        plt.text(34.3,0.999,u'±1‰')
    else:
        plt.text(34.3,1.00085,u'+1‰')
        plt.text(34.3,0.99885,u'–1‰')
    plt.ylim(llim,ulim)
    handles = [rwgt]
    labels  = [h.get_label() for h in handles]
    legend = ax.legend(handles, labels, title=r'No. of points = {}'.format(nobs), loc='upper left', numpoints=1, handlelength=0)

    fignam = scen+'.'+proc+'.'+sums[i]+'-'+ecms+'.'+obsv+comb+'.'+sclenam[scle]+'_interpol_ratio'+'.png'
    plt.savefig(fignam)
#plt.show()

# Asymmetry plots
    fig = plt.figure(figsize=(16,12))
    ax  = fig.gca()
    plt.title(r'Asymmetry: {} {} scale {} {} for {} at {} ({})'.format(proc, sums[i], scle, labl, obsv, ecms, scen), fontsize='x-large', fontweight='bold', loc='left')
    plt.xlabel('Observable bin index', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
    plt.ylabel('(fastNLO-NNLOJET)/(fastNLO+NNLOJET)', horizontalalignment='right', x=1.0, verticalalignment='bottom', y=1.0)
    plt.axhline(y=-0.001, linestyle='--', linewidth=1.0, color='black')
    plt.axhline(y=+0.001, linestyle='--', linewidth=1.0, color='black')
    plt.fill_between([0.,nobs+1],-0.001,0.001, color='black', alpha=0.1)
    plt.text(34,+0.001,u'+1‰')
    plt.text(34,-0.001,u'–1‰')

    xc = np.arange(nobs)+1

    awgt = plt.errorbar(xc, a_f2nb[i], yerr=0.*a_f2nb[i], marker='<', linestyle='none', label='$\mu$ $\pm$ $\Delta\mu$', color='orange')
    plt.errorbar(xc, a_f2nb[i], yerr=0.*a_f2nb[i], marker='.', linestyle='none', color='orange')
    plt.xlim(0.,nobs+1)
    plt.ylim(-0.01,0.01)
    handles = [awgt]
    labels  = [h.get_label() for h in handles]
    legend = ax.legend(handles, labels, title=r'No. of points = {}'.format(nobs), loc='upper left', numpoints=1, handlelength=0)

    fignam = scen+'.'+proc+'.'+sums[i]+'-'+ecms+'.'+obsv+comb+'.'+sclenam[scle]+'_interpol_asymm'+'.png'
#    plt.savefig(fignam)
#plt.show()

exit(0)
