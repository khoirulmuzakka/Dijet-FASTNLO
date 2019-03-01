#!/usr/bin/env python2
#-*- coding:utf-8 -*-
import argparse
import fastnlo
import glob, os, pylab, sys
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import pandas as pd
# from copy import deepcopy
from matplotlib import cm
from fastnlo import fastNLOLHAPDF
from fastnlo import SetGlobalVerbosity
from fnlo_parsers import FNLOYodaOutputParser
import re
from StringIO import StringIO

# Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('-a','--afile', default='fastnlo', required=True,
                    help='Base name for fastNLO interpolation result, either with log or tab.gz extension.')
parser.add_argument('-b','--bfile', default='theory', required=True,
                    help='Theory result, either with dat, out, or yoda extension.')
parser.add_argument('-c', '--cfile', default='', required=False, type=str,
                    help='Base name addition for comparisons. ' 'Default: none')
parser.add_argument('-o', '--outfiles', default='fnlo-compare', required=False, type=str,
                    help='Start of output filename for plots. ' 'Default: fnlo-compare')

#parse arguments
args   = vars(parser.parse_args())
namesp = parser.parse_args()

# Plot style settings
# Location of matplotlibrc
# print mpl.matplotlib_fname()
# List active style settings
# print mpl.rcParams
# Needs matplotlib > 1.3
# plt.style.use('ggplot')
# plt.style.use('presentation')
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (16, 12),
          'axes.labelsize':  'x-large',
          'axes.titlesize':  'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large'}
pylab.rcParams.update(params)

# Input files
afile = args['afile']
bfile = args['bfile']
cfile = args['cfile']
lgtit = afile
ext1 = [ 'log', 'tab.gz' ]
ext2 = [ 'dat', 'out', 'yoda' ]
aext = ''
bext = ''
cext = ''
for ext in ext1:
    tmp = afile+'.'+ext
    print 'Searching for afile ',tmp
    if os.path.isfile(tmp):
        afile = tmp
        aext  = ext
        break
for ext in ext2:
    tmp = bfile+'.'+ext
    print 'Searching for bfile ',tmp
    if os.path.isfile(tmp):
        bfile = tmp
        bext  = ext
        break
for ext in ext2:
    tmp = cfile+'.'+ext
    print 'Searching for cfile ',tmp
    if os.path.isfile(tmp):
        cfile = tmp
        cext  = ext
        break

if aext:
    print 'Reading fastNLO results from:',afile
else:
    print 'ERROR! No fastNLO file found for arguments -a = ', args['afile'],'-c = ',args['cfile']
    exit(1)
if cext:
    print 'Reading statistical uncertainties for fastNLO from ',cfile
else:
    print 'WARNING! No file with statistical uncertainty found for argument -a = ', args['afile']
if bext:
    print 'Reading theory results from:',bfile
else:
    print 'ERROR! No theory file found for argument -b = ', args['bfile']
    exit(1)

# Output files
outfil = args['outfiles']

# Read fastNLO log file
if aext=='log':
    acols = FNLOYodaOutputParser(afile)
    atab  = acols.get_table()
#    print atab

xsa = np.array(atab["cross_section"])
dxsalr = np.array(atab["lower_uncertainty"])
dxsaur = np.array(atab["upper_uncertainty"])
xsal = xsa * (1+dxsalr)
xsau = xsa * (1+dxsaur)

# Read corresponding statistical uncertainty from NNLOJET dat file
if cext=='dat':
    ccols = np.loadtxt(cfile,usecols=range(3,5))
xs_dat  = np.array(ccols[:,0])
dxs_dat = np.array(ccols[:,1])
dxsar   = np.divide(dxs_dat, xs_dat, out=np.ones_like(dxs_dat), where=xs_dat!=0)
dxsa    = dxsar*xsa

# Read Joao's output
if bext=='out':
    bcols = np.genfromtxt(bfile,comments='#',names=['xl','xm','xu','xs','dxs','xsl','xsu'])
#    print bcols

xsb = bcols['xs']
dxsb = bcols['dxs']
xsbl = bcols['xsl']
xsbu = bcols['xsu']
dxsblr = (xsbl-xsb)/xsb
dxsbur = (xsbu-xsb)/xsb

# Calculate ratios
r_abl = np.divide(xsal, xsbl, out=np.ones_like(xsal), where=xsbl!=0)
r_ab  = np.divide(xsa,  xsb , out=np.ones_like(xsa), where=xsb!=0)
r_abu = np.divide(xsau, xsbu, out=np.ones_like(xsau), where=xsbu!=0)
nobs  = r_ab.size
dxsbr = dxsb/xsb
dxsr_ab = np.sqrt(dxsar*dxsar + dxsbr*dxsbr)

#print dxsar
#print dxsbr
#print dxsar/dxsbr

# Plots
titwgt = 'bold'
limfs = 'x-large'

# Calculate statistical uncertainties
#dr_ablo  = np.multiply(r_ablo[i],dsta)
#dr_abnlo = np.multiply(r_abnlo[i],dsta)
#r_data   = np.divide(xsd, xsb_nlo[i], out=np.ones_like(xsd), where=xsb_nlo[i]!=0)
#dr_data  = np.multiply(r_data,dstd)

fig = plt.figure()
gs  = gridspec.GridSpec(2,1,height_ratios=[3, 1],hspace=0.08)
ax  = plt.subplot(gs[0])
axr = plt.subplot(gs[1])

xaxe = 'bins' # x axis with bin numbers ('bins') or physics observable
if xaxe=='bins':
    xscl = 'linear'
    xlab = 'Observable bin index'
    xmin = 0
    xmax = nobs+1
else:
    xscl = 'log'
    xlab = xaxe
    xmin = xl[0]
    xmax = xu[nobs-1]

ax.set_xscale(xscl)
axr.set_xscale(xscl)
ax.set_yscale('log')
axr.set_yscale('linear')

axr.set_xlabel(xlab, horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
axr.set_ylabel('Grid/NNLOJET', horizontalalignment='center', verticalalignment='center', labelpad=20)
ax.set_ylabel(r'$\bf d\sigma/dX$ (pb/[X])', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, labelpad=20)
ax.set_xticklabels([])

ax.set_title(r'Interpolation vs. direct calculation', fontweight=titwgt, y=1.05)
axr.set_title(r'Ratio for {}'.format(cfile), fontweight=titwgt)
plt.axhline(y=1.01, linestyle='--', linewidth=1.0, color='black')
plt.axhline(y=0.99, linestyle='--', linewidth=1.0, color='black')
plt.fill_between([0.0,nobs+1],0.99,1.01, color='black', alpha=0.1)
#    plt.text(nobs+1.5,1.0,u'$\pm$1‰',fontsize=limfs)
plt.text(nobs+1.5,1.0,u'$\pm$1%',fontsize=limfs)

dx = 1./4.
x  = np.arange(1   , nobs+1.e-6)
xa = np.arange(1-dx, nobs-dx+1.e-6)
xb = np.arange(1+dx, nobs+dx+1.e-6)

xcb = ax.errorbar(x, xsb, yerr=dxsb, marker='s', markersize=12, linestyle='none', fillstyle='none', label=u'$\sigma(\mu_0)\,\pm\,$stat. (dir.)', color='red')
xca = ax.errorbar(x, xsa, yerr=dxsa, marker='x', markersize=12, linestyle='none', label=u'$\sigma(\mu_0)\,\pm\,$stat. (int.)', color='black')
xcsc = ax.fill_between(x,xsb-1.*(xsb-xsbl),xsb+1.*(xsbu-xsb),alpha=0.5,edgecolor='#ffffff',facecolor='pink',label='scale envelope (dir.)')

rl = axr.errorbar(xa, r_abl, yerr=dxsar*r_abl, marker='v', linestyle='none', label=u'$\sigma(\mu_\mathrm{l,6P})\,\pm\,$stat. (int.)', color='blue')
rc = axr.errorbar(x,  r_ab,  yerr=dxsr_ab*r_ab , marker='d', linestyle='none', label=u'$\sigma(\mu_0)\,\pm\,$stat. (tot.)', color='green')
ru = axr.errorbar(xb, r_abu, yerr=dxsbr*r_abu, marker='^', linestyle='none', label=u'$\sigma(\mu_\mathrm{u,6P})\,\pm\,$stat. (dir.)', color='red')
#axr.fill_between(x,xsbl/xsb,xsbu/xsb)

ax.set_xlim(0.0,nobs+1)
axr.set_xlim(0.0,nobs+1)
axr.set_ylim(0.90,1.10)

handles = [xcb,xcsc,xca,rl,rc,ru]
labels  = [h.get_label() for h in handles]

#legend = ax.legend()
legend = ax.legend(handles, labels, title=r'{}'.format(lgtit), loc='best', numpoints=1, handlelength=1.5)
legend.get_title().set_fontsize(limfs)

fignam = outfil+'.png'
plt.savefig(fignam)

exit(0)
