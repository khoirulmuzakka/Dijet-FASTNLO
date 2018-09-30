#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import yoda
import fastnlo

import argparse
import re
import os

parser = argparse.ArgumentParser(description='Compares cross sections calculated with fastnlo with those calculated by rivet.')
parser.add_argument('basename',default='yb0_ystar0_zpt',help='This is the basename of the table, it will be used to find the corresponding table file name and the yoda histogram in the yoda file provided')
parser.add_argument('--path',default='.',help='Path that will be used to search for tables')
parser.add_argument('--yoda',default='',help='File name of the yoda file that should contain cross sections calculated by rivet')
parser.add_argument('--pdf',default='CT14nlo',help='Name of the pdf to use in the cross section calculation with fastnlo. Note this MUST be the same pdf that was used by rivet to calculate the yoda histograms')
parser.add_argument('--out',default='',help='Output file name')

args = parser.parse_args()

for f in os.listdir(args.path):
    if re.search(args.basename+'.*\.tab.gz',f):
        table = args.path+'/'+f

fnlo = fastnlo.fastNLOLHAPDF(table,args.pdf,0)

if args.yoda:
    yodafile = args.yoda
else:
    yodafile = args.path+'/Rivet.yoda'

for key, value in yoda.read(yodafile).items():
    if re.search(args.basename, key):
        rivet = value

fnlo.SetContributionON(fastnlo.kFixedOrder, fastnlo.kLeading, True)
fnlo.SetContributionON(fastnlo.kFixedOrder, fastnlo.kNextToLeading, True)
fnlo.CalcCrossSection()

xs_rivet = np.array([ b.height for b in rivet.bins ])
xs_fnlo  = np.array(fnlo.GetCrossSection())

print xs_fnlo
print xs_rivet

N = fnlo.GetNDim0Bins()

fig, ax = plt.subplots()
ax.set(ylabel='$\mathrm{xs_{fastNLO}/xs_{Rivet}}$',xlabel='$p_{T}^{Z}$ bin index',title='fastNLO cross check')

ax.yaxis.label.set_size(18)
ax.set(ylim=[0.99,1.01])
ax.grid(b=True,which='major')

ax.set_yticks([1-0.001,1+0.001],minor=True)
ax.set_yticklabels([u' -1 \u2030',u'+1 \u2030'],minor=True)
ax.tick_params(axis='y',which='minor',labelleft=False,labelright=True,)
ax.grid(b=True,which='minor',linestyle='--',color='k')
ax.axhspan(1-0.001,1+0.001,color='b',alpha=.1)

ax.plot( range(1,N+1), xs_fnlo/xs_rivet, 'bo')

if args.out:
    out = args.out
else:
    out = args.path+'/'+args.basename+'_check.png'
fig.savefig(out)
