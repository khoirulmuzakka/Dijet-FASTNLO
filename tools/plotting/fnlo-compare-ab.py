#!/usr/bin/env python2
#-*- coding:utf-8 -*-
#
# Make comparison plots between two predictions
#
# created by K. Rabbertz: 14.09.2018
#
#-----------------------------------------------------------------------
#
# Use matplotlib with Cairo offline backend for png, eps, or svg output
import matplotlib as mpl
mpl.use('Cairo')
import argparse
import glob
import os
import re
import sys
from StringIO import StringIO
import matplotlib.lines as mpllines
import matplotlib.gridspec as gridspec
import matplotlib.patches as mplpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import (
    FormatStrFormatter, ScalarFormatter, AutoMinorLocator, MultipleLocator)
from matplotlib import cm
# numpy
import numpy as np
# pandas
#import pandas as pd
#from fnlo_parsers import FNLOCppReadOutputParser
#from fnlo_parsers import FNLOYodaOutputParser
# fastNLO
#import fastnlo
#from fastnlo import fastNLOLHAPDF
#from fastnlo import SetGlobalVerbosity

# Redefine ScalarFormatter


class ScalarFormatterForceFormat(ScalarFormatter):
    # Override function that finds format to use.
    def _set_format(self, vmin, vmax):
        self.format = "%1.2f"  # Give format here

# Action class to allow comma-separated list in options


class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))


# Some global definitions
_formats = {'eps': 0, 'png': 1, 'svg': 2}
_debug = False

########################################################################################################################


def main():
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
              'axes.labelsize':  'x-large',
              'axes.titlesize':  'x-large',
              #'axes.linewidth':  2, #increase default linewidth
              'xtick.labelsize': 'x-large',
              'ytick.labelsize': 'x-large'}
    mpl.rcParams.update(params)

    # Define arguments & options
    parser = argparse.ArgumentParser(
        epilog='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Positional arguments
    #    parser.add_argument('file', type=str,
    #                       help='Blafasel.')
    # Optional arguments
    parser.add_argument('-a', '--afile', default=None, required=True,
                        help='1st fastNLO result.')
    parser.add_argument('-b', '--bfile', default=None, required=True,
                        help='2nd fastNLO result.')
    parser.add_argument('-o', '--outfiles', default='fnlo-compare', required=False, type=str,
                        help='Start of output filename for plots. ' 'Default: fnlo-compare')

    # parse arguments
    args = vars(parser.parse_args())
    namesp = parser.parse_args()

    # Input files
    afile = args['afile']+'.log'
    bfile = args['bfile']+'.log'
    astat = args['afile']+'.dat'
    bstat = args['bfile']+'.dat'

    # Output files
    outfil = args['outfiles']

    exit(0)

    # Extract no. of columns from files
    print 'Reading no. of columns from 1st file ', afile
    cols = np.loadtxt(afile, comments=['#', ' #', 'C', 'L'])
    nacol = cols.shape[1]
    print 'Found ', nacol, ' columns'
    print 'Reading no. of columns from 2nd file ', bfile
    cols = np.loadtxt(bfile, comments=['#', ' #', 'C', 'L'])
    nbcol = cols.shape[1]
    print 'Found ', nbcol, ' columns'
    if nacol != nbcol:
        print 'Different no. of columns found, aborted!'
        exit(1)
        ncol = nacol

    # Read files columnwise
    if ncol == 9:
        iobsa, bwa, id0a, xla, xua, xavea, xsa_lo, xsa_nlo, ka_nlo = np.loadtxt(afile, comments=['#', ' #', 'C', 'L'], dtype=[(
            'f0', int), ('f1', float), ('f2', int), ('f3', float), ('f4', float), ('f5', float), ('f6', float), ('f7', float), ('f8', float)], unpack=True)
        iobsb, bwb, id0b, xlb, xub, xaveb, xsb_lo, xsb_nlo, kb_nlo = np.loadtxt(bfile, comments=['#', ' #', 'C', 'L'], dtype=[(
            'f0', int), ('f1', float), ('f2', int), ('f3', float), ('f4', float), ('f5', float), ('f6', float), ('f7', float), ('f8', float)], unpack=True)
    else:
        print 'ncol different from 9 not yet implemented. Aborted!'
        exit(2)

    exit(0)


# Split columns according to no. of observable evaluations
nvala = iobsa.size
nvalb = iobsb.size
nobsa = max(iobsa)
nobsb = max(iobsb)
print 'Found ', nvala, ' result lines for ', nobsa, ' observable bins in file ', afile
print 'Found ', nvalb, ' result lines for ', nobsb, ' observable bins in file ', bfile
if nvala != nvalb or nobsa != nobsb:
    print 'nval or nobs different between the two files to compare, aborted!'
    exit(3)
nobs = nobsa
nval = nvala
nrep = nval/nobs
xl = np.split(xla, nrep)
xu = np.split(xla, nrep)
xsa_lo = np.split(xsa_lo, nrep)
xsb_lo = np.split(xsb_lo, nrep)
xsa_nlo = np.split(xsa_nlo, nrep)
xsb_nlo = np.split(xsb_nlo, nrep)
ka_nlo = np.split(ka_nlo, nrep)
kb_nlo = np.split(kb_nlo, nrep)

xsd = []
dstd = []
xsd.append(np.loadtxt(astat, usecols=(2,), comments=['#']))
dstd.append(np.loadtxt(astat, usecols=(3,), comments=['#']))
xsa = []
dxsa = []
xsa.append(np.loadtxt(astat, usecols=(4,), comments=['#']))
dxsa.append(np.loadtxt(astat, usecols=(5,), comments=['#']))
dsta = np.divide(dxsa, xsa, out=np.ones_like(dxsa), where=xsa != 0)

# Calculate ratios
r_ablo = np.divide(xsa_lo, xsb_lo, out=np.ones_like(xsa_lo), where=xsb_lo != 0)
r_abnlo = np.divide(xsa_nlo, xsb_nlo, out=np.ones_like(
    xsa_nlo), where=xsb_nlo != 0)

# Plots
titwgt = 'bold'
limfs = 'x-large'
for i in range(nrep):

    # Calculate statistical uncertainties
    dr_ablo = np.multiply(r_ablo[i], dsta)
    dr_abnlo = np.multiply(r_abnlo[i], dsta)
    r_data = np.divide(xsd, xsb_nlo[i], out=np.ones_like(
        xsd), where=xsb_nlo[i] != 0)
    dr_data = np.multiply(r_data, dstd)

    fig = plt.figure()
    ax = fig.gca()
    plt.title(r'Ratio a over b for scale var {}'.format(i), fontweight=titwgt)
    plt.xlabel('Observable bin index', horizontalalignment='right',
               x=1.0, verticalalignment='top', y=1.0)
    plt.ylabel('Ratio', horizontalalignment='right', x=1.0,
               verticalalignment='top', y=1.0, labelpad=20)
    plt.axhline(y=1.01, linestyle='--', linewidth=1.0, color='black')
    plt.axhline(y=0.99, linestyle='--', linewidth=1.0, color='black')
    plt.fill_between([0.0, nobs+1], 0.99, 1.01, color='black', alpha=0.1)
#    plt.text(nobs+1.5,1.0,u'$\pm$1‰',fontsize=limfs)
    plt.text(nobs+1.5, 1.0, u'$\pm$1%', fontsize=limfs)

    dx = 1./4.
    x = np.arange(1, nobs+1.e-6)
    xa = np.arange(1-dx, nobs-dx+1.e-6)
    xb = np.arange(1+dx, nobs+dx+1.e-6)

    # Don't show statistical uncertainties for NNLOJET grid/NNLOJET with identical events
    rlo = plt.errorbar(x, r_ablo[i], yerr=dr_ablo[0], marker='<',
                       linestyle='none', label=u'LO$\,\pm\,$stat. unc.', color='blue')
    rnlo = plt.errorbar(x, r_abnlo[i], yerr=dr_abnlo[0], marker='s',
                        linestyle='none', label=u'NLO$\,\pm\,$stat. unc.', color='green')
    rdat = plt.errorbar(x, r_data[0], yerr=dr_data[0]/100., marker='o',
                        linestyle='none', label=u'data$\,\pm\,$stat. unc.', color='red')

    plt.xlim(0.0, nobs+1)
    plt.ylim(0.50, 1.50)

    handles = [rlo, rnlo, rdat]
    labels = [h.get_label() for h in handles]

    legend = ax.legend(handles, labels, loc='upper left',
                       numpoints=1, handlelength=0)
    legend = ax.legend(handles, labels, title=r'{}'.format(
        outfil), loc='upper left', numpoints=1, handlelength=0)
    legend.get_title().set_fontsize(limfs)

    fignam = outfil+'_scalevar-no-'+str(i)+'.png'
    plt.savefig(fignam)

exit(0)

###################################################################################

if __name__ == "__main__":
    main()
