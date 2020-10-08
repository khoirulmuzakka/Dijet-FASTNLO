#!/usr/bin/env python3
#-*- coding:utf-8 -*-
#
########################################################################
#
# Plot the statistical uncertainty of all channels
#
# Created by K. Rabbertz, 07.10.2020
#
########################################################################
#
# python2 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import argparse
import glob
import os
import re
import string
import sys
import timeit
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (FormatStrFormatter, LogFormatter, NullFormatter, ScalarFormatter, AutoMinorLocator, MultipleLocator)
from matplotlib import cm
# We do not want any interactive plotting! Figures are saved to files instead.
# This also avoids the ANNOYANCE of frequently missing Tkinter/tkinter (python2/3) GUI backends!
# To produce scalable graphics for publication use eps, pdf, or svg as file format.
# For this to work we try the Cairo backend, which can do all of these plus the raster format png.
# If this is not usable, we fall back to the Agg backend capable only of png for nice web plots.
#ngbackends = mpl.rcsetup.non_interactive_bk
#print('[fastnnlo_statunc]: Non GUI backends are: ', ngbackends)
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
#        print('[fastnnlo_statunc]: Can not use cairo backend :-(')
#        print('                   cairocffi or pycairo are required to be installed')
    else:
        if cairo.version_info < (1, 11, 0):
            # Introduced create_for_data for Py3.
            usecairo = False
#            print('[fastnnlo_statunc]: Can not use cairo backend :-(')
#            print('                   cairo {} is installed; cairo>=1.11.0 is required'.format(cairo.version))
if usecairo:
    mpl.use('cairo')
else:
    backend = 'agg'
    useagg = True
    try:
        mpl.use(backend, force=True)
        print('[fastnnlo_statunc]: Warning! Could not import cairo backend :-( Using agg instead for raster plots only!')
    except:
        useagg = False
        print('[fastnnlo_statunc]: Can not use agg backend :-(')
        raise ImportError('[fastnnlo_statunc]: Neither cairo nor agg backend found :-( Cannot produce any plots. Good bye!')
    mpl.use('agg')
import matplotlib.pyplot as plt
# numpy
import numpy as np
# fastNLO for direct evaluation of interpolation grids
# ATTENTION: fastNLO python extension is required for Python 3!
import fastnlo
from fastnlo import fastNLOLHAPDF
from fastnlo import SetGlobalVerbosity
#import warnings
#warnings.filterwarnings("error")


# Redefine ScalarFormatter
class ScalarFormatterForceFormat(ScalarFormatter):
    # Override function that finds format to use.
    def _set_format(self, vmin, vmax):
        self.format = "%1.2f"  # Give format here


# Action class to allow comma-separated (or empty) list in options
class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values:
            setattr(namespace, self.dest, values[0].split(','))
        else:
            setattr(namespace, self.dest, [''])


# Some global definitions
_channels = [ 'LO', 'R', 'V', 'RRa', 'RRb', 'RV', 'VV' ]
_orders = [ 'LO', 'NLO', 'NNLO', 'NLO_only', 'NNLO_only' ]
_fntrans = str.maketrans({'[': '', ']': '', '(': '', ')': '', ',': ''}) # Filename translation table
_formats = {'eps': 0, 'pdf': 1, 'png': 2, 'svg': 3}
_text_to_order = {'LO': 0, 'NLO': 1, 'NNLO': 2, 'NLO_only': 3, 'NNLO_only': 4}
_order_to_text = {0: 'LO', 1: 'NLO', 2: 'NNLO', 3: 'NLO_only', 4: 'NNLO_only'}
_order_color = {'LO': 'g', 'NLO': 'b', 'NNLO': 'r', 'NLO_only': 'tab:olive', 'NNLO_only': 'tab:brown'}
_order_symbol = {'LO': '^', 'NLO': 's', 'NNLO': 'o', 'NLO_only': 'tri_up', 'NNLO_only': 'plus'}
_channel_color = {'LO': 'tab:green', 'R': 'tab:cyan', 'V': 'tab:blue', 'RRa': 'tab:red', 'RRb': 'tab:orange', 'RV': 'tab:pink', 'VV': 'tab:purple'}
_channel_symbol = {'LO': 'o', 'R': 's', 'V': 'D', 'RRa': 'v', 'RRb': '^', 'RV': '<', 'VV': '>'}
_colors = ['tab:orange', 'tab:green', 'tab:purple', 'tab:blue', 'tab:brown']
_symbols = ['s', 'X', 'o', '^', 'v']
_hatches = ['', '//', '\\', '|', '-']
_scale_to_text = {0: 'kScale1', 1: 'kScale2', 2: 'kQuadraticSum', 3: 'kQuadraticMean', 4: 'kQuadraticSumOver4',
                  5: 'kLinearMean', 6: 'kLinearSum', 7: 'kScaleMax', 8: 'kScaleMin', 9: 'kProd',
                  10: 'kS2plusS1half', 11: 'kPow4Sum', 12: 'kWgtAvg', 13: 'kS2plusS1fourth', 14: 'kExpProd2', 15: 'kExtern'}
_pdfbasenames = ['ABM11', 'ABMP15', 'ABMP16', 'CJ12', 'CJ15', 'CT10', 'CT14', 'HERAPDF20', 'JR14',
                 'MMHT2014', 'MSTW2008', 'NNPDF23', 'NNPDF30', 'NNPDF31']
_debug = False

#####################################################################################


#
# Function plotting statistical uncertainties for each channel with respect to
# the total at a chosen order
#
def plotting(x_axis, xmin, xmax, dxsr_ch, dxsr_or, xlabel, ylabel, title, tablename, order, given_filename, scale_name, nice_scale_name, formats):

    fig = plt.figure(figsize=(7, 7))
    gs = gridspec.GridSpec(1, 1, fig)
    ax1 = plt.subplot(gs[0, :])
    ax1.set_autoscalex_on(True)
    plt.setp(ax1.get_xticklabels(), visible=True)

    # For plotting shifted results, 'next to each other', handling via shift from bincenter
    # Always seven _channels
    shift_list = [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]

    ax1.set_xscale('log', nonposx='clip')
    ax1.set_yscale('log', nonposy='clip')
    ax1.set_ylim(0.00001, 0.1)
    ax1.set_xlabel(r'%s' %xlabel, horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
    ax1.set_ylabel(r'%s(%s)' % (ylabel, order), horizontalalignment='right', x=1.0,
                   verticalalignment='top', y=1.0, rotation=90, labelpad=24)
    ax1.set_title('%s' % title, loc='left')
#    ax1.text(0.03, 0.15, 'Normalised to: %s' % order, horizontalalignment='left',
#             verticalalignment='bottom', transform=ax1.transAxes)

    # Plot all normalised statistical uncertainties in channel list
    xs_index = -1
    for channel_item, shift in zip(_channels, shift_list):
        xs_index += 1
        yerror = 0*dxsr_ch[xs_index, :]
        ax1.errorbar(x_axis*shift, dxsr_ch[xs_index], yerr=yerror, elinewidth=1, linewidth=0.0,
                     ms=6, marker=_channel_symbol[channel_item], color=_channel_color[channel_item], fmt='.', label=channel_item)
    ax1.errorbar(x_axis*shift, dxsr_or[_text_to_order[order]], yerr=yerror, elinewidth=1, linewidth=0.0,
                 ms=6, marker='X', color='k', fmt='.', label='total')
    ax1.legend(fontsize=10, numpoints=1)

    fig.tight_layout()

    if given_filename is not None:
        filename = '%s.statunc-%s.%s' % (given_filename, order, scale_name)
    else:
        filename = '%s.statunc-%s.%s' % (tablename, order, scale_name)
        filename = filename+'.stat'

    # Do not use characters defined in _fntrans for filenames
    filename = filename.translate(_fntrans)

    for fmt in formats:
        figname = '%s.%s' % (filename, fmt)
        fig.savefig(figname)
        print('[fastnnlo_statunc]: Plot saved as:', figname)
#        print(plt.get_fignums())
    plt.close(fig)


########################################################################

def main():
    # Start timer
    # just for measuring wall clock time - not necessary
    start_time = timeit.default_timer()
    # Define arguments & options
    parser = argparse.ArgumentParser(
        epilog='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Positional arguments
    parser.add_argument('table', type=str, nargs='+',
                        help='Filename glob of fastNLO tables to be evaluated. This must be specified!')
    # Optional arguments
    parser.add_argument('-d', '--datfiles', required=False, nargs='?', type=str, action=SplitArgs,
                        help='Comma-separated or empty list of NNLOJET dat files with statistical uncertainties for each channel. If empty, dat files matching to table name are used.')
    parser.add_argument('-f', '--filename', default=None, type=str,
                        help='Output filename (optional).')
    parser.add_argument('--format', required=False, nargs=1, type=str, action=SplitArgs,
                        help='Comma-separated list of plot formats to use: eps, pdf, png, svg. If nothing is chosen, png is used.')
    parser.add_argument('-o', '--order', default='LO', type=str,
                        help='Order to normalise to: LO, NLO, or NNLO. If nothing is chosen, LO is used.')
    parser.add_argument('-s', '--scale', default=0, required=False, nargs='?', type=int,
                        choices=list(range(16)), metavar='[0-15]',
                        help='For flexible-scale tables define central scale choice for MuR and MuF by selection enum fastNLO::ScaleFunctionalForm ("0"=kScale1, "1"=kScale2, "2"=kQuadraticSum), ...')
    parser.add_argument('--scalename', default=None, type=str,
                        help='Replace default scale name by given string.')
    parser.add_argument('--title', default=None, type=str,
                        help='Replace table name as default title by given string.')
    parser.add_argument('-v', '--verbose', action="store_true",
                        help="Increase output verbosity.")
    parser.add_argument('--xlabel', default=None, type=str,
                        help='Replace x axis default label by given string.')
    parser.add_argument('--ylabel', default=None, type=str,
                        help='Replace y axis default label by given string.')

    # Print header
    print("\n###########################################################################################")
    print("# fastnnlo_statunc:")
    print("# Plot the statistical uncertainty composition")
    print("###########################################################################################\n")

    # Parse arguments
    args = vars(parser.parse_args())

    # List of table names
    files = args['table']
    print('\n')
    print('[fastnnlo_statunc]: Analysing table list:')
    for file in files:
        print('[fastnnlo_statunc]:   ', file)

    # Order to normalise to
    order = args['order']
    if order not in _orders:
        print('[fastnnlo_statunc]: Illegal order chosen, aborted! order =', order)
        exit(1)
    print('[fastnnlo_statunc]: Normalising statistical uncertainties to', order)

    # Check existence of NNLOJET dat files for all channels and normalising order
    datfilenames = []
    if args['datfiles'] is None or args['datfiles'][0] == '':
        print('[fastnnlo_statunc]: Automatic filename matching is used to load statistical uncertainties from NNLOJET.')
    else:
        for datfile in args['datfiles']:
            lstat = os.path.isfile(datfile)
            if lstat:
                datfilenames.append(datfile)
            else:
                print('[fastnnlo_statunc]: Given file ', datfile,
                      'for statistical uncertainties not found, aborted!')
                exit(1)

    # Scale choice
    scale_choice = args['scale']

    # Scale name
    nice_scale_name = args['scalename']

    # Given filename
    given_filename = args['filename']

    # Plot formats to use
    formats = args['format']
    if formats is None:
        formats = ['png']
    for fmt in formats:
        if fmt not in _formats:
            print('[fastnnlo_statunc]: Illegal format specified, aborted!')
            print('[fastnnlo_statunc]: Format list:', args['format'])
            exit(1)
        elif fmt != 'png' and not usecairo:
            print('[fastnnlo_statunc]: Vector format plots not possible without cairo backend, aborted!')
            print('[fastnnlo_statunc]: Format list:', args['format'])
            exit(1)

    # Plot labelling
    nice_title = args['title']
    nice_xlabel = args['xlabel']
    nice_ylabel = args['ylabel']

    # Verbosity
    verb = args['verbose']
    if verb:
        print('[fastnnlo_statunc]: Using matplotlib version ', mpl.__version__)
        print('                     from location ', mpl.__file__)




    # Loop over table list
    for table in files:
        print('[fastnnlo_statunc]: Analysing table: ', table)
        # Get rid of extensions (.tab.gz or .tab)
        tablename = os.path.splitext(os.path.basename(table))[0]
        tablename = os.path.splitext(tablename)[0]

        ###################### Start EVALUATION with fastNLO library ###################################################
        # SetGlobalVerbosity(0) # Does not work since changed to default in the following call
        # Only take the general information (bin_bounds, labels, order_existence, etc.) from table.
        fnlo = fastnlo.fastNLOLHAPDF(table, 'CT14nnlo', 0)

        # Get labeling for the x-axis
        # Dimensionality of the table:
        ndim = fnlo.GetNumDiffBin()
        if verb:
            print('\n')
            print('[fastnnlo_statunc]: Table Dimensions: ', ndim)

        # Labels of all the dimensions:
        labels = fnlo.GetDimLabels()
        if verb:
            print('[fastnnlo_statunc]: Labels:', labels)

        # x label of first dimension from table:
        xlabel = fnlo.GetDimLabel(0)
        if verb:
            print('[fastnnlo_statunc]: x-label:', xlabel)

        # Generic y label
        ylabel = '$\Delta\sigma_\mathrm{stat}\,/\,\sigma$'

        # Creating x-axis
        bin_bounds = np.array(fnlo.GetObsBinsBounds(0))
        if verb:
            print('[fastnnlo_statunc]: bin_bounds.T: \n', bin_bounds.T, '\n')
            print('[fastnnlo_statunc]: bin_bounds.flatten()',
                  bin_bounds.flatten(), '\n')

        x_axis = (bin_bounds.T[0]+bin_bounds.T[1]) / 2.  # this is a list of bin centers
        xmin = 0.95*min(bin_bounds.ravel())
        xmax = 1.05*max(bin_bounds.ravel())
        if verb:
            print('[fastnnlo_statunc]: xmin=%s, xmax=%s. \n' % (xmin, xmax))

        # Preparing x-errors (via bin_bounds) --> x_errors[0, :] are initially negative (positive via -1*), x_errors[1, :] positive
        x_errors = np.array(
            [-1*(bin_bounds.T[0]-x_axis), bin_bounds.T[1]-x_axis])
        if verb:
            print('[fastnnlo_statunc]: \n x_errors: ', x_errors, '\n')

        # Check existence of orders in table
        lflex = fnlo.GetIsFlexibleScaleTable()
        scale_name = 'kScale1'
        o_existence = [False, False, False]
        cnt_order = -1
        max_order = 0
        for i in range(3):
            o_existence[i] = fnlo.SetContributionON(
                fastnlo.kFixedOrder, i, True)
            if o_existence[i]:
                max_order = i
                if not lflex:
                    if scale_choice != 0:
                        print(
                            '[fastnnlo_statunc]: Invalid choice of scale = ', scale_choice, ' Aborted!')
                        print(
                            '[fastnnlo_statunc]: For fixed-scale tables only the default=0 is allowed.')
                        exit(1)
                    else:
                        scale_name = fnlo.GetScaleDescription(i, 0)
                else:
                    if scale_choice < 2:
                        scale_name = fnlo.GetScaleDescription(i, scale_choice)
                    else:
                        scl0 = fnlo.GetScaleDescription(i, 0)
                        scl1 = fnlo.GetScaleDescription(i, 1)
                        scale_name = _scale_to_text[scale_choice] + \
                            '_'+scl0+'_'+scl1
                if cnt_order == i-1:
                    cnt_order += 1
        if verb:
            print('[fastnnlo_statunc]: Table has continuous orders up to',
                  cnt_order, 'and a maximal order of', max_order)

        # Set iordmax to maximum found in table
        iordmax = max_order

        if iordmax > cnt_order:
            print('[fastnnlo_statunc]: Invalid choice of orders. Aborted!')
            print('[fastnnlo_statunc]: Highest order requested is',
                  _order_to_text[iordmax], 'but continuous orders are available only up to', _order_to_text[cnt_order])
            exit(1)

        # Read in normalisation cross section (LO, NLO, or NNLO)
        sep = '.'
        parts = tablename.split(sep)
        parts[1] = order
        datfile = sep.join(parts) + '.dat'
        print('[fastnnlo_statunc]: Taking normalisation cross section from', datfile)
        cols = np.loadtxt(datfile, usecols=list(range(3, 5)))
        xs_norm = np.array(cols[:, 0])

        # Read in cross section and statistical uncertainty for each order and exclusive order
        xs      = []
        dxs     = []
        dxsr_or = []
        for lorder in _orders:
            parts = tablename.split(sep)
            parts[1] = lorder
            datfile = sep.join(parts) + '.dat'
            print('[fastnnlo_statunc]: Reading statistical uncertainties from', datfile)
            cols = np.loadtxt(datfile, usecols=list(range(3, 5)))
            xs_dat   = np.array(cols[:, 0])
            dxs_dat  = np.array(cols[:, 1])
            dxsr_dat = np.divide(dxs_dat, xs_norm, out=np.ones_like(dxs_dat), where=xs_norm != 0)
            xs.append(xs_dat)
            dxs.append(dxs_dat)
            dxsr_or.append(dxsr_dat)

        # Read in statistical uncertainty for each channel, order, and exclusive order
        dxsr_ch = []
        if args['datfiles'] is None or args['datfiles'][0] == '':
            for channel in _channels:
                parts = tablename.split(sep)
                parts[1] = channel
                datfile = sep.join(parts) + '.dat'
                datfilenames.append(datfile)

        lstat = (len(datfilenames) > 0)
        if lstat and len(datfilenames) != len(_channels):
            print('[fastnnlo_statunc]: Mismatch between required channels and no. of filenames for statistical uncertainties, aborted!')
            exit(1)

        dxsr = []
        for fname in datfilenames:
            print('[fastnnlo_statunc]: Taking statistical uncertainties from', fname)
            cols = np.loadtxt(fname, usecols=list(range(3, 5)))
            dxs_dat = np.array(cols[:, 1])
            dxsr_dat = np.divide(dxs_dat, xs_norm, out=np.ones_like(dxs_dat), where=xs_norm != 0)
            dxsr.append(dxsr_dat)
        dxsr_ch = abs(np.array(dxsr))
        # Empty list for use with next table in automatic mode
        if args['datfiles'] is None or args['datfiles'][0] == '':
            datfilenames = []



        ############################## Do the plotting ####################################################

        if nice_title is None:
            title = tablename
        else:
            title = nice_title
        if nice_scale_name is None:
            nice_scale_name = scale_name
        if nice_xlabel is not None:
            xlabel = nice_xlabel
        if nice_ylabel is not None:
            ylabel = nice_ylabel

        # Without statistical uncertainties create zero array of proper dimensions here
        if _debug:
            print('dxsr_ch', dxsr_ch)

        #        plotting(x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, dxsr_ch, nostat, xlabel, ylabel, title, tablename,
        #                 order_list, given_filename, scale_name, nice_scale_name, variation_type, formats)
        plotting(x_axis, xmin, xmax, dxsr_ch, dxsr_or, xlabel, ylabel, title, tablename,
                 order, given_filename, scale_name, nice_scale_name, formats)

        stop_time = timeit.default_timer()
        timediff = stop_time-start_time
        print('fastnnlo_statunc: Elapsed time: %s sec = %s min' %
              (timediff, round(timediff/60, 2)))


if __name__ == '__main__':
    main()
