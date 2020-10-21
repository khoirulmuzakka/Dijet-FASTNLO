#!/usr/bin/env python3
#-*- coding:utf-8 -*-
#
########################################################################
#
# Plot the scale uncertainty
#
# Created by B. Schillinger, 09.10.2018
# Modified by K. Rabbertz, 31.10.2018
# Prepared for python3 by K. Rabbertz, 07.03.2020
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
#print('[fastnnlo_scaleunc]: Non GUI backends are: ', ngbackends)
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
#        print('[fastnnlo_scaleunc]: Can not use cairo backend :-(')
#        print('                   cairocffi or pycairo are required to be installed')
    else:
        if cairo.version_info < (1, 11, 0):
            # Introduced create_for_data for Py3.
            usecairo = False
#            print('[fastnnlo_scaleunc]: Can not use cairo backend :-(')
#            print('                   cairo {} is installed; cairo>=1.11.0 is required'.format(cairo.version))
if usecairo:
    mpl.use('cairo')
else:
    backend = 'agg'
    useagg = True
    try:
        mpl.use(backend, force=True)
        print('[fastnnlo_scaleunc]: Warning! Could not import cairo backend :-( Using agg instead for raster plots only!')
    except:
        useagg = False
        print('[fastnnlo_scaleunc]: Can not use agg backend :-(')
        raise ImportError('[fastnnlo_scaleunc]: Neither cairo nor agg backend found :-( Cannot produce any plots. Good bye!')
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
_fntrans = str.maketrans({'[': '', ']': '', '(': '', ')': '', ',': ''}) # Filename translation table
_formats = {'eps': 0, 'pdf': 1, 'png': 2, 'svg': 3}
_text_to_order = {'LO': 0, 'NLO': 1, 'NNLO': 2}
_order_to_text = {0: 'LO', 1: 'NLO', 2: 'NNLO'}
_order_color = {'LO': 'g', 'NLO': 'b', 'NNLO': 'r'}
_order_symbol = {'LO': '^', 'NLO': 's', 'NNLO': 'o'}
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


# Function plotting cross sections for a list of orders into one figure including scale uncertainties
# The ratio is done always with respect to the first order appearing in the order_list
# given e.g. via '-o NLO,NNLO', i.e. the first evaluated cross section, here NLO,
# that is stored in xs_all[0].
def plotting(x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, dxsr_cn, nostat, xlabel, ylabel, title, tablename, order_list, given_filename, scale_name, nice_scale_name, pdfset, variation_type, formats):
    if variation_type == 'Scale uncertainty (2P)':
        vartype = '2P'
    elif variation_type == 'Scale uncertainty (6P)':
        vartype = '6P'

    pdfnicename = 'Undefined'
    for pdfn in _pdfbasenames:
        if pdfn in pdfset:
            pdfnicename = pdfn

    gs = gridspec.GridSpec(3, 3, hspace=0)
    fig = plt.figure(figsize=(7, 7))
    ax1 = plt.subplot(gs[:-1, :])
    ax1.set_autoscalex_on(False)
    plt.setp(ax1.get_xticklabels(), visible=False)

    # For plotting various results, max = 3, 'next to each other', handling via shift from bincenter
    if len(order_list) == 1:
        shift_list = [1.00]
    elif len(order_list) == 2:
        shift_list = [0.98, 1.02]
    elif len(order_list) == 3:
        shift_list = [0.96, 1.00, 1.04]
    else:
        print('[fastnnlo_scaleunc]: Too many orders to plot simultaneously. Aborted!')
        print('[fastnnlo_scaleunc]: Current maximum is 3. The order list is', order_list)
        sys.exit(1)

    # Upper subplot setup
    #
    # Set limits on x axis on coupled axis from lower plot
    #    ax1.set_xlim(left=xmin, right=xmax)
    # TODO: 'minor_thresholds' needs matplotlib > 1.5.0
    #    axfmt = LogFormatter(labelOnlyBase=False, minor_thresholds=(2, 0.4))
    #    ax1.get_xaxis().set_minor_formatter(axfmt)
    #        ax1.get_xaxis().set_minor_formatter(NullFormatter())
    ax1.set_xscale('log', nonposx='clip')
    ax1.set_yscale('log', nonposy='clip')
    # Set label on x axis on coupled axis from lower plot
    #        ax1.set_xlabel(r'%s' %xlabel, horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
    ax1.set_ylabel(r'%s' % ylabel, horizontalalignment='right', x=1.0,
                   verticalalignment='top', y=1.0, rotation=90, labelpad=24)
    # TODO: Commented out for now, since loc attribute not defined in prehistoric Centos7 matplotlib version 1.2.0
    ax1.set_title('%s' % title, loc='left')
    #    ax1.set_title('%s' % title)
    ax1.text(0.03, 0.15, 'PDF set: %s' % pdfnicename, horizontalalignment='left',
             verticalalignment='bottom', transform=ax1.transAxes)
    ax1.text(0.03, 0.08, 'Scale: %s' % nice_scale_name, horizontalalignment='left',
             verticalalignment='bottom', transform=ax1.transAxes)
    ax1.text(0.03, 0.03, '%s' % variation_type, horizontalalignment='left',
             verticalalignment='bottom', transform=ax1.transAxes)

    #        Only for publication
    # H1
    #    ax1.text(0.35, 0.90, r'$30 < Q^2 < 42\,\mathrm{GeV}^2$',
    #             horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)
    # ZEUS
    #    ax1.text(0.35, 0.90, r'$500 < Q^2 < 1000\,\mathrm{GeV}^2$',
    #             horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)

    # Plot all x section results in order_list
    xs_index = -1
    for order_item, shift in zip(order_list, shift_list):
        xs_index += 1
        yerror = 0*xs_all[xs_index, :]
        if not nostat:
            yerror = np.multiply(
                xs_all[xs_index, :], dxsr_cn[xs_index, :])
        else:
            yerror = abs(abs_scale_unc[xs_index])
        ax1.errorbar(x_axis*shift, xs_all[xs_index], yerr=yerror, elinewidth=1, linewidth=0.0,
                     ms=6, marker=_order_symbol[order_item], color=_order_color[order_item], fmt='.', label=order_item)
        ax1.fill_between(x_axis*shift, xs_all[xs_index] + xs_all[xs_index]*rel_scale_unc[xs_index, 2, :],
                         xs_all[xs_index] + xs_all[xs_index]*rel_scale_unc[xs_index, 1, :], color=_order_color[order_item], hatch=_hatches[xs_index], alpha=0.30)
    ax1.legend(fontsize=10, numpoints=1)


    # Lower subplot setup
    #
    ax2 = plt.subplot(gs[2, :], sharex=ax1)
    # Set common x axis bounds for both
    ax2.set_xlim(left=xmin, right=xmax)
    #    ax2.set_xbound(lower=xmin, upper=xmax)
    ax2.set_yscale('linear', nonposy='clip')
    ax2.set_xlabel(r'%s' % xlabel, horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
    ax2.set_ylabel(r'Ratio to %s' % order_list[0], horizontalalignment='center', x=1.0, verticalalignment='top', y=0.5, rotation=90, labelpad=24)
    #        ax2.legend(fontsize=10, numpoints=1)
    #        ax1.set_xticklabels([])
    #        ax1.set_xticks([])
    #        ax1.get_xaxis().set_visible(False)
    ax2.axhline(y=1, xmin=0, xmax=1, color='k', linestyle='dotted', linewidth=1.6, alpha=0.2)

    # Ratio subplot with relative scale uncertainties; denominator in ratio = first order in order_list
    xs_index = -1
    ordernames = ''
    for item in order_list:
        xs_index += 1
        # Do not show error bars
        yerror = 0*xs_all[xs_index, :]
        mymarker = _order_symbol[item]
        mymksize = 4
        if nostat:
            # Error bars = scale uncertainty
            #            yerror = abs(abs_scale_unc[xs_index])
            #            mymksize = 4
            # No error bars, only small symbols
            yerror = 1*yerror
            #            mymarker = '.'
        else:
            # Error bars = statistical uncertainty if known
            yerror = np.multiply(
                xs_all[xs_index, :], dxsr_cn[xs_index, :])
        ordernames += '_%s' % item
        ax2.errorbar(x_axis, xs_all[xs_index]/xs_all[0], yerr=yerror/xs_all[0], barsabove=True, elinewidth=1, linewidth=0.0,
                     ms=mymksize, marker=mymarker, color=_order_color[item], fmt='.', label=item)
        ax2.fill_between(x_axis, (xs_all[xs_index]*(1+rel_scale_unc[xs_index, 2, :])/xs_all[0]),
                         (xs_all[xs_index]*(1+rel_scale_unc[xs_index, 1, :])/xs_all[0]), color=_order_color[item], hatch=_hatches[xs_index], alpha=0.30)

    fig.tight_layout()

    if given_filename is not None:
        filename = '%s.scaleunc-%s.%s.%s' % (given_filename,
                                             vartype, ordernames[1:], scale_name)
    else:
        filename = '%s.scaleunc-%s.%s.%s.%s' % (
            tablename, vartype, ordernames[1:], pdfset, scale_name)
        if not nostat:
            filename = filename+'.stat'

    # Do not use characters defined in _fntrans for filenames
    filename = filename.translate(_fntrans)

    for fmt in formats:
        figname = '%s.%s' % (filename, fmt)
        fig.savefig(figname)
        print('[fastnnlo_scaleunc]: Plot saved as:', figname)
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
    parser.add_argument('-a', '--asymmetric', action="store_true",
                        help='If -a is chosen, use asymmetric (6P) scale variations; otherwise use symmetric ones (2P).')
    parser.add_argument('-d', '--datfiles', required=False, nargs='?', type=str, action=SplitArgs,
                        help='Comma-separated or empty list of NNLOJET dat files with statistical uncertainties for each order to show. If empty, dat files matching to table name are used. If not set, statistical uncertainties are ignored.')
    parser.add_argument('-f', '--filename', default=None, type=str,
                        help='Output filename (optional).')
    parser.add_argument('--format', required=False, nargs=1, type=str, action=SplitArgs,
                        help='Comma-separated list of plot formats to use: eps, pdf, png, svg. If nothing is chosen, png is used.')
    parser.add_argument('-m', '--member', default=0, type=int,
                        help='Member of PDFset, default is 0.')
    parser.add_argument('-o', '--order', required=False, nargs=1, type=str, action=SplitArgs,
                        help='Comma-separated list of orders to show: LO, NLO, and/or NNLO. If nothing is chosen, show all orders available in table.')
    parser.add_argument('-p', '--pdfset', default='CT14nlo', type=str,
                        help='PDFset to evaluate fastNLO table.')
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
    print("# fastnnlo_scaleunc:")
    print("# Plot the scale uncertainty")
    print("###########################################################################################\n")

    # Parse arguments
    args = vars(parser.parse_args())

    # List of table names
    files = args['table']
    print('\n')
    print('[fastnnlo_scaleunc]: Analysing table list:')
    for file in files:
        print('[fastnnlo_scaleunc]:   ', file)

    # PDF set name
    pdfset = os.path.basename(args['pdfset'])
    print('[fastnnlo_scaleunc]: Using PDF set', pdfset)

    # Orders to be shown
    iorders = []
    iordmin = _text_to_order['LO']
    iordmax = _text_to_order['NNLO']
    if args['order'] is None:
        print('[fastnnlo_scaleunc]: Evaluate table up to highest available order.')
    else:
        for ord in args['order']:
            if ord in _text_to_order:
                iorders.append(_text_to_order[ord])
            else:
                print('[fastnnlo_scaleunc]: Illegal order specified, aborted!')
                print('[fastnnlo_scaleunc]: Order list:', args['order'])
                exit(1)
        iordmin = min(iorders)
        iordmax = max(iorders)
        print('[fastnnlo_scaleunc]: Evaluate table up to order(s)', args['order'])

    # Check existence of NNLOJET dat file for each order if desired
    datfilenames = []
    nostat = False
    if args['datfiles'] is None:
        nostat = True
        print('[fastnnlo_scaleunc]: No statistical uncertainties requested.')
    elif args['datfiles'][0] == '':
        print('[fastnnlo_scaleunc]: Automatic filename matching is used to load statistical uncertainties from NNLOJET.')
    else:
        for datfile in args['datfiles']:
            lstat = os.path.isfile(datfile)
            if lstat:
                datfilenames.append(datfile)
            else:
                print('[fastnnlo_scaleunc]: Given file ', datfile,
                      'for statistical uncertainties not found, aborted!')
                exit(1)

    # Type of scale variation (symmetric vs asymmetric)
    if args['asymmetric']:
        scale_var_type = fastnlo.kAsymmetricSixPoint
        variation_type = 'Scale uncertainty (6P)'
    else:
        scale_var_type = fastnlo.kSymmetricTwoPoint
        variation_type = 'Scale uncertainty (2P)'

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
            print('[fastnnlo_scaleunc]: Illegal format specified, aborted!')
            print('[fastnnlo_scaleunc]: Format list:', args['format'])
            exit(1)
        elif fmt != 'png' and not usecairo:
            print('[fastnnlo_scaleunc]: Vector format plots not possible without cairo backend, aborted!')
            print('[fastnnlo_scaleunc]: Format list:', args['format'])
            exit(1)

    # Plot labelling
    nice_title = args['title']
    nice_xlabel = args['xlabel']
    nice_ylabel = args['ylabel']

    # Verbosity
    verb = args['verbose']
    if verb:
        print('[fastnnlo_scaleunc]: Using matplotlib version ', mpl.__version__)
        print('                     from location ', mpl.__file__)

    # Loop over table list
    for table in files:
        # Table name
        tablepath = os.path.split(table)[0]
        if not tablepath:
            tablepath = '.'
        tablename = os.path.split(table)[1]
        if tablename.endswith('.tab.gz'):
            tablename = tablename.replace('.tab.gz', '', 1)
        elif tablename.endswith('.tab'):
            tablename = tablename.replace('.tab', '', 1)
        else:
            print('[fastnnlo_scaleunc]: Error! Wrong extension for table: ', table)
            exit(1)
        print('[fastnnlo_scaleunc]: Analysing table: ', table)

        ###################### Start EVALUATION with fastNLO library ###################################################
        # SetGlobalVerbosity(0) # Does not work since changed to default in the following call
        # Take the general information (bin_bounds, labels, order_existence, etc.) from given pdfset.
        fnlo = fastnlo.fastNLOLHAPDF(table, args['pdfset'], args['member'])

        # Get labeling for the x-axis
        # Dimensionality of the table:
        ndim = fnlo.GetNumDiffBin()
        if verb:
            print('\n')
            print('[fastnnlo_scaleunc]: Table Dimensions: ', ndim)

        # Labels of all the dimensions:
        labels = fnlo.GetDimLabels()
        if verb:
            print('[fastnnlo_scaleunc]: Labels:', labels)

        # x label of first dimension from table:
        xlabel = fnlo.GetDimLabel(0)
        if verb:
            print('[fastnnlo_scaleunc]: x-label:', xlabel)

        # Generic y label
        ylabel = '$\sigma \pm \Delta\sigma(\mu_R,\mu_F)$'

        # Creating x-axis
        bin_bounds = np.array(fnlo.GetObsBinsBounds(0))
        if verb:
            print('[fastnnlo_scaleunc]: bin_bounds.T: \n', bin_bounds.T, '\n')
            print('[fastnnlo_scaleunc]: bin_bounds.flatten()',
                  bin_bounds.flatten(), '\n')

        x_axis = (bin_bounds.T[0]+bin_bounds.T[1]) / 2.  # this is a list of bin centers
        xmin = 0.95*min(bin_bounds.ravel())
        xmax = 1.05*max(bin_bounds.ravel())
        if verb:
            print('[fastnnlo_scaleunc]: xmin=%s, xmax=%s. \n' % (xmin, xmax))

        # Preparing x-errors (via bin_bounds) --> x_errors[0, :] are initially negative (positive via -1*), x_errors[1, :] positive
        x_errors = np.array(
            [-1*(bin_bounds.T[0]-x_axis), bin_bounds.T[1]-x_axis])
        if verb:
            print('[fastnnlo_scaleunc]: \n x_errors: ', x_errors, '\n')

        # Check existence of orders in table
        lflex = fnlo.GetIsFlexibleScaleTable()
        scale_name = 'scale1'
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
                            '[fastnnlo_scaleunc]: Invalid choice of scale = ', scale_choice, ' Aborted!')
                        print(
                            '[fastnnlo_scaleunc]: For fixed-scale tables only the default=0 is allowed.')
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
            print('[fastnnlo_scaleunc]: Table has continuous orders up to',
                  cnt_order, 'and a maximal order of', max_order)

        # If previously undefined set iordmax to maximum found in table
        if args['order'] is None:
            iordmax = max_order

        if iordmax > cnt_order:
            print('[fastnnlo_scaleunc]: Invalid choice of orders. Aborted!')
            print('[fastnnlo_scaleunc]: Highest order requested is',
                  _order_to_text[iordmax], 'but continuous orders are available only up to', _order_to_text[cnt_order])
            exit(1)

        order_list = []
        if args['order'] is None:
            for iord in range(cnt_order+1):
                order_list.append(_order_to_text[iord])
        else:
            order_list = args['order']

        print('[fastnnlo_scaleunc]: List of requested orders:', order_list)

        # Read in statistical uncertainty for each order if requested
        dxsr_cn = []
        sep = '.'
        if not nostat:
            if args['datfiles'][0] == '':
                for order in order_list:
                    parts = tablename.split(sep)
                    parts[1] = order
                    datfile = tablepath + '/' + sep.join(parts) + '.dat'
                    datfilenames.append(datfile)

            lstat = (len(datfilenames) > 0)
            if lstat and len(datfilenames) != len(order_list):
                print('[fastnnlo_scaleunc]: Mismatch between no. of requested orders and no. of filenames for statistical uncertainties, aborted!')
                exit(1)

            dxsr = []
            for fname in datfilenames:
                print(
                    '[fastnnlo_scaleunc]: Taking statistical uncertainties from', fname)
                cols = np.loadtxt(fname, usecols=list(range(3, 5)))
                xs_dat = np.array(cols[:, 0])
                dxs_dat = np.array(cols[:, 1])
                dxsr_dat = np.divide(dxs_dat, xs_dat, out=np.ones_like(
                    dxs_dat), where=xs_dat != 0)
                dxsr.append(dxsr_dat)
            dxsr_cn = abs(np.array(dxsr))
            # Empty list for use with next table in automatic mode
            if args['datfiles'][0] == '':
                datfilenames = []

        # For flexible-scale tables set scale to user choice (default is 0)

        if lflex:
            print(
                '[fastnnlo_scaleunc]: Setting requested scale choice for flexible-scale table:', scale_choice)
            fnlo.SetMuRFunctionalForm(scale_choice)
            fnlo.SetMuFFunctionalForm(scale_choice)
        else:
            if scale_choice == 0:
                print(
                    '[fastnnlo_scaleunc]: Evaluating fixed-scale table. Scale choice must be', scale_choice)
            else:
                print(
                    '[fastnnlo_scaleunc]: No scale choice possible for fixed-scale table. Aborted!')
                print('[fastnnlo_scaleunc]: scale_choice = ', scale_choice)
                exit(1)

        # Now evaluate fastNLO table having a look at scale uncertainties
        xs_list = []  # will contain total cross section for selected orders out of LO, NLO, NNLO
        # list for relative scale uncertainties (low, high) for selected orders
        rel_unc_list = []
        for n in order_list:
            for j in range(0, max_order+1):
                if j <= _text_to_order[n]:
                    fnlo.SetContributionON(fastnlo.kFixedOrder, j, True)
                else:
                    fnlo.SetContributionON(fastnlo.kFixedOrder, j, False)
            if verb:
                print('[fastnnlo_scaleunc]: \n')
                print('[fastnnlo_scaleunc]: Calculate XS for order: %s' % n, '\n')
                print(
                    '[fastnnlo_scaleunc]: ----  ----  ----  ----  ----  ----  ----  ----')
            fnlo.CalcCrossSection()
            xs_list.append(fnlo.GetCrossSection())

            ### Get scale uncertainties ###
            if verb:
                print('[fastnnlo_scaleunc]: Used scale factor MuF: ',
                      fnlo.GetScaleFactorMuF())
                print('[fastnnlo_scaleunc]: Used scale factor MuR: ',
                      fnlo.GetScaleFactorMuR(), '\n')
                print('[fastnnlo_scaleunc]: Calculate scale uncertainties \n')
            # RELATIVE scale uncertainty with chosen type of scale variation (symmetric or asymmetric)
            # Up to NLO, it is possible to use HOPPET with fixed-scale tables
            #                fnlo.UseHoppetScaleVariations(True)
            # Calculate this already for all accessible orders in any case
            rel_scale_unc_item = np.array(
                fnlo.GetScaleUncertaintyVec(scale_var_type))
            rel_unc_list.append(rel_scale_unc_item)
            if verb:
                print('[fastnnlo_scaleunc]: \n')
                print('[fastnnlo_scaleunc]: Relative scale uncertainty in %s: \n' % n)
                # 3 entries: central value (xs), unc_low, unc_high
                print(rel_scale_unc_item, '\n')
                print(
                    '---------------------------------------------------------------------------------------')
                print(
                    '---------------------------------------------------------------------------------------')

        xs_all = np.array(xs_list)
        rel_scale_unc = np.array(rel_unc_list)
        ##########
        # structure of rel_scale_unc:
        # rel_scale_unc[0,:,:] means LO, rel_scale_unc[1,:,:] means NLO, and rel_scale_unc[2,:,:] in NNLO
        # rel_scale_unc[0,0,:] means xs in LO
        # rel_scale_unc[0,1,:] means rel. uncertainty upwards (in LO)
        # rel_scale_unc[0,2,:] means rel. uncertainty downwards (in LO)
        #########

        if verb:
            print('[fastnnlo_scaleunc]: Cross section xs_all uses ', order_list)

        if verb:
            print (xs_all, '\n \n')

        # ABSOLUTE scale uncertainty
        # length of axis 0 in xs_all equals number of orders
        num_orders = np.size(xs_all, 0)
        abs_scale_unc = np.empty([num_orders, 2, len(x_axis)])
        for k in range(0, len(xs_all)):
            abs_scale_unc[k, 0, :] = rel_scale_unc[k, 2, :] * \
                xs_all[k]  # absolute uncertainties downwards (low)
            abs_scale_unc[k, 1, :] = rel_scale_unc[k, 1, :] * \
                xs_all[k]  # absolute uncertainties upwards (high)

        if verb:
            print(
                '[fastnnlo_scaleunc]: Absolute Scale uncertainties downwards, upwards (order by order): \n')
            print(abs_scale_unc, '\n')

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
        if len(dxsr_cn) == 0:
            dxsr_cn = np.zeros((xs_all.shape[0], xs_all.shape[1]))
        if _debug:
            print('dxsr_cn', dxsr_cn)

        plotting(x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, dxsr_cn, nostat, xlabel, ylabel, title, tablename,
                 order_list, given_filename, scale_name, nice_scale_name, pdfset, variation_type, formats)

        stop_time = timeit.default_timer()
        timediff = stop_time-start_time
        print('fastnnlo_scaleunc: Elapsed time: %s sec = %s min' %
              (timediff, round(timediff/60, 2)))


if __name__ == '__main__':
    main()
