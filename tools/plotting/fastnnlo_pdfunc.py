#!/usr/bin/env python2
#-*- coding:utf-8 -*-

###########################################
#
# Plot the PDF uncertainty
#
#
# Created by B.Schillinger, 20.11.2018
# Modified by K. Rabbertz, 16.05.2019
#
###########################################
#
import argparse
import glob
import os
import re
import string
import sys
import timeit
# Use matplotlib with Cairo offline backend for eps, pdf, png, or svg output
import matplotlib as mpl
# mpl.use('Agg')
mpl.use('Cairo')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import (FormatStrFormatter, LogFormatter,
                               NullFormatter, ScalarFormatter, AutoMinorLocator, MultipleLocator)
from matplotlib import cm
# numpy
import numpy as np
# fastNLO for direct evaluation of interpolation grid
import fastnlo
from fastnlo import fastNLOLHAPDF
from fastnlo import SetGlobalVerbosity

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
_formats = {'eps': 0, 'pdf': 1, 'png': 2, 'svg': 3}
_text_to_order = {'LO': 0, 'NLO': 1, 'NNLO': 2}
_order_to_text = {0: 'LO', 1: 'NLO', 2: 'NNLO'}
_order_color = {'LO': 'g', 'NLO': 'b', 'NNLO': 'r'}
#_colors = ['tab:orange', 'tab:green', 'tab:purple', 'tab:blue', 'tab:brown']
#_colors         = ['darkorange', 'limegreen', 'mediumpurple', 'steelblue', 'saddlebrown']
_colors         = ['orange', 'green', 'purple', 'blue', 'brown']
#_symbols = ['s', 'X', 'o', '^', 'v']
_symbols = ['s', 'x', 'o', '^', 'v']
_hatches = ['-', '//', '\\', '|', '.']
_scale_to_text = {0: 'kScale1', 1: 'kScale2', 2: 'kQuadraticSum', 3: 'kQuadraticMean', 4: 'kQuadraticSumOver4',
                  5: 'kLinearMean', 6: 'kLinearSum', 7: 'kScaleMax', 8: 'kScaleMin', 9: 'kProd',
                  10: 'kS2plusS1half', 11: 'kPow4Sum', 12: 'kWgtAvg', 13: 'kS2plusS1fourth', 14: 'kExpProd2', 15: 'kExtern'}
_pdfbasenames = ['ABM11', 'ABMP15', 'ABMP16', 'CJ12', 'CJ15', 'CT10', 'CT14', 'HERAPDF20', 'JR14',
                 'MMHT2014', 'MSTW2008', 'NNPDF23', 'NNPDF30', 'NNPDF31']
_debug = False

#####################################################################################


# Function plotting cross sections for a list of PDF sets into one figure including PDF uncertainties
# The ratio is done always with respect to the first PDF appearing in the pdfsets list
# given e.g. via '-p CT14nnlo,MMHT2014nnlo68cl', i.e. the first evaluated cross section, here CT14nnlo,
# that is stored in xs_all[0].
# If multiple orders are requested, one plot per order is created.
# Otherwise all orders are plotted into one figure.
def plotting(x_axis, xmin, xmax, xs_all, rel_pdf_unc, abs_pdf_unc, dxsr_cn, nostat, xlabel, ylabel, title, tablename, order_list, given_filename, scale_name, nice_scale_name, pdfsets, formats, logx, logy):

    pdfnicenames = []
    for pdf in pdfsets:
        nicename = 'Undefined'
        for pdfn in _pdfbasenames:
            if pdfn in pdf:
                nicename = pdfn
        pdfnicenames.append(nicename)

    gs = gridspec.GridSpec(3, 3)

    # For plotting various results, max = 5, 'next to each other', handling via shift from bincenter
    if len(pdfsets) == 1:
        shift_list = [1.00]
    elif len(pdfsets) == 2:
        shift_list = [0.985, 1.015]
    elif len(pdfsets) == 3:
        shift_list = [0.97, 1.00, 1.03]
    elif len(pdfsets) == 4:
        shift_list = [0.955, 0.985, 1.015, 1.045]
    elif len(pdfsets) == 5:
        shift_list = [0.94, 0.97, 1.00, 1.03, 1.06]
    else:
        print('[fastnnlo_pdfunc]: Too many PDFs to plot simultaneously. Aborted!')
        print('[fastnnlo_pdfunc]: Current maximum is 5. The PDF list is', pdfsets)
        sys.exit()

    # Plot all x section results in order_list
    ord_index = -1
    ordernames = ''
    pdffilenames = ''
    for pdf in pdfsets:
        pdffilenames += '_%s' % pdf

    # Produce one plot per order with all PDFs
    for order_item in order_list:
        ord_index += 1

        fig = plt.figure(figsize=(7, 7))
        ax1 = plt.subplot(gs[:-1, :])

        print('[fastnnlo_pdfunc]: Producing %s plot.' % order_item)

        # Loop over PDF sets
        pdf_index = -1
        for pdf, shift in zip(pdfnicenames, shift_list):
            pdf_index += 1
            yerror = 0*xs_all[pdf_index, ord_index, :]
            if pdf_index == 0 and not nostat:
                yerror = np.multiply(
                    xs_all[pdf_index, ord_index, :], dxsr_cn[ord_index, :])
            elif nostat:
                yerror = abs(abs_pdf_unc[pdf_index, ord_index, :, :])
            ax1.errorbar(x_axis*shift, xs_all[pdf_index, ord_index, :], yerr=yerror,
                         elinewidth=1, linewidth=1.0, ms=6, marker=_symbols[pdf_index], color=_colors[pdf_index], fmt='.', label=pdf)
            ax1.fill_between(x_axis*shift, xs_all[pdf_index, ord_index, :] + xs_all[0, ord_index, :]*rel_pdf_unc[pdf_index, ord_index, 2, :],
                             xs_all[pdf_index, ord_index, :] + xs_all[0,
                                                                      ord_index, :]*rel_pdf_unc[pdf_index, ord_index, 1, :],
                             color=_colors[pdf_index], alpha=0.3, hatch=_hatches[pdf_index])

        # 'minor_thresholds' needs matplotlib > 1.5.0
        #        axfmt = LogFormatter(labelOnlyBase=False, minor_thresholds=(2, 0.9))
        ax1.set_xlim([xmin, xmax])
#                        if logx: ax1.set_xscale('log', nonposx='clip')
#                        else: ax1.set_xscale('linear')
#                        if logy: ax1.set_yscale('log', nonposy='clip')
#                        else: ax1.set_yscale('linear')
        ax1.set_xscale('log', nonposx='clip')
        #        ax1.get_xaxis().set_minor_formatter(axfmt)
#                        ax1.get_xaxis().set_minor_formatter(NullFormatter())
        ax1.set_yscale('log', nonposy='clip')
#                        ax1.set_xlabel(r'%s' %xlabel, horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
        ax1.set_ylabel(r'%s' % ylabel, horizontalalignment='right', x=1.0,
                       verticalalignment='top', y=1.0, rotation=90, labelpad=24)
        ax1.legend(fontsize=10, numpoints=1)
        ax1.text(0.03, 0.15, 'Reference PDF: %s' %
                 pdfnicenames[0], horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)
        ax1.text(0.03, 0.08, 'Scale: %s' % nice_scale_name, horizontalalignment='left',
                 verticalalignment='bottom', transform=ax1.transAxes)
        ax1.text(0.03, 0.03, 'PDF uncertainty at %s' % order_item, horizontalalignment='left',
                 verticalalignment='bottom', transform=ax1.transAxes)
        ax1.set_title('%s' % title, loc='left')

#        Only for publication
# H1
#        ax1.text(0.35, 0.90, r'$30 < Q^2 < 42\,\mathrm{GeV}^2$',
#                 horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)
# ZEUS
#        ax1.text(0.35, 0.90, r'$500 < Q^2 < 1000\,\mathrm{GeV}^2$',
#                 horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)

        # Ratio subplot with relative pdf uncertainties; denominator in ratio = first PDF in pdfsets list for requested order
        ax2 = plt.subplot(gs[2, :], sharex=ax1)
#        ax2.get_xaxis().set_minor_formatter(axfmt)
        ax2.set_yscale('linear', nonposy='clip')
        ax2.set_xlabel(r'%s' % xlabel, horizontalalignment='right',
                       x=1.0, verticalalignment='top', y=1.0)

        patches = []  # Later needed for legend
        for pdf_index in range(0, len(pdfsets)):
            # Divide xs for each PDF by first given PDF (xs)
            # Do not show error bars
            yerror = 0*xs_all[pdf_index, ord_index, :]
            mymarker = _symbols[pdf_index]
            mymksize = 4
            if nostat:
                # Error bars = PDF uncertainty
                #            yerror = abs(abs_pdf_unc[pdf_index, ord_index, :, :])
                #            mymksize = 4
                # No error bars, only small symbols
                yerror = 1*yerror
                #                mymarker = '.'
            else:
                # Error bars = statistical uncertainty if known
                # Only show for first PDF
                if pdf_index == 0:
                    yerror = np.multiply(
                        xs_all[pdf_index, ord_index, :], dxsr_cn[ord_index, :])
            ax2.errorbar(x_axis, xs_all[pdf_index, ord_index, :]/xs_all[0, ord_index, :], yerr=yerror/xs_all[0, ord_index, :], barsabove=True,
                         elinewidth=1, linewidth=1.0, ms=mymksize, marker=mymarker, color=_colors[pdf_index], fmt='.', label=pdf)
            ax2.fill_between(x_axis, (xs_all[pdf_index, ord_index, :]*(1+rel_pdf_unc[pdf_index, ord_index, 2, :])/xs_all[0, ord_index, :]),
                             (xs_all[pdf_index, ord_index, :]*(1+rel_pdf_unc[pdf_index,
                                                                             ord_index, 1, :])/xs_all[0, ord_index, :]),
                             color=_colors[pdf_index], alpha=0.3, hatch=_hatches[pdf_index])
            patches.append(mpl.patches.Rectangle(
                (0, 0), 0, 0, color=_colors[pdf_index], label=pdfnicenames[pdf_index], alpha=0.4))
            ax2.add_patch(patches[pdf_index])

        ax2.set_ylabel(r'Ratio to ref. PDF', horizontalalignment='center',
                       x=1.0, verticalalignment='top', y=0.5, rotation=90, labelpad=24)
        ax2.axhline(y=1, xmin=0, xmax=1, color='k',
                    linestyle='dotted', linewidth=1.6, alpha=0.2)

        fig.tight_layout()
        ordernames += '_%s' % order_item

        if given_filename is not None:
            filename = '%s.pdfunc-%s.%s' % (given_filename,
                                            order_item, pdffilenames[1:])
        else:
            filename = '%s.pdfunc-%s.%s.%s' % (tablename,
                                               order_item, pdffilenames[1:], scale_name)
            if not nostat:
                filename = filename+'.stat'
        for fmt in formats:
            figname = '%s.%s' % (filename, fmt)
            fig.savefig(figname)
            print('[fastnnlo_pdfunc]: Plot saved as:', figname)
        plt.close(fig)


#####################################################################################

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
                        help='Comma-separated list of NNLOJET dat files with statistical uncertainties for each order to show. If set to "auto", dat files matching to table name are used. If nothing is chosen, statistical uncertainties are ignored.')
    parser.add_argument('-f', '--filename', default=None, type=str,
                        help='Output filename (optional).')
    parser.add_argument('--format', required=False, nargs='?', type=str, action=SplitArgs,
                        help='Comma-separated list of plot formats to use: eps, pdf, png, svg. If nothing is chosen, png is used.')
    parser.add_argument('--logx', default=True, required=False, nargs=1, type=bool,
                        help='Switch between linear and logarithmic x axis. NOT working yet!')
    parser.add_argument('--logy', default=True, required=False, nargs=1, type=bool,
                        help='Switch between linear and logarithmic y axis. NOT working yet!')
    parser.add_argument('-m', '--member', default=0, type=int,
                        help='Member of PDFset, default is 0.')
    parser.add_argument('-o', '--order', required=False, nargs='?', type=str, action=SplitArgs,
                        help='Comma-separated list of orders to show: LO, NLO, and/or NNLO. If nothing is chosen, show all orders available in table.')
    parser.add_argument('-p', '--pdfset', required=False, nargs='?', type=str, action=SplitArgs,
                        default=['CT14nnlo'],
                        help='Comma-separated list of PDF sets to use.')
    parser.add_argument('-s', '--scale', default=0, required=False, nargs='?', type=int,
                        choices=range(16), metavar='[0-15]',
                        help='For flexible-scale tables define central scale choice for MuR and MuF by selection enum fastNLO::ScaleFunctionalForm ("0"=kScale1, "1"=kScale2, "2"=kQuadraticSum), ...')
    parser.add_argument('--scalename', default=None, type=str,
                        help='Replace default scale name by given string.')
    parser.add_argument('--title', default=None, type=str,
                        help='Replace table name as default title by given string.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Increase output verbosity.')
    parser.add_argument('--xlabel', default=None, type=str,
                        help='Replace x axis default label by given string.')
    parser.add_argument('--ylabel', default=None, type=str,
                        help='Replace y axis default label by given string.')

    # Parse arguments
    args = vars(parser.parse_args())

    # List of table names
    files = args['table']
    print('\n')
    print('[fastnnlo_pdfunc]: Analysing table list: ')
    for file in files:
        print('[fastnnlo_pdfunc]:   ', file)

    # PDF set name
    pdfsets = []
    if args['pdfset'] is None:
        pdfsets = ['CT14nnlo']
        print('[fastnnlo_pdfunc]: Using default PDF set: ', pdfsets[0])
    else:
        for pdf in args['pdfset']:
            pdfsets.append(pdf)
        print('[fastnnlo_pdfunc]: Using PDF set(s): ', pdfsets)

    # Orders to be shown
    iorders = []
    iordmin = _text_to_order['LO']
    iordmax = _text_to_order['NNLO']
    if args['order'] is None:
        print('[fastnnlo_pdfunc]: Evaluate table up to highest available order.')
    else:
        for ord in args['order']:
            if _text_to_order.has_key(ord):
                iorders.append(_text_to_order[ord])
            else:
                print('[fastnnlo_pdfunc]: Illegal order specified, aborted!')
                print('[fastnnlo_pdfunc]: Order list:', args['order'])
                exit(1)
        iordmin = min(iorders)
        iordmax = max(iorders)
        print('[fastnnlo_pdfunc]: Evaluate table up to order(s)', args['order'])

    # Check existence of NNLOJET dat file for each order if desired
    datfilenames = []
    nostat = False
    if args['datfiles'] is None:
        nostat = True
        print('[fastnnlo_pdfunc]: No statistical uncertainties requested.')
    elif args['datfiles'][0] == 'auto':
        print('[fastnnlo_pdfunc]: Automatic filename matching is used to load statistical uncertainties from NNLOJET.')
    else:
        for datfile in args['datfiles']:
            lstat = os.path.isfile(datfile)
            if lstat:
                datfilenames.append(datfile)
            else:
                print('[fastnnlo_pdfunc]: Given file ', datfile, 'for statistical uncertainties not found, aborted!')
                exit(1)

    # Given filename
    given_filename = args['filename']

    # Plot formats to use
    formats = args['format']
    if formats is None:
        formats = ['png']
    for fmt in formats:
        if not _formats.has_key(fmt):
            print('[fastnnlo_pdfunc]: Illegal format specified, aborted!')
            print('[fastnnlo_pdfunc]: Format list:', args['format'])
            exit(1)

    # Logarithmic or linear x axis
    logx = args['logx']

    # Logarithmic or linear x axis
    logy = args['logy']

    # Scale choice
    scale_choice = args['scale']

    # Scale name
    nice_scale_name = args['scalename']

    # Plot labelling
    nice_title = args['title']
    nice_xlabel = args['xlabel']
    nice_ylabel = args['ylabel']

    # Verbosity
    verb = args['verbose']

    # Loop over table list
    for table in files:
        print('[fastnnlo_pdfunc]: Analysing table: ', table)
        # Get rid of extensions (.tab.gz or .tab)
        tablename = os.path.splitext(os.path.basename(table))[0]
        tablename = os.path.splitext(tablename)[0]

        ###################### Start EVALUATION with fastNLO library ###################################################
        # SetGlobalVerbosity(0) # Does not work since changed to default in the following call
        # Take the general information (bin_bounds, labels, order_existence, etc.) from first given pdfset.
        fnlo = fastnlo.fastNLOLHAPDF(table, pdfsets[0], args['member'])

        # Get labeling for the x-axis
        # Dimensionality of the table:
        ndim = fnlo.GetNumDiffBin()
        if verb:
            print('\n')
            print('[fastnnlo_pdfunc]: Table Dimensions: ', ndim)

        # Labels of all the dimensions:
        dimlabels = fnlo.GetDimLabels()
        if verb:
            print('[fastnnlo_pdfunc]: Dimension Labels: ', dimlabels)

        # x label of first dimension from table:
        xlabel = fnlo.GetDimLabel(0)
        if verb:
            print('[fastnnlo_pdfunc]: x-label: ', xlabel)

        # Generic y label
        ylabel = '$\sigma \pm \Delta\sigma(\mathrm{PDF})$'

        # Creating x-axis
        bin_bounds = np.array(fnlo.GetObsBinsBounds(0))
        if verb:
            print('[fastnnlo_pdfunc]: bin_bounds.T: \n', bin_bounds.T, '\n')
            print('[fastnnlo_pdfunc]: bin_bounds.flatten(): \n', bin_bounds.flatten(), '\n')

        x_axis = (bin_bounds.T[0]+bin_bounds.T[1]) / \
            2.  # this is a list of bin centers
        xmin = 0.95*min(bin_bounds.ravel())
        xmax = 1.05*max(bin_bounds.ravel())
        if verb:
            print('[fastnnlo_pdfunc]: xmin=%s, xmax=%s. \n' % (xmin, xmax))

        # Preparing x-errors (via bin_bounds) --> x_errors[0, :] are initially negative (positive via -1*), x_errors[1, :] positive
        x_errors = np.array(
            [-1*(bin_bounds.T[0]-x_axis), bin_bounds.T[1]-x_axis])
        if verb:
            print('[fastnnlo_pdfunc]: x_errors: \n', x_errors, '\n')

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
                        print('[fastnnlo_pdfunc]: Invalid choice of scale = ', scale_choice, ' Aborted!')
                        print('[fastnnlo_pdfunc]: For fixed-scale tables only the default=0 is allowed.')
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
            print('[fastnnlo_pdfunc]: Table has continuous orders up to', cnt_order, 'and a maximal order of', max_order)

        if iordmax > cnt_order:
            print('[fastnnlo_pdfunc]: Invalid choice of orders. Aborted!')
            print('[fastnnlo_pdfunc]: Highest order requested is', _order_to_text[iordmax], 'but orders are available only up to', cnt_order)
            exit(1)

        order_list = []
        if args['order'] is None:
            for iord in range(cnt_order+1):
                order_list.append(_order_to_text[iord])
        else:
            order_list = args['order']

        print('[fastnnlo_pdfunc]: List of requested orders:', order_list)

        # Read in statistical uncertainty for each order if requested
        dxsr_cn = []
        sep = '.'
        if not nostat:
            if args['datfiles'][0] == 'auto':
                for order in order_list:
                    parts = tablename.split(sep)
                    parts[1] = order
                    datfile = sep.join(parts) + '.dat'
                    datfilenames.append(datfile)

            lstat = (len(datfilenames) > 0)
            if lstat and len(datfilenames) != len(order_list):
                print('[fastnnlo_pdfunc]: Mismatch between no. of requested orders and no. of filenames for statistical uncertainties, aborted!')
                exit(1)

            dxsr = []
            for fname in datfilenames:
                print('[fastnnlo_pdfunc]: Taking statistical uncertainties from', fname)
                cols = np.loadtxt(fname, usecols=range(3, 5))
                xs_dat = np.array(cols[:, 0])
                dxs_dat = np.array(cols[:, 1])
                dxsr_dat = np.divide(dxs_dat, xs_dat, out=np.ones_like(
                    dxs_dat), where=xs_dat != 0)
                dxsr.append(dxsr_dat)
            dxsr_cn = abs(np.array(dxsr))
            # Empty list for use with next table in automatic mode
            if args['datfiles'][0] == 'auto':
                datfilenames = []

        # For flexible-scale tables set scale to user choice (default is 0)
        if lflex:
            print('[fastnnlo_pdfunc]: Setting requested scale choice for flexible-scale table:', scale_choice)
            fnlo.SetMuRFunctionalForm(scale_choice)
            fnlo.SetMuFFunctionalForm(scale_choice)
        else:
            if scale_choice == 0:
                print('[fastnnlo_pdfunc]: Evaluating fixed-scale table. Scale choice must be', scale_choice)
            else:
                print('[fastnnlo_pdfunc]: No scale choice possible for fixed-scale table. Aborted!')
                print('[fastnnlo_pdfunc]: scale_choice = ', scale_choice)
                exit(1)

        # Now evaluate fastNLO table for each of the given PDF sets
        # will contain for each given pdf the total cross section for LO, NLO, NNLO (or requested orders)
        xs_list = []
        # list for relative pdf uncertainties (low, high) for each PDF in LO, NLO, NNLO
        rel_unc_list = []

        # Outermost loop: go through all the given pdf sets
        # Evaluate fastNLO table focusing on PDF uncertainties
        for pdf in pdfsets:
            xs_list_tmp = []  # will contain total cross section for selected orders out of LO, NLO, NNLO for single pdf
            # list for relative PDF uncertainties (low, high) for selected orders
            rel_unc_list_tmp = []
            print('-----------------------------------------------------------------------------------------------')
            print('#############################  %s  ########################################' % pdf)
            print('-----------------------------------------------------------------------------------------------')
            print('[fastnnlo_pdfunc]: Calculate XS and uncertainty for %s \n' % pdf)
            fnlo.SetLHAPDFFilename(pdf)
            fnlo.SetLHAPDFMember(0)

            for n in order_list:
                for j in range(0, max_order+1):
                    if j <= _text_to_order[n]:
                        fnlo.SetContributionON(fastnlo.kFixedOrder, j, True)
                    else:
                        fnlo.SetContributionON(fastnlo.kFixedOrder, j, False)
                if verb:
                    print('[fastnnlo_pdfunc]: Calculate XS for order: %s' % n)
                    print('[fastnnlo_pdfunc]: ---- ---- ---- ---- ---- ---- ---- \n')

                fnlo.CalcCrossSection()
                xs_list_tmp.append(fnlo.GetCrossSection())

                ### Get PDF uncertainties ###
                rel_pdf_unc_item = np.array(fnlo.GetPDFUncertaintyVec(
                    fastnlo.kLHAPDF6))  # Now calculated for order n
                rel_unc_list_tmp.append(rel_pdf_unc_item)
                if verb:
                    print('[fastnnlo_pdfunc]: \n')
                    print('[fastnnlo_pdfunc]: Relative PDF uncertainty in %s: \n' % n)
                    # has 3 entries: central value (xs), unc_high, unc_low
                    print(rel_pdf_unc_item, '\n')
                    print('-----------------------------------------------------------------------------------------------')
                    print('-----------------------------------------------------------------------------------------------')
            xs_list.append(xs_list_tmp)
            rel_unc_list.append(rel_unc_list_tmp)
        xs_all = np.array(xs_list)
        # Remember: both arrays here do only contain the CHOSEN orders! (or per default LO, NLO, (NNLO))
        rel_pdf_unc = np.array(rel_unc_list)
        #########                                                                                               #########
        # structure of rel_pdf_unc:                                                                                     #
        # rel_pdf_unc[0,:,:] refers to the order stored in order_list[0] --> this (0th item) is not necessarily LO!     #
        # The same applies to rel_pdf_unc[1,:,:] and so on.                                                             #
        # rel_pdf_unc[0,0,:] corresponds to the cross section in the order stored in order_list[0].                     #
        # rel_pdf_unc[0,1,:] corresponds to the rel. uncertainty upwards (in the order order_list[0])                   #
        # rel_pdf_unc[0,2,:] corresponds to the rel. uncertainty downwards (in the order order_list[0])                 #
        ########                                                                                                #########

        if verb:
            print('[fastnnlo_pdfunc]: Cross section summary. ')
            print('[fastnnlo_pdfunc]: For each pdf %s xs_all contains XS corresponding to %s.' % (pdfsets, order_list))
            print(xs_all, '\n \n')
            print('Size xs_all: ', np.shape(xs_all))  # test
            print('Size rel_pdf_unc: ', np.shape(rel_pdf_unc))  # test

        # ABSOLUTE pdf uncertainty
        # length of axis 0 in xs_all equals number of (investigated) pdfsets
        num_pdfsets = np.size(xs_all, 0)
        # length of axis 1 in xs_all equals number of (investigated) orders
        num_orders = np.size(xs_all, 1)
        abs_pdf_unc = np.empty([num_pdfsets, num_orders, 2, len(x_axis)])
        for p in range(0, num_pdfsets):
            for k in range(0, num_orders):  # k is order index
                # absolute uncertainties downwards (low)
                abs_pdf_unc[p, k, 0, :] = rel_pdf_unc[p,
                                                      k, 2, :]*xs_all[p, k, :]
                # absolute uncertainties upwards (high)
                abs_pdf_unc[p, k, 1, :] = rel_pdf_unc[p,
                                                      k, 1, :]*xs_all[p, k, :]

        if verb:
            print('[fastnnlo_pdfunc]: Absolute PDF uncertainties downwards, upwards for %s in %s: \n' % (pdfsets, order_list))
            print(abs_pdf_unc, '\n')

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
            dxsr_cn = np.zeros((xs_all.shape[0], xs_all.shape[2]))
        if _debug:
            print('dxsr_cn', dxsr_cn)

        plotting(x_axis, xmin, xmax, xs_all, rel_pdf_unc, abs_pdf_unc, dxsr_cn, nostat, xlabel, ylabel, title, tablename,
                 order_list, given_filename, scale_name, nice_scale_name, pdfsets, formats, logx, logy)

    stop_time = timeit.default_timer()
    timediff = stop_time-start_time
    print('fastnnlo_pdfunc: Elapsed time: %s sec = %s min' % (timediff, round(timediff/60., 2)))


if __name__ == '__main__':
    main()
