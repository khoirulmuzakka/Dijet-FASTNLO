#!/usr/bin/env python2
#-*- coding:utf-8 -*-

#########################################
#
# Plotting of the pdf uncertainty.
# Updated for tables with
# multidimensional binning.
#
#
# Created by B.Schillinger, 20.11.2018
# Last modified: 18.01.2019
# Updated to 3dim tables: 08.03.2019
#
#
#
########################################

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.pyplot import cm
import matplotlib.ticker
import os
import sys
import timeit
#import matplotlib.style
# matplotlib.style.use('classic')


import fastnlo
from fastnlo import SetGlobalVerbosity

_text_to_order = {'LO': 0, 'NLO': 1, 'NNLO': 2}
_order_to_text = {0: 'LO', 1: 'NLO', 2: 'NNLO'}
_order_color = {'LO': 'g', 'NLO': 'b', 'NNLO': 'r'}
_order_hatch = {'LO': '-', 'NLO': '//', 'NNLO': '\\'}  # just for testing

#######################################################################################################
# function for handling tables with binning in mutliple dimensions (e.g. y, p_T,...)


def dimensionhandler(fnlo, verb):

    # Get labeling for the x-axis
    # Dimensionality of the table:
    ndim = fnlo.GetNumDiffBin()
    if verb:
        print '\n'
        print '[fastnnlo_pdfunc]: Table Dimensions: ', ndim

    # Labels of all the dimensions:
    dimlabels = fnlo.GetDimLabels()
    if verb:
        print '[fastnnlo_pdfunc]: Dimension Labels: ', dimlabels

    # choose binning dimension (innermost, fnlo name: 'IODimI'. For number of observable bin (starting from 1) see for instance fnlo-tk-cppread output and code)
    # number of dimensions minus 1 (starting indexing at 0)
    binningdim = ndim-1

    # Label of binning dimension:
    xlabel = fnlo.GetDimLabel(binningdim)
    if verb:
        print '[fastnnlo_pdfunc]: x-label: ', xlabel

    # Depending on table dimension, a different amount of plots will be needed
    # --> check if lower bin bound is lower than previous lower bin bound (indicates new phase space slice in table)
    # --> create several x-axes accordingly
    bbDimO = []  # won't be used again if table is 1dim, otherwise it is set anew below
    bbDimM = []  # ------------------------- " -------------------------------------
    if (ndim == 2):
        bbDimO = np.array(fnlo.GetObsBinsBounds(0))
    elif (ndim == 3):
        bbDimO = np.array(fnlo.GetObsBinsBounds(0))
        bbDimM = np.array(fnlo.GetObsBinsBounds(1))

    # list in which bins show up repeatedly if ndim>1
    total_binbounds = np.array(fnlo.GetObsBinsBounds(binningdim))
    #print "TOTAL BINBOUNDS", total_binbounds
    #print "LEN(TOTALBINBOUNDS)", len(total_binbounds[:, 0])
    tmp_binbounds = total_binbounds
    tmp_len = 0  # keeps track of the index of the last lower bin bound
    nplot = 0  # counts how many plots are needed for given table
    plots_binbounds_list = []
    x_axis_list = []
    # contains the index of the "last lower bin bound" for each plot to be produced
    tmp_len_list = []
    while (len(tmp_binbounds) > 0):
        # go through all lower binbounds
        for i in range(1, len(tmp_binbounds[:, 0])):
            # check: when does "next x-axis" start? // when do lower bin bounds get smaller again
            if (tmp_binbounds[i, 0] <= tmp_binbounds[(i-1), 0]):
                # only first time adding i-1 (because first index in total_binbounds etc. is 0)
                #...otherwise add i
                if nplot == 0:
                    tmp_len += (i-1)
                else:
                    # x-axis of plot goes up to index(!) tmp_len (later, this is necessary for plotting routine)
                    tmp_len += i
                plots_binbounds_list.append(tmp_binbounds[0:i, :])
                tmp_len_list.append(tmp_len)

                # update temporary binbounds-array
                newlo = tmp_len+1
                tmp_binbounds = total_binbounds[newlo:, :]

                # increase number of plots to be produced by 1
                nplot += 1
                break

            # How to end this loop? --> identify and cut out last tmp_binbounds
            # happens if the current tmp_binbounds is the last array (= lower bin bound does not decrease anymore.)
            elif (i == (len(tmp_binbounds[:, 0])-1)):
                # append last plot's bin bounds:
                plots_binbounds_list.append(tmp_binbounds)
                # because of elif-condition here, index had no chance to rise one step further
                tmp_len += (i+1)
                tmp_len_list.append(tmp_len)
                nplot += 1
                tmp_binbounds = []
                tmp_binbounds = np.array(tmp_binbounds)

    # this is a numpy array containing lists of (usually) different length.
    all_bin_bounds = np.array(plots_binbounds_list)
    # contains one list with bin bounds (lower, upper) for each plot.

    if verb:
        # , all_bin_bounds.T, '\n'
        print '[fastnnlo_pdfunc]: Bin bounds for all plots: \n'
        #print '[fastnnlo_pdfunc]: bin_bounds.flatten()', all_bin_bounds.flatten(), '\n'

    # Creating x-axis and 'x-errors'
    xmin_list = []
    xmax_list = []
    x_errors = []  # stays list because the contained arrays are of different length
    for n in range(0, nplot):
        if verb:
            print '[fastnnlo_pdfunc]: ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----'
            print 'shape of all_bin_bounds: ', np.shape(all_bin_bounds)
            print '%s. plot (%s bins) bin bounds (lower, upper): \n' % (
                n, len(all_bin_bounds[n]))
            # print nth list in all_bin_bounds array
            print 'all_bin_bounds[n].T \n', all_bin_bounds[n].T, '\n'
        # this is a list which contains one list of bincenters for each plot
        x_axis_list.append((all_bin_bounds[n].T[0]+all_bin_bounds[n].T[1])/2.)
        xmin_list.append(0.95*min(all_bin_bounds[n].ravel()))
        xmax_list.append(1.05*max(all_bin_bounds[n].ravel()))
        if verb:
            print '[fastnnlo_pdfunc]: xmin=%s, xmax=%s. \n' % (
                xmin_list[n], xmax_list[n])

        # Preparing x-errors (via bin bounds) --> x_errors[n, 0, :] are initially negative (positive via -1*), x_errors[n, 1, :] positive
        x_errors_tmp = np.array(
            [-1*(all_bin_bounds[n].T[0]-x_axis_list[n]), all_bin_bounds[n].T[1]-x_axis_list[n]])
        x_errors.append(x_errors_tmp)

        if verb:
            print 'fastnnlo_pdfunc]: x_errors for plot %s: \n' % n, x_errors[n], '\n'

    return ndim, binningdim, xlabel, all_bin_bounds, x_axis_list, xmin_list, xmax_list, x_errors, nplot, dimlabels, tmp_len_list, bbDimO, bbDimM


#################################################################################################################


def plotting(pdfsets, order_list, xs_chosen, abs_pdf_unc, rel_pdf_unc, xlabel, tablename, given_filename, x_axis, xmin, xmax, nameDimO, nameDimM, nameDimI, dimlabels, plotnr):

    gs = gridspec.GridSpec(3, 3)
    fig = plt.figure(figsize=(7, 7))
    ax1 = plt.subplot(gs[:-1, :])

    # specification of dimensions for file name:
    if len(dimlabels) == 3:
        dimspeci = '%s.%s.%s' % (dimlabels[0], dimlabels[1], dimlabels[2])
    elif len(dimlabels) == 2:
        dimspeci = '%s.%s' % (dimlabels[0], dimlabels[1])
    else:
        dimspeci = '%s' % dimlabels[0]  # ==xlabel

    # Check whether only one PDF set is being looked at or several PDF sets.
    if (len(pdfsets) > 1):  # one plot contains all chosen pdf sets for one specific order + summary plot if 3 orders in order_list
        # Several PDFs will be plotted 'next to each other' --> need tiny shift to be displayed
        if len(pdfsets) == 2:
            shift_list = [0.99, 1.01]
        elif len(pdfsets) == 3:
            shift_list = [0.98, 1.00, 1.02]
        elif len(pdfsets) == 4:
            shift_list = [0.97, 0.99, 1.01, 1.03]
        elif len(pdfsets) == 5:
            shift_list = [0.96, 0.98, 1.00, 1.02, 1.04]
        else:
            print '[fastnnlo_pdfunc]: Too many PDFs to produce combined plot. Aborted!'
            print '[fastnnlo_pdfunc]: Max. possible: 5. Number of given PDF Sets: %s' % len(
                pdfsets)
            sys.exit()

        order_index = 0  # necessary, because even though LO might not even be investigated, we need index 0, which could then correspond to NLO or NNLO
        ordernames = ''
        pdfnames = ''
        for pdf in pdfsets:
            pdfnames += '_%s' % pdf

        # producing one plot per order (comparing all the given pdfs)
        for order_item in order_list:
            fig = plt.figure(figsize=(7, 7))
            ax1 = plt.subplot(gs[:-1, :])

            print '[fastnnlo_pdfunc]: Producing %s plot.' % order_item
            color = iter(cm.rainbow(np.linspace(0, 1, len(pdfsets))))
            patches = []  # later needed for legend

            pdf_index = 0
            for pdf, shift in zip(pdfsets, shift_list):  # go through all the pdfs
                c = next(color)
                ax1.errorbar(x_axis*shift, xs_chosen[pdf_index, order_index, :], yerr=abs(abs_pdf_unc[pdf_index, order_index, :, :]), elinewidth=1,
                             linewidth=0.0, ms=4, color=c, fmt='.', label=pdf)
                pdf_index += 1

            ax1.set_xscale('log', nonposx='clip')
            ax1.set_yscale('log', nonposy='clip')
            ax1.set_xlabel('%s' % xlabel)
            ax1.set_ylabel('XS with abs. pdf_unc', rotation=90)
            ax1.legend(fontsize=10, numpoints=1)
            ax1.text(0.06, 0.06, '%s' % order_item, fontsize=12, horizontalalignment='left',
                     verticalalignment='bottom', transform=ax1.transAxes)

            ax1.text(0.06, 0.12, '%s' % nameDimO, fontsize=12, horizontalalignment='left',
                     verticalalignment='bottom', transform=ax1.transAxes)
            ax1.text(0.06, 0.18, '%s' % nameDimM, fontsize=12, horizontalalignment='left',
                     verticalalignment='bottom', transform=ax1.transAxes)

            ax1.set_title('%s' % tablename)

            ax2 = plt.subplot(gs[2, :])
            ax2.set_xlim([xmin, xmax])

            # ratioplot normalised to first given pdf.
            color = iter(cm.rainbow(np.linspace(0, 1, len(pdfsets))))
            for p in range(0, len(pdfsets)):  # for pdf in pdfsets
                c = next(color)
                # divide xs for each PDF by first given PDF (xs)
                ax2.semilogx(x_axis, xs_chosen[p, order_index, :]/xs_chosen[0, order_index, :], '.', ms=4, ls='dashed', linewidth=0.2,
                             color=c, label=pdfsets[p])

                ax2.fill_between(x_axis, (xs_chosen[p, order_index, :]/xs_chosen[0, order_index, :])+rel_pdf_unc[p, order_index, 2, :],
                                 (xs_chosen[p, order_index, :]/xs_chosen[0, order_index, :])+rel_pdf_unc[p, order_index, 1, :], color=c,
                                 alpha=0.3)
                # patch for current pdf (needed for labeling/legend)
                patches.append(matplotlib.patches.Rectangle(
                    (0, 0), 0, 0, color=c, label=pdfsets[p], alpha=0.4))
                ax2.add_patch(patches[p])

            ax2.set_ylabel('rel. pdf_unc')
            ax2.text(0.02, 0.92, 'Reference PDF: %s' %
                     pdfsets[0], horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes)
            ax2.axhline(y=1, xmin=0, xmax=1, linewidth=0.4,
                        color='k', linestyle='dotted')

            fig.tight_layout()

            order_index += 1
            ordernames += '_%s' % order_item

            if given_filename is not None:
                filename = '%s.pdfunc-%s.plot%s.%s.%s' % (
                    given_filename, order_item, plotnr, dimspeci, pdfnames[1:])
            else:
                filename = '%s.pdfunc-%s.plot%s.%s%s' % (
                    tablename, order_item, plotnr, dimspeci, pdfnames[1:])

            fig.savefig('%s.png' % filename)
            print '[fastnnlo_pdfunc]: %s plot saved as: %s.png' % (
                order_item, filename)
            plt.close()

        plt.close()
        # Summary plot if 3 orders have been looked at
        if len(order_list) == 3:  # could adjust this to arbitrary length of order
            gs2 = gridspec.GridSpec(2, 2)
            # for more whitespace between the subplots
            gs2.update(wspace=0.4, hspace=0.4)
            ##fig2 = plt.figure(figsize=(7,7))
            # different size compared to other plots that are produced
            fig2 = plt.figure(figsize=(12, 10))
            ax0 = plt.subplot(gs2[0, 0])
            ax1 = plt.subplot(gs2[0, 1])
            ax2 = plt.subplot(gs2[1, 0])
            ax3 = plt.subplot(gs2[1, 1])

            # subplots containing ratioplot for each order (comparing all pdfs, ratio to first given pdf)
            # and one subplot for first pdf (comparing all orders, ratio to 1st given order)
            # for ax, order in zip([ax0, ax1, ax2], ['LO', 'NLO', 'NNLO']): #here ax0 will always contain LO, ax1 NLO and ax2 NNLO
            order_index = 0  # to go through order_list item by item, keeping in mind that order_list is not sorted
            patches = []  # later needed for legend
            # here it can happen, that NNLO is plotted in ax0, depending on content of order_list
            for ax, order in zip([ax0, ax1, ax2], order_list):
                patches = []  # currently necessary, not needed if one only makes one legend for all axes
                color = iter(cm.rainbow(np.linspace(0, 1, len(pdfsets))))
                for p in range(0, len(pdfsets)):  # for pdf in pdfsets
                    c = next(color)
                    ax.semilogx(x_axis, xs_chosen[p, order_index, :]/xs_chosen[0, order_index, :], '.', ms=4, ls='dashed', linewidth=0.2,
                                color=c)  # , label=pdfsets[p])
                    ax.fill_between(x_axis, (xs_chosen[p, order_index, :]/xs_chosen[0, order_index, :])+rel_pdf_unc[p, order_index, 2, :],
                                    (xs_chosen[p, order_index, :]/xs_chosen[0, order_index, :])+rel_pdf_unc[p, order_index, 1, :], color=c,
                                    alpha=0.3)
                    # Alternative: create the patches only once (for ax0) --> use the same for ax1 and ax2
                    # if ax==ax0:
                    ###     patches.append(matplotlib.patches.Rectangle((0, 0), 0, 0, color=c, label=pdfsets[p], alpha=0.3))
                    # ax.add_patch(patches[p])
                    patches.append(matplotlib.patches.Rectangle(
                        (0, 0), 0, 0, color=c, label=pdfsets[p], alpha=0.3))
                    ax.add_patch(patches[p])

                ax.set_xlabel('%s' % xlabel)
                ax.set_ylabel('rel. pdf_unc')
                ax.text(0.02, 0.12, 'Reference PDF: %s' %
                        pdfsets[0], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
                ax.text(0.02, 0.06, 'Order: %s' %
                        _order_to_text[order_index], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)

                ax.text(0.06, 0.18, '%s' % nameDimO, fontsize=12, horizontalalignment='left',
                        verticalalignment='bottom', transform=ax1.transAxes)
                ax.text(0.06, 0.24, '%s' % nameDimM, fontsize=12, horizontalalignment='left',
                        verticalalignment='bottom', transform=ax1.transAxes)

                ax.axhline(y=1, xmin=0, xmax=1, linewidth=0.4,
                           color='k', linestyle='dotted')

                order_index += 1
                ax.legend()
            #plt.legend(handles=patches, mode='expand', bbox_to_anchor=(1.02, 1), loc='upper left') #'upper left' refers to the legend box ######

            order_index = 0
            for order_item in order_list:
                ax3.semilogx(x_axis, xs_chosen[0, order_index, :]/xs_chosen[0, 0, :], '.', ms=4, ls='dashed', linewidth=0.2,
                             color=_order_color[order_item], label=order_item)

                ax3.fill_between(x_axis, (xs_chosen[0, order_index, :]/xs_chosen[0, 0, :])+rel_pdf_unc[0, order_index, 2, :],
                                 (xs_chosen[0, order_index, :]/xs_chosen[0, 0, :])+rel_pdf_unc[0, order_index, 1, :], color=_order_color[order_item],
                                 alpha=0.46)  # , hatch=_order_hatch[order_item])

                order_index += 1
            ax3.set_title('%s' % pdfsets[0])
            ax3.set_xlabel('%s' % xlabel)
            ax3.set_ylabel('rel. pdf_unc')
            ax3.set_ylim([0.82, None])  # to avoid problems with text
            ax3.axhline(y=1, xmin=0, xmax=1, linewidth=0.4,
                        color='k', linestyle='dotted')
            ax3.text(0.02, 0.08, 'Reference order: %s' %
                     order_list[0], horizontalalignment='left', verticalalignment='bottom', transform=ax3.transAxes)
            ax3.text(0.02, 0.02, 'PDF Set: %s' %
                     pdfsets[0], horizontalalignment='left', verticalalignment='bottom', transform=ax3.transAxes)
            ax3.legend()

            # adjust fontsize to smaller size? (looks okay so far)
            ax3.text(0.02, 0.14, '%s' % nameDimO, fontsize=12, horizontalalignment='left',
                     verticalalignment='bottom', transform=ax1.transAxes)
            ax3.text(0.02, 0.20, '%s' % nameDimM, fontsize=12, horizontalalignment='left',
                     verticalalignment='bottom', transform=ax1.transAxes)

            fig2.suptitle('%s' % tablename, fontsize=16)

            if given_filename is not None:
                filename = '%s.pdfunc-summary.plot%s.%s%s' % (
                    given_filename, plotnr, dimspeci, pdfnames[1:])
            else:
                filename = '%s.pdfunc-summary.plot%s.%s%s' % (
                    tablename, plotnr, dimspeci, pdfnames[1:])

            fig2.savefig('%s.png' % filename, bbox_inches='tight')
            print '[fastnnlo_pdfunc]: Summary plot saved as: %s.png' % filename

        else:
            print '[fastnnlo_pdfunc]: For less than 3 orders --> no summary plot produced.'

    elif (len(pdfsets) == 1):  # plotting with given or default pdf, one plot contains all given/default orders
        # For plotting several orders 'next to each other', handling via shift from bincenter
        if len(order_list) == 1:
            shift_list = [1.00]
        elif len(order_list) == 2:
            shift_list = [0.98, 1.02]
        elif len(order_list) == 3:
            shift_list = [0.98, 1.00, 1.02]
        else:
            print '[fastnnlo_pdfunc]: Too many orders to plot. Aborted!'
            print '[fastnnlo_pdfunc]: Requested orders: %s' % order_list
            sys.exit()

        # this is necessary because the orders in order_list are unsorted (order_index is not tied to specific order)
        order_index = 0
        for order_item, shift in zip(order_list, shift_list):
            ax1.errorbar(x_axis*shift, xs_chosen[0, order_index, :], yerr=abs(abs_pdf_unc[0, order_index, :, :]), elinewidth=1,
                         linewidth=0.0, ms=4, color=_order_color[order_item], fmt='.', label=order_item)
            order_index += 1

        ax1.set_xscale('log', nonposx='clip')
        ax1.set_yscale('log', nonposy='clip')
        ax1.set_xlabel('%s' % xlabel)
        ax1.set_ylabel('XS with abs. pdf_unc', rotation=90)
        ax1.legend(fontsize=10, numpoints=1)
        ax1.text(0.02, 0.06, '%s' % pdfsets[0], horizontalalignment='left',
                 verticalalignment='bottom', transform=ax1.transAxes)

        ax1.text(0.02, 0.12, '%s' % nameDimO, fontsize=12, horizontalalignment='left',
                 verticalalignment='bottom', transform=ax1.transAxes)
        ax1.text(0.02, 0.18, '%s' % nameDimM, fontsize=12, horizontalalignment='left',
                 verticalalignment='bottom', transform=ax1.transAxes)

        ax1.set_title('%s' % tablename)

        # Ratio subplot with relative pdf uncertainties, denominator in ratio = first order in order_list (not necessarly LO!)
        ax2 = plt.subplot(gs[2, :])
        ax2.set_xlim([xmin, xmax])

        ordernames = ''
        order_index = 0
        for order_item in order_list:
            ax2.semilogx(x_axis, xs_chosen[0, order_index, :]/xs_chosen[0, 0, :], '.', ms=4, ls='dashed', linewidth=0.2,
                         color=_order_color[order_item], label=order_item)
            # divide xs fo each order by xs of first given order! (--> ratio)
            ax2.fill_between(x_axis, (xs_chosen[0, order_index, :]/xs_chosen[0, 0, :])+rel_pdf_unc[0, order_index, 2, :],
                             (xs_chosen[0, order_index, :]/xs_chosen[0, 0, :])+rel_pdf_unc[0, order_index, 1, :], color=_order_color[order_item],
                             alpha=0.46)  # , hatch=_order_hatch[order_item])

            order_index += 1
            ordernames += '_%s' % order_item

        ax2.set_ylabel('rel. pdf_unc')
        ax2.text(0.02, 0.86, 'Reference order: %s' %
                 order_list[0], horizontalalignment='left', verticalalignment='bottom', transform=ax2.transAxes)
        ax2.axhline(y=1, xmin=0, xmax=1, linewidth=0.4,
                    color='k', linestyle='dotted')

        fig.tight_layout()

        if given_filename is not None:
            filename = '%s.pdfunc-%s.plot%s.%s%s' % (
                given_filename, pdfsets[0], plotnr, dimspeci, ordernames[1:])
        else:
            filename = '%s.pdfunc-%s.plot%s.%s%s' % (
                tablename, pdfsets[0], plotnr, dimspeci, ordernames[1:])

        fig.savefig('%s.png' % filename)
        print '[fastnnlo_pdfunc]: Saved as: %s.png' % filename


######################################################################################################################################

def main():
    # np.seterr(over='raise') #for debugging
    # just for measuring wall clock time - not necessary
    start_time = timeit.default_timer()
    # Input and Options
    parser = argparse.ArgumentParser(epilog='')

    # table is always required.
    parser.add_argument('table', type=str,
                        help='fastNLO table that shall be evaluated.')
    parser.add_argument('-p', '--pdfset', default=['CT14nlo'], type=str, nargs='*',
                        help='PDFset(s) to evaluate fastNLO table. First given PDFset will be used as reference PDF in ratio subplot.')
    parser.add_argument('-m', '--member', default=0, type=int,
                        help='Member of PDFset, default is 0.')
    parser.add_argument('-o', '--order', required=False, nargs='+', type=str, choices=['LO', 'NLO', 'NNLO'],  # should 0, 1, 2 also be made possible?
                        help='Order up to which pdf uncertainty will be plotted. Valid choices: LO, NLO, NNLO. '
                        'Per default: Plot all orders that are available in given table.')
    parser.add_argument('-f', '--filename', default=None, type=str,
                        help='Output filename (optional).')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Increase output verbosity.')

    # Parse arguments
    args = vars(parser.parse_args())

    # Table name
    tablename = os.path.splitext(os.path.basename(args['table']))[0]
    # to get rid of extension (.tab.gz or .tab)
    tablename = os.path.splitext(tablename)[0]
    print '\n'
    print '[fastnnlo_pdfunc]: Analysing table:      ', tablename

    # PDF set name
    pdfsets = []
    for ipdf in args['pdfset']:
        pdfsets.append(os.path.basename(ipdf))
    print '[fastnnlo_pdfunc]: Using PDF set(s):     ', pdfsets

    # Orders that will be plotted
    iorders = []
    iordmin = 0
    iordmax = _text_to_order['NNLO']
    if args['order'] is None:
        print '[fastnnlo_pdfunc]: Evaluate table up to highest available order.'
    else:
        for order in args['order']:
            iorders.append(_text_to_order[order])
        iordmin = min(iorders)
        iordmax = max(iorders)
        print '[fastnnlo_pdfunc]: Evaluate table in the following order(s): ', args[
            'order']

    # Given output filename
    given_filename = args['filename']

    # Verbosity
    verb = args['verbose']

    ######################################## Start EVALUATION with fastNLO library #######################################
    # Take the general information (bin_bounds, labels, order_existence, etc.) from first given pdfset.
    fnlo = fastnlo.fastNLOLHAPDF(args['table'], pdfsets[0], args['member'])

    # Get labeling, dimensionality, bin_bounds of table etc. via dimensionhandler()
    ndim, binningdim, xlabel, all_bin_bounds, x_axis_list, xmin_list, xmax_list, x_errors, nplot, dimlabels, tmp_len_list, bbDimO, bbDimM = dimensionhandler(
        fnlo, verb)

    # Check existence of orders in table
    # dictionary that contains the input orders (int) and whether they exist or not (bool)
    _order_existence = {}
    order_list = []  # list that will contain strings of orders to be investigated
    if args['order'] is not None:
        for i in range(0, len(iorders)):  # check existence for all input orders
            # append chosen order to dictionary for order existence (for instance: '0':'True' for existing LO)
            # remember that the input orders are unsorted and could therefore as well be [NNLO, LO]!
            _order_existence[iorders[i]] = fnlo.SetContributionON(
                fastnlo.kFixedOrder, iorders[i], True)
            if _order_existence[iorders[i]] == False:
                sys.exit('[fastnnlo_pdfunc]: Chosen order %s is not available in given table. Aborted! \n' %
                         _order_to_text[iorders[i]])
            else:
                order_list.append(_order_to_text[iorders[i]])

    elif args['order'] is None:
        for i in [0, 1, 2]:  # assuming we have no NNNLO tables yet, could of course be kept more general
            _order_existence[i] = fnlo.SetContributionON(
                fastnlo.kFixedOrder, i, True)
            if _order_existence[i] == False:
                print '[fastnnlo_pdfunc]: Given table does not contain %s. \n' % _order_to_text[i]
            elif _order_existence[i] == True:
                if verb:
                    print '[fastnnlo_pdfunc]: Given table contains %s. \n' % _order_to_text[i]
                order_list.append(_order_to_text[i])
                iordmax = i  # update highest existing order

    print '[fastnnlo_pdfunc]: List of orders that will be investigated: ', order_list, '\n \n'

    # Start individual evaluation for each of the given PDF sets
    # will contain for each given pdf the total cross section for LO, NLO, NNLO (or requested orders)
    xs_list = []
    # list for relative pdf uncertainties (low, high) for each PDF in LO, NLO, NNLO
    rel_unc_list = []

    # Outermost loop: go through all the given pdf sets
    # Evaluate fastNLO table focusing on PDF uncertainties
    for pdf in pdfsets:
        xs_list_tmp = []  # temporary xs list for single pdf
        rel_unc_list_tmp = []  # only temporary list. will be renewed for each pdf
        print '-----------------------------------------------------------------------------------------------'
        print '#############################  %s  ########################################' % pdf
        print '-----------------------------------------------------------------------------------------------'
        print '[fastnnlo_pdfunc]: Calculate XS and uncertainty for %s \n' % pdf
        fnlo = fastnlo.fastNLOLHAPDF(args['table'], pdf, args['member'])

        # in earlier versions: calculated cross sections for ALL available orders (not only the selected ones!) and filled array xs_all
        # in this version: only calculate the cross section for the chosen orders, produce array xs_chosen.
        # This can cause the following situation: if order_list=[LO, NNLO], the nlo cross section won't be calculated at all!
        # Keep this in mind, for instance when printing xs_chosen!
        for n in order_list:
            for j in range(0, iordmax+1):
                if j <= _text_to_order[n]:
                    fnlo.SetContributionON(fastnlo.kFixedOrder, j, True)
                else:
                    fnlo.SetContributionON(fastnlo.kFixedOrder, j, False)
            if verb:
                print '[fastnnlo_pdfunc]: Calculate XS for order: %s' % n
                print '[fastnnlo_pdfunc]: ---- ---- ---- ---- ---- ---- ---- \n'

            fnlo.CalcCrossSection()
            xs_list_tmp.append(fnlo.GetCrossSection())

            ### Get PDF uncertainties ###
            # RELATIVE pdf uncertainty
            rel_pdf_unc_item = np.array(fnlo.GetPDFUncertaintyVec(
                fastnlo.kLHAPDF6))  # now calculated for order n
            # rel_pdf_unc_item = np.ones((3,33))                  #####this is just for testing the plotting function without big calculation
            # rel_pdf_unc_item[2,:]=(-1)*rel_pdf_unc_item[2,:]    #####only for testing
            rel_unc_list_tmp.append(rel_pdf_unc_item)
            if verb:
                print '[fastnnlo_pdfunc]: Cross section in %s: \n' % n, np.array(
                    xs_list_tmp)[-1], '\n'  # print most recent xs
                print '[fastnnlo_pdfunc]: Relative PDF uncertainty in %s: \n' % n
                # has 3 entries: central value (xs), unc_high, unc_low
                print rel_pdf_unc_item, '\n'
                print '-----------------------------------------------------------------------------------------------'
                print '-----------------------------------------------------------------------------------------------'
        xs_list.append(xs_list_tmp)
        rel_unc_list.append(rel_unc_list_tmp)
    xs_chosen = np.array(xs_list)
    # Remember: both arrays here do only contain the CHOSEN orders! (or per default LO, NLO, (NNLO))
    rel_pdf_unc = np.array(rel_unc_list)
    ############### Structure is the same in case of multidimensional tables, the array is just longer and this will be taken care of when plotting #######
    #########                                                                                               #########
    # structure of rel_pdf_unc:                                                                                     #
    # rel_pdf_unc[0,:,:] refers to the order stored in order_list[0] --> this (0th item) is not necessarily LO!     #
    # The same applies to rel_pdf_unc[1,:,:] and so on.                                                             #
    # rel_pdf_unc[0,0,:] corresponds to the cross section in the order stored in order_list[0].                     #
    # rel_pdf_unc[0,1,:] corresponds to the rel. uncertainty upwards (in the order order_list[0])                   #
    # rel_pdf_unc[0,2,:] corresponds to the rel. uncertainty downwards (in the order order_list[0])                 #
    ########                                                                                                #########

    if verb:
        print '[fastnnlo_pdfunc]: Cross section summary. '
        print '[fastnnlo_pdfunc]: For each pdf %s xs_chosen contains XS corresponding to %s.' % (
            pdfsets, order_list)
        print xs_chosen, '\n \n'
        # print 'Size xs_chosen: ', np.shape(xs_chosen) ###test
        # print 'Size rel_pdf_unc: ', np.shape(rel_pdf_unc) ###test

    # ABSOLUTE pdf uncertainty
    # length of axis 0 in xs_chosen equals number of (investigated) pdfsets
    num_pdfsets = np.size(xs_chosen, 0)
    # length of axis 1 in xs_chosen equals number of (investigated) orders
    num_orders = np.size(xs_chosen, 1)
    # abs_pdf_unc = np.empty([num_pdfsets, num_orders, 2, len(x_axis)]) #old version for 1dim tables
    abs_pdf_unc = np.empty([num_pdfsets, num_orders, 2, np.size(xs_chosen, 2)])

    for p in range(0, num_pdfsets):
        for k in range(0, num_orders):  # k is order index
            # absolute uncertainties downwards (low)
            abs_pdf_unc[p, k, 0, :] = rel_pdf_unc[p,
                                                  k, 2, :]*xs_chosen[p, k, :]
            abs_pdf_unc[p, k, 1, :] = rel_pdf_unc[p, k, 1, :] * \
                xs_chosen[p, k, :]  # absolute uncertainties upwards (high)

    if verb:
        print '[fastnnlo_pdfunc]: Absolute PDF uncertainties downwards, upwards for %s in %s: \n' % (
            pdfsets, order_list)
        print abs_pdf_unc, '\n'

    ################################## Start PLOTTING of the pdf uncertainties ############################################################################

    # Plot all required plots by taking the relevant information from the big/full arrays that contain the results for all(!) bins of the investigated table.
    # This is done by using the indices obtained with the dimensionhandler() method.

    low_index = 0
    up_index = 0
    #print 'tmp_len_list: ', tmp_len_list
    for n in range(0, nplot):
        # tmp_len_list contains the last index (full arrays) for each bin of outer(!) dimension
        up_index = tmp_len_list[n]
        #(= last index for each x-axis); originally used in dimensionhandler()
        xmin = xmin_list[n]
        xmax = xmax_list[n]
        x_axis = x_axis_list[n]

        #print 'lower index: ', low_index
        #print 'upper index: ', up_index

        xs_chosen_n = xs_chosen[:, :, low_index:(up_index+1)]
        abs_pdf_unc_n = abs_pdf_unc[:, :, :, low_index:(up_index+1)]
        rel_pdf_unc_n = rel_pdf_unc[:, :, :, low_index:(up_index+1)]

        if (ndim >= 2):
            bbDimO_low = bbDimO[low_index, 0]
            bbDimO_up = bbDimO[up_index, 1]
            print 'bin bounds DimO: low=%s, up=%s. ' % (bbDimO_low, bbDimO_up)

            if (ndim == 3):
                bbDimM_low = bbDimM[low_index, 0]
                bbDimM_up = bbDimM[up_index, 1]
                print 'bin bounds DimM: low=%s, up=%s. ' % (
                    bbDimM_low, bbDimM_up)

        # outputs for debugging
        #print 'plot %s'%n
        #print 'xmin, xmax, x_axis, xs_chosen, abspdf, relpdf, lowind, upind \n'
        #print xmin
        #print xmax
        #print x_axis
        #print xs_chosen_n, 'len', len(xs_chosen_n[0, 0, :])
        #print abs_pdf_unc_n
        #print rel_pdf_unc_n
        #print low_index
        #print up_index
        #print tmp_len_list[n]
        #... and so on --> then use plotting function

        low_index = up_index+1

        # so far, plots are simply numbered --> see label with observable range / phase space information in plot itself
        plotnr = '%s' % n

        # the following is used to name the plots
        if (ndim > 1):
            if (ndim == 2):
                # name and range of outermost binning dimension
                nameDimO = fnlo.GetDimLabel(
                    0)+': %s to %s' % (bbDimO_low, bbDimO_up)
                nameDimM = ''  # no middle dimension
                # name of innermost dimension (== xlabel)
                nameDimI = fnlo.GetDimLabel(1)

            elif (ndim == 3):
                # name of outermost binning dimension
                nameDimO = fnlo.GetDimLabel(
                    0)+': %s to %s' % (bbDimO_low, bbDimO_up)
                # name of middle binning dimension
                nameDimM = fnlo.GetDimLabel(
                    1)+': %s to %s' % (bbDimM_low, bbDimM_up)
                # name of innermost dimension (== xlabel)
                nameDimI = fnlo.GetDimLabel(2)
        else:
            nameDimO = ''
            nameDimM = ''
            nameDimI = fnlo.GetDimLabel(0)  # same as xlabel

        plotting(pdfsets, order_list, xs_chosen_n, abs_pdf_unc_n, rel_pdf_unc_n, xlabel, tablename, given_filename,
                 x_axis, xmin, xmax, nameDimO, nameDimM, nameDimI, dimlabels, plotnr)

    stop_time = timeit.default_timer()
    timediff = stop_time-start_time
    print 'Time: %s sec = %s min' % (timediff, round(timediff/60., 2))


if __name__ == '__main__':
    main()
