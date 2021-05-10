#!/usr/bin/env python2
#-*- coding:utf-8 -*-

###############################################
#
# Comparison of different values for alpha_s
#
#
# Created by B.Schillinger, 12.12.2018
# Last modified: 21.12.2018
#
#
###############################################

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


####################################################################################


def plotting(alphas_list, xs_chosen, order_list, x_axis, xlabel, tablename, pdfset, xmin, xmax, given_filename):

    gs = gridspec.GridSpec(3, 3)
    fig = plt.figure(figsize=(7, 7))
    ax1 = plt.subplot(gs[:-1, :])

    # XS for several values for alpha_s will be plotted 'next to each other',
    # need tiny shift from bincenter to display these.
    if len(alphas_list) == 1:
        shift_list = [1.00]
    elif len(alphas_list) == 2:
        shift_list = [0.91, 1.01]
    elif len(alphas_list) == 3:
        shift_list = [0.98, 1.00, 1.02]
    elif len(alphas_list) == 4:
        shift_list = [0.97, 0.99, 1.01, 1.03]
    elif len(alphas_list) == 5:
        shift_list = [0.96, 0.98, 1.00, 1.02, 1.04]
    else:
        print '[fastnnlo_alphas]: Too many alpha_s values to produce combined plot. Aborted!'
        print '[fastnnlo_alphas]: Max. possible: 5. Number of given alpha_s values: %s' % len(
            alphas_list)
        sys.exit()

    order_index = 0  # necessary, because even though LO might not even be investigated, we need index 0, which could then correspond to NLO or NNLO
    ordernames = ''
    alphanames = ''

    # producing one plot per order (comparing all given alpha_s values)
    for order_item in order_list:
        fig = plt.figure(figsize=(7, 7))
        ax1 = plt.subplot(gs[:-1, :])

        print '[fastnnlo_alphas]: Producing %s plot.' % order_item
        color = iter(cm.rainbow(np.linspace(0, 1, len(alphas_list))))
        #patches = []

        alphas_index = 0  # starting with first given alpha_s
        # go through all given values for alpha_s
        for alpha, shift in zip(alphas_list, shift_list):
            c = next(color)
            ax1.errorbar(x_axis*shift, xs_chosen[alphas_index, order_index, :], elinewidth=1,
                         linewidth=0.0, ms=4, color=c, fmt='.', label=r'$\alpha_s = %s$' % alpha)  # uncertainties?
            alphaname = str(alpha).replace('.', '')  # for filename
            alphanames += '_%s' % alphaname
            alphas_index += 1

        ax1.set_xscale('log', nonposx='clip')
        ax1.set_yscale('log', nonposy='clip')
        ax1.set_xlabel('%s' % xlabel)
        ax1.set_ylabel(r'XS for different $\alpha_s$',
                       fontsize=12, rotation=90)
        ax1.legend(fontsize=10, numpoints=1)
        ax1.text(0.04, 0.08, '%s' % order_item, fontsize=12, horizontalalignment='left',
                 verticalalignment='bottom', transform=ax1.transAxes)
        ax1.text(0.04, 0.02, '%s' % pdfset, fontsize=12, horizontalalignment='left',
                 verticalalignment='bottom', transform=ax1.transAxes)

        ax1.set_title('%s' % tablename)

        ax2 = plt.subplot(gs[2, :])
        ax2.set_xlim([xmin, xmax])
        # ax2.set_ylim([0.97, 1.03]) #for testing purposes

        # Ratioplot normalised to first given alpha_s value
        color = iter(cm.rainbow(np.linspace(0, 1, len(alphas_list))))
        for a in range(0, len(alphas_list)):
            c = next(color)
            # divide xs for each alpha_s by xs of first given alpha_s
            ax2.semilogx(x_axis, xs_chosen[a, order_index, :]/xs_chosen[0, order_index, :], '.', ms=4, ls='dashed', linewidth=0.2,
                         color=c, label=alphas_list[a])

            # and the uncertainties?

        ax2.set_ylabel(
            r'Ratio $\alpha_{s, i}$ to $\alpha_{s, 0}$', fontsize=12, rotation=90)
        ax2.text(0.02, 0.92, r'Reference $\alpha_s$: %s' %
                 alphas_list[0], horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes)
        ax2.axhline(y=1, xmin=0, xmax=1, linewidth=0.4,
                    color='k', linestyle='dotted')

        plt.grid(b=True, which='both', axis='both',
                 color='k', linestyle='steps--', alpha=0.2)

        fig.tight_layout()

        order_index += 1
        ordernames += order_item

        if given_filename is not None:
            filename = '%s.alphas-%s.%s' % (given_filename,
                                            order_item, alphanames[1:])
        else:
            filename = '%s.alphas-%s.%s' % (tablename,
                                            order_item, alphanames[1:])

        fig.savefig('%s.png' % filename)
        print '[fastnnlo_alphas]: %s plot saved as: %s.png' % (
            order_item, filename)
        plt.close()


####################################################################################

def main():

    # just for measuring wall clock time - not necessary
    start_time = timeit.default_timer()

    # Input and Options
    parser = argparse.ArgumentParser(epilog='')

    # table is always required.
    parser.add_argument('table', type=str,
                        help='fastNLO table that shall be evaluated.')
    parser.add_argument('-p', '--pdfset', default='CT14nlo', type=str, nargs='?',
                        help='PDFset to evaluate fastNLO table.')
    parser.add_argument('-m', '--member', default=0, type=int,
                        help='Member of PDFset, default is 0.')
    parser.add_argument('-o', '--order', required=False, nargs='+', type=str, choices=['LO', 'NLO', 'NNLO'],
                        help='Order in which cross section will be calculated. Valid choices: LO, NLO, NNLO. '
                        'Per default: Plot all order that are available in given table.')
    parser.add_argument('-a', '--alphas', default=[0.118], nargs='*', type=float,
                        help='Alpha_s values that will be used.')
    parser.add_argument('-f', '--filename', default=None, type=str,
                        help='Output filename (optional).')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Increase output verbosity.')

    # Parse arguments
    args = vars(parser.parse_args())

    # Table name
    table = args['table']
    if table.endswith('.tab.gz'):
        tablename = table.replace('.tab.gz', '', 1)
    elif table.endswith('.tab'):
        tablename = table.replace('.tab', '', 1)
    else:
        print '[fastnnlo_alphas]: Error! Wrong extension for table: ', table
        exit(1)
    print '\n'
    print '[fastnnlo_alphas]: Analysing table: ', table

    # PDF set name
    pdfset = os.path.basename(args['pdfset'])
    print '[fastnnlo_alphas]: Using PDF set: ', pdfset

    # Orders to be shown
    iorders = []
    iordmin = 0
    iordmax = _text_to_order['NNLO']
    if args['order'] is None:
        print '[fastnnlo_alphas]: Evaluate table up to highest available order.'
    else:
        for order in args['order']:
            iorders.append(_text_to_order[order])
        iordmin = min(iorders)
        iordmax = max(iorders)
        print '[fastnnlo_alphas]: Evaluate table in the following order(s): ', args[
            'order']

    # Alpha_s values to be investigated
    alphas_list = args['alphas']
    print '[fastnnlo_alphas]: Alpha_s values that will be used: ', alphas_list

    # Given output filename
    given_filename = args['filename']

    # Verbosity
    verb = args['verbose']

    ###################################### Start EVALUATION with fastNLO library #################################################
    # Take the general information (bin_bounds, labels, order_existence, etc.) from first given pdfset.
    # use fastNLOALphas() instead of fastNLOLHAPDF() in order to change alpha_s
    fnlo = fastnlo.fastNLOAlphas(table, pdfset, args['member'])

    # Get labeling for the x-axis
    # Dimensionality of the table:
    ndim = fnlo.GetNumDiffBin()
    if verb:
        print '\n'
        print '[fastnnlo_alphas]: Table Dimensions: ', ndim

    # Labels of all the dimensions:
    dimlabels = fnlo.GetDimLabels()
    if verb:
        print '[fastnnlo_alphas]: Dimension Labels: ', dimlabels

    # Label of first dimension:
    xlabel = fnlo.GetDimLabel(0)
    if verb:
        print '[fastnnlo_alphas]: x-label: ', xlabel

    # Creating x-axis
    bin_bounds = np.array(fnlo.GetObsBinsBounds(0))
    if verb:
        print '[fastnnlo_alphas]: bin_bounds.T: \n', bin_bounds.T, '\n'
        print '[fastnnlo_alphas]: bin_bounds.flatten(): \n', bin_bounds.flatten(), '\n'

    x_axis = (bin_bounds.T[0]+bin_bounds.T[1]) / \
        2.  # this is a list of bin centers
    xmin = 0.95*min(bin_bounds.ravel())
    xmax = 1.05*max(bin_bounds.ravel())
    if verb:
        print '[fastnnlo_alphas]: xmin=%s, xmax=%s. \n' % (xmin, xmax)

    # Preparing x-errors (via bin_bounds) --> x_errors[0, :] are initially negative (positive via -1*), x_errors[1, :] positive
    x_errors = np.array([-1*(bin_bounds.T[0]-x_axis), bin_bounds.T[1]-x_axis])
    if verb:
        print '[fastnnlo_alphas]: x_errors: \n', x_errors, '\n'

    # Check existence of orders in table
    # dictionary that contains the input orders (int) and whether they exist or not (bool)
    _order_existence = {}
    order_list = []  # list that will contain strings of orders to be inverstigated
    if args['order'] is not None:
        for i in range(0, len(iorders)):  # check existence for all input orders
            # append chosen order to dictionary for order existence (for instance: '0':True for existing LO)
            # remember that the input orders are unsorted and could therefore as well be [NNLO, LO]!
            _order_existence[iorders[i]] = fnlo.SetContributionON(
                fastnlo.kFixedOrder, iorders[i], True)
            if _order_existence[iorders[i]] == False:
                sys.exit('[fastnnlo_alphas]: Chosen order %s is not available in given table. Aborted! \n' %
                         _order_to_text[iorders[i]])
            else:
                order_list.append(_order_to_text[iorders[i]])

    elif args['order'] is None:
        for i in [0, 1, 2]:  # assuming we have no NNNLO tables yet, could of course be kept more general
            _order_existence[i] = fnlo.SetContributionON(
                fastnlo.kFixedOrder, i, True)
            if _order_existence[i] == False:
                print '[fastnnlo_alphas]: Given table does not contain %s. \n' % _order_to_text[i]
            elif _order_existence[i] == True:
                if verb:
                    print '[fastnnlo_alphas]: Given table contains %s. \n' % _order_to_text[i]
                order_list.append(_order_to_text[i])
                iordmax = i  # update highest existing order

    print '[fastnnlo_alphas]: List of orders that will be investigated: ', order_list, '\n \n'
    #print fnlo.Print(0)

    # Start individual evaluation for different values for alpha_s
    # will contain for each given alpha_s value the corresponding cross section in requested (or available) orders
    xs_list = []

    # Outermost loop: go through all the given alpha_s values
    # Evaluate fastNLO table accordingly
    for alpha in alphas_list:
        xs_list_tmp = []  # temporary xs list for one alpha_s value
        fnlo.SetAlphasMz(alpha, True)
        print '-----------------------------------------------------------------------------------------'
        print '################################ alpha_s = %s ################################' % alpha
        print '-----------------------------------------------------------------------------------------'
        print '[fastnnlo_alphas]: Alpha_s has now been set to alpha_s=%s' % fnlo.GetAlphasMz()
        print '[fastnnlo_alphas]: Calculate XS for alpha_s = %s \n' % alpha

        # calculate XS only for chosen orders --> xs_chosen
        for n in order_list:
            for j in range(0, iordmax+1):
                if j <= _text_to_order[n]:
                    fnlo.SetContributionON(fastnlo.kFixedOrder, j, True)
                else:
                    fnlo.SetContributionON(fastnlo.kFixedOrder, j, False)
            if verb:
                print '[fastnnlo_alphas]: Calculate XS for order: %s' % n
                print '[fastnnlo_alphas]: ---- ---- ---- ---- ---- ---- ---- \n'

            fnlo.CalcCrossSection()
            xs_list_tmp.append(fnlo.GetCrossSection())

            if verb:
                print '[fastnnlo_alphas]: Cross section in %s: \n' % n, np.array(
                    xs_list_tmp)[-1], '\n'  # print most recent xs

        xs_list.append(xs_list_tmp)
    xs_chosen = np.array(xs_list)

    if verb:
        print '[fastnnlo_alphas]: Cross section summary. '
        print '[fastnnlo_alphas]: For each alpha_s %s xs_chosen contains XS corresponding to %s.' % (
            alphas_list, order_list)
        print xs_chosen, '\n \n'

    ######################################################## Start PLOTTING ######################################################################

    plotting(alphas_list, xs_chosen, order_list, x_axis, xlabel,
             tablename, pdfset, xmin, xmax, given_filename)

    stop_time = timeit.default_timer()
    timediff = stop_time-start_time
    print 'Time: %s sec = %s min' % (timediff, round(timediff/60., 2))


if __name__ == '__main__':
    main()
