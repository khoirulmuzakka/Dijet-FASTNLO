#!/usr/bin/env python2
#-*- coding:utf-8 -*-

##############################################
#
# Plotting of the scale uncertainty.
#
#
# Created by B.Schillinger, 09.10.2018
# Modified by K. Rabbertz, 31.10.2018
# Last modified: 01.11.2018
#
#############################################

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.pyplot import cm
import matplotlib.pylab as pylab
import matplotlib.ticker #to fix x-axis ticks
import os, sys

import fastnlo
from fastnlo import SetGlobalVerbosity

_text_to_order = {'LO':0, 'NLO':1, 'NNLO':2}
_order_to_text = {0:'LO', 1:'NLO', 2:'NNLO'}
_order_color   = {'LO':'g', 'NLO':'b', 'NNLO':'r'}


########################################################################################################################


# Function for plotting a list of orders into one figure
# The ratio is done always with respect to the first order appearing in the order_list
# given e.g. via '-o NLO NNLO', i.e. the first evaluated cross section, here NLO,
# that is stored in xs_all[0].
def plotting(x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_name):
        if variation_type=='Scale uncertainty (2P)':
                vartype='2P'
        elif variation_type=='Scale uncertainty (6P)':
                vartype='6P'

        gs = gridspec.GridSpec(3,3)
        fig = plt.figure(figsize=(7,7))
        ax1 = plt.subplot(gs[:-1,:])

        # For plotting various results, max = 3, 'next to each other', handling via shift from bincenter
        if len(order_list)==1:
                shift_list=[1.00]
        elif len(order_list)==2:
                shift_list=[0.98, 1.02]
        elif len(order_list)==3:
                shift_list=[0.98, 1.00, 1.02]
        else:
                print '[fastnnlo_scaleunc]: List of items to plot too long. Aborted!'
                print '[fastnnlo_scaleunc]: The order list is', order_list
                exit(1)

        # Plot all x section results in order_list
        xs_index = -1
        for order_item, shift in zip(order_list, shift_list):
                xs_index += 1
                ax1.errorbar(x_axis*shift, xs_all[xs_index], yerr=abs(abs_scale_unc[xs_index]), elinewidth=1, linewidth=0.0, ms=4, color=_order_color[order_item], fmt='.', label=order_item)

        ax1.set_xlim([xmin, xmax])
        ax1.set_xscale('log', nonposx='clip')
        ax1.set_yscale('log', nonposy='clip')
        ax1.set_xlabel('%s' %xlabel)
        ax1.set_ylabel('XS with abs_scale_unc', rotation=90)
        ax1.legend(fontsize=10, numpoints=1)
        ax1.text(0.02, 0.10, '%s' %variation_type, horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)
        ax1.text(0.02, 0.05, '%s' %scale_name, horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)
        ax1.set_title('%s' %tablename)


        # Do subplot with ratios and relative scale uncertainties; denominator in ratio = first order in order_list
        ax2 = plt.subplot(gs[2, :])
        ax2.set_xlim([xmin, xmax])
        xs_index = -1
        ordernames = ''
        for item in order_list:
                xs_index   += 1
                ordernames += '_%s' %item
                ax2.semilogx(x_axis, xs_all[xs_index]/xs_all[0], ls='dashed', lw=1.0, color=_order_color[item], label=item)
                ax2.fill_between(x_axis, (xs_all[xs_index]/xs_all[0])+rel_scale_unc[xs_index, 2, :],
                                 (xs_all[xs_index]/xs_all[0])+rel_scale_unc[xs_index, 1, :], color=_order_color[item], alpha=0.50)
                ax2.set_ylabel('rel. scale unc')

        ax2.axhline(y=1, xmin=0, xmax=1, color='k', linestyle='dotted', linewidth=1.6, alpha=0.2)
        fig.tight_layout()

        if given_filename is not None:
                filename = '%s.scaleunc-%s.%s.%s' %(given_filename, vartype, ordernames[1:], scale_name)
        else:
                filename = '%s.scaleunc-%s.%s.%s' %(tablename, vartype, ordernames[1:], scale_name)

        fig.savefig('%s.png' %filename)
        print '[fastnnlo_scaleunc]: Plot saved as: %s.png' %filename
        fig.savefig('%s.eps' %filename)
        print '[fastnnlo_scaleunc]: Plot saved as: %s.eps' %filename


########################################################################################################################


def main():
        # Define arguments & options
        parser = argparse.ArgumentParser(epilog='')

        parser.add_argument('table', type=str, help='fastNLO table that shall be evaluated. To avoid mixing up the table argument with options taking multiple inputs like -o|--order, table must come first or must be separated from the -o option by another option!') #table is always required

        parser.add_argument('-p', '--pdfset', default='CT14nlo',
                                help='PDFset to evaluate fastNLO table.')
        parser.add_argument('-m', '--member', default=0, type=int,
                                help='Member of PDFset, default is 0.')
        parser.add_argument('-o', '--order', required=False, nargs='+', type=str, choices=['LO', 'NLO', 'NNLO'],
                            help='Blank-separated list of orders up to which the scale uncertainty is shown: LO, NLO, and/or NNLO. If nothing is chosen, successively produce plots for all orders that are available in table.')
        parser.add_argument('-s', '--scale', default=1, required=False, nargs='?', type=int, choices=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],
                                help='For flexible-scale tables define central scale choice for MuR and MuF by selection enum fastNLO::ScaleFunctionalForm ("1"=kScale1, "2"=kScale2, "3"=kQuadraticSum).')
        parser.add_argument('-a', '--asymmetric', action="store_true",
                            help='If -a is chosen, use asymmetric (6P) scale variations; otherwise use symmetric ones (2P).')
        parser.add_argument('-f', '--filename', default=None, type=str,
                                help='Output filename (optional).')
        parser.add_argument('-v', '--verbose', action="store_true",
                            help="Increase output verbosity")

        # Parse arguments
        args = vars(parser.parse_args())

        # Table name
        tablename = os.path.splitext(os.path.basename(args['table']))[0]
        tablename = os.path.splitext(tablename)[0] #to get rid of extension (.tab.gz or .tab)
        print '\n'
        print '[fastnnlo_scaleunc]: Analysing table', tablename

        # PDF set name
        pdfset = os.path.basename(args['pdfset'])
        print '[fastnnlo_scaleunc]: Using PDF set', pdfset

        # Orders to be shown
        iorders = []
        iordmin = 0
        iordmax = _text_to_order['NNLO']
        if args['order'] is None:
                print '[fastnnlo_scaleunc]: Evaluate table up to highest available order.'
        else:
                for ord in args['order']:
                        iorders.append(_text_to_order[ord])
                iordmin = min(iorders)
                iordmax = max(iorders)
                print '[fastnnlo_scaleunc]: Evaluate table up to order(s)', args['order']

        # Type of scale variation (symmetric vs asymmetric)
        if args['asymmetric']:
                scale_var_type = fastnlo.kAsymmetricSixPoint
                variation_type = 'Scale uncertainty (6P)'
        else:
                scale_var_type = fastnlo.kSymmetricTwoPoint
                variation_type = 'Scale uncertainty (2P)'

        # Scale choice
        scale_choice = args['scale']

        # Given filename
        given_filename = args['filename']

        # Verbosity
        verb = args['verbose']

        ###################### Start EVALUATION with fastNLO library ###################################################
        # SetGlobalVerbosity(0) # Does not work since changed to default in the following call
        fnlo = fastnlo.fastNLOLHAPDF(args['table'], args['pdfset'], args['member'])

        # Get labeling for the x-axis
        # Dimensionality of the table:
        ndim = fnlo.GetNumDiffBin()
        if verb:
                print '[fastnnlo_scaleunc]: Dimensions:', ndim

        # Labels of all the dimensions:
        labels = fnlo.GetDimLabels()
        if verb:
                print '[fastnnlo_scaleunc]: Labels:', labels

        # Label of first dimension:
        xlabel = fnlo.GetDimLabel(0)
        if verb:
                print '[fastnnlo_scaleunc]: x-label:', xlabel

        # Creating x-axis
        bin_bounds = np.array(fnlo.GetObsBinsBounds(0))
        if verb:
                print '[fastnnlo_scaleunc]: bin_bounds.T: \n', bin_bounds.T, '\n'
                print '[fastnnlo_scaleunc]: bin_bounds.flatten()', bin_bounds.flatten(), '\n'

        x_axis = (bin_bounds.T[0]+bin_bounds.T[1])/2. #this is a list of bin centers
        xmin = 0.95*min(bin_bounds.ravel())
        xmax = 1.05*max(bin_bounds.ravel())
        if verb:
                print '[fastnnlo_scaleunc]: xmin =', xmin, ', xmax =', xmax

        # Preparing x-errors (via bin_bounds) --> x_errors[0, :] are initially negative (positive via -1*), x_errors[1, :] positive
        x_errors = np.array([-1*(bin_bounds.T[0]-x_axis), bin_bounds.T[1]-x_axis])
        if verb:
                print '[fastnnlo_scaleunc]: \n x_errors: ', x_errors, '\n'

        # Check existence of orders in table
        lflex = fnlo.GetIsFlexibleScaleTable()
        scale_name = 'scale1'
        o_existence = [False, False, False]
        cont_order = -1
        max_order  = 0
        for i in [0, 1, 2]:
                o_existence[i] = fnlo.SetContributionON(fastnlo.kFixedOrder, i, True)
                if o_existence[i]:
                        max_order   = i
                        if scale_choice == 1:
                                scale_name = fnlo.GetScaleDescription(i,0)
                        elif lflex:
                                scale_name = fnlo.GetScaleDescription(i,1)
                        if cont_order == i-1:
                                cont_order += 1
        if verb:
                print '[fastnnlo_scaleunc]: Table has continuous orders up to', cont_order, 'and a maximal order of', max_order

        if iordmax > cont_order:
                print '[fastnnlo_scaleunc]: Invalid choice of orders. Aborted!'
                print '[fastnnlo_scaleunc]: Highest order requested is', _order_to_text[iordmax], 'but orders are available only up to', cont_order
                exit(1)

        order_list = []
        if args['order'] is None:
                for iord in range(cont_order+1):
                        order_list.append(_order_to_text[iord])
        else:
                order_list = args['order']

        print '[fastnnlo_scaleunc]: List of requested orders:', order_list

        # For flexible-scale tables set scale to user choice (default is 1)

        if lflex:
                print '[fastnnlo_scaleunc]: Setting requested scale choice for flexible-scale table:', scale_choice
                fnlo.SetMuRFunctionalForm(scale_choice-1)
                fnlo.SetMuFFunctionalForm(scale_choice-1)
        else:
                if scale_choice == 1:
                        print '[fastnnlo_scaleunc]: Evaluating fixed-scale table. Scale choice must be', scale_choice
                else:
                        print '[fastnnlo_scaleunc]: No scale choice possible for fixed-scale table. Aborted!'

        # Now evaluate fastNLO table having a look at scale uncertainties
        xs_list = [] #will contain total cross section for LO, NLO, NNLO (note: could also be handled via GetScaleUnertaintyVec()[0] )
        rel_unc_list = [] #list for relative scale uncertainties (low, high) for LO, NLO, NNLO
        for n in order_list:
                for j in range(0, max_order+1):
                        if j <= _text_to_order[n]:
                                fnlo.SetContributionON(fastnlo.kFixedOrder, j, True)
                        else:
                                fnlo.SetContributionON(fastnlo.kFixedOrder, j, False)
                if verb:
                        print '[fastnnlo_scaleunc]: \n'
                        print '[fastnnlo_scaleunc]: Calculate XS for order: %s' %n, '\n'
                        print '[fastnnlo_scaleunc]: ----  ----  ----  ----  ----  ----  ----  ----'
                fnlo.CalcCrossSection()
                xs_list.append(fnlo.GetCrossSection())

                ### Get scale uncertainties ###
                if verb:
                        print '[fastnnlo_scaleunc]: Used scale factor MuF: ', fnlo.GetScaleFactorMuF()
                        print '[fastnnlo_scaleunc]: Used scale factor MuR: ', fnlo.GetScaleFactorMuR(), '\n'
                        print '[fastnnlo_scaleunc]: Calculate scale uncertainties \n'
                ## RELATIVE scale uncertainty with chosen type of scale variation (symmetric or asymmetric)
                rel_scale_unc_item = np.array(fnlo.GetScaleUncertaintyVec(scale_var_type)) #calculate this already for all accessible orders in any case
                rel_unc_list.append(rel_scale_unc_item)
                if verb:
                        print '[fastnnlo_scaleunc]: \n'
                        print '[fastnnlo_scaleunc]: Relative scale uncertainty in %s: \n'%n
                        print rel_scale_unc_item, '\n' #3 entries: central value (xs), unc_low, unc_high
                        print '---------------------------------------------------------------------------------------'
                        print '---------------------------------------------------------------------------------------'

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
                print '[fastnnlo_scaleunc]: Cross section xs_all uses ', order_list

        if verb:
                print xs_all, '\n \n'

        ## ABSOLUTE scale uncertainty
        num_orders = np.size(xs_all, 0) #length of axis 0 in xs_all equals number of orders
        abs_scale_unc = np.empty([num_orders, 2, len(x_axis)])
        for k in range(0, len(xs_all)):
                abs_scale_unc[k, 0, :] = rel_scale_unc[k, 2, :]*xs_all[k] #absolute uncertainties downwards (low)
                abs_scale_unc[k, 1, :] = rel_scale_unc[k, 1, :]*xs_all[k] #absolute uncertainties upwards (high)

        if verb:
                print '[fastnnlo_scaleunc]: Absolute Scale uncertainties downwards, upwards (order by order): \n'
                print abs_scale_unc, '\n'

        ############################## Do the plotting ####################################################

        plotting(x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_name)

if __name__ == '__main__':
        main()
