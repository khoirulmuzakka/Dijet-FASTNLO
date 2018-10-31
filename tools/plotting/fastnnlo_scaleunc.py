#!/usr/bin/env python2
#-*- coding:utf-8 -*-

##############################################
#
# Plotting of the scale uncertainty.
#
#
#
# Created by B.Schillinger, 09.10.2018
# Last modified: 26.10.2018
#
#############################################


import argparse
import numpy as np
import os, sys
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.pyplot import cm
import matplotlib.pylab as pylab
import matplotlib.ticker #to fix x-axis ticks

import fastnlo

_text_to_order = {'LO':0, 'NLO':1, 'NNLO':2}
_order_to_text = {0:'LO', 1:'NLO', 2:'NNLO'}
_order_color   = {'LO':'g', 'NLO':'b', 'NNLO':'r'}

#function for plotting uncertainties for one single order
def plotting_single(order_index, x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_choice):
        scale_name={0:'kScale1', 1:'kScale2'}
        if variation_type=='scale uncertainty (2P)':
                vartype='2P'
        elif variation_type=='scale uncertainty (6P)':
                vartype='6P'

        gs = gridspec.GridSpec(3,3)
        fig = plt.figure(figsize=(7,7))
        ax1 = plt.subplot(gs[:-1,:])

        order_name = _order_to_text[order_index]
        #ax1.loglog(x_axis, xs_all[order_index], 'gd', alpha=0.1)
        ax1.errorbar(x_axis, xs_all[order_index], #xerr=x_errors,
                        yerr=abs(abs_scale_unc[order_index]), elinewidth=1, linewidth=0.0, ms=4, color=_order_color[order_name], fmt='.', alpha=0.9, label=order_name)

        #ax1.axis()
        ax1.set_xscale('log', nonposx='clip')
        ax1.set_yscale('log', nonposy='clip')
        ax1.set_xlim([xmin, xmax])
        ax1.set_xlabel('%s' %xlabel)
        ax1.set_ylabel('XS with abs_scale_unc', rotation=90)
        ax1.legend(fontsize=10, numpoints=1)
        ax1.text(0.02, 0.10, '%s' %variation_type, horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)
        ax1.text(0.02, 0.05, '%s' %scale_name[scale_choice], horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)

        ax1.set_title('%s' %tablename)

        #subplot with relative scale uncertainties
        ax2 = plt.subplot(gs[2, :])
        ax2.set_xlim([xmin, xmax])

        ax2.semilogx(x_axis, xs_all[order_index]/xs_all[order_index], ls='solid', linewidth=0.02, color='red', alpha=0.8) #is this necessary?

        ax2.fill_between(x_axis, 1+rel_scale_unc[order_index, 2, :], 1+rel_scale_unc[order_index, 1, :], color=_order_color[order_name], alpha=0.6, label=order_name)
        ax2.axhline(y=1, xmin=0, xmax=xmax, color='k', linestyle='dotted')
        #ax2.axis([xmin, xmax, 0.4, 1.6])
        # more flexible axis settings
        ax2.axis([xmin, xmax, 1+1.1*min(rel_scale_unc[order_index, 2, :]), 1+1.1*max(rel_scale_unc[order_index, 1, :])])
        ax2.set_ylabel('rel. scale_unc')
        # .................. #
        plt.tight_layout()

        if given_filename is not None:
                filename = '%s.scaleunc-%s.%s' %(given_filename, vartype, order_name)
        else:
                filename = '%s.scaleunc-%s.%s' %(tablename, vartype, order_name)
        plt.savefig('%s.png' %filename)
        print '[fastnnlo_scaleunc]: saved as: %s.png' %filename

#function for plotting multiple orders into one figure
def plotting_multiple(lowest_order, x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_choice):
        scale_name={0:'kScale1', 1:'kScale2'}
        if variation_type=='scale uncertainty (2P)':
                vartype='2P'
        elif variation_type=='scale uncertainty (6P)':
                vartype='6P'

        gs = gridspec.GridSpec(3,3)
        fig = plt.figure(figsize=(7,7))
        ax1 = plt.subplot(gs[:-1,:])

        #for plotting the different orders 'next to each other', handling via shift from bincenter
        if len(order_list)==2:
                shift_list=[0.98, 1.02]
        elif len(order_list)==3:
                shift_list=[0.98, 1, 1.02]
        else:
                shift_list=np.zeros(len(order_list)) #just in case, usually this should not happen
        xs_index = -1
        for order_item, shift in zip(order_list, shift_list):
                xs_index += 1
                ax1.errorbar(x_axis*shift, xs_all[xs_index], yerr=abs(abs_scale_unc[xs_index]), elinewidth=1, linewidth=0.0, ms=4, color=_order_color[order_item], fmt='.', label=order_item)

        ax1.set_xscale('log', nonposx='clip')
        ax1.set_yscale('log', nonposy='clip')
        ax1.set_xlim([xmin, xmax])
        ax1.set_xlabel('%s' %xlabel)
        ax1.set_ylabel('XS with abs_scale_unc', rotation=90)
        ax1.legend(fontsize=10, numpoints=1)
        ax1.text(0.02, 0.10, '%s' %variation_type, horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)
        ax1.text(0.02, 0.05, '%s' %scale_name[scale_choice], horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)


        ax1.set_title('%s' %tablename)


        #subplot with relative scale uncertainties, denominator in ratio = lowest given order
        ax2 = plt.subplot(gs[2, :])
        ax2.set_xlim([xmin, xmax])

        ordernames = ''
        xs_index = -1
        for order_item in order_list:
                xs_index += 1
                order_index = _text_to_order[order_item]
                ax2.semilogx(x_axis, xs_all[xs_index]/xs_all[lowest_order], '.', ms=4, ls='solid', linewidth=0.2, color=_order_color[order_item], label=order_item)
                #ax2.fill_between(xs_axis, 1+rel_scale_unc[xs_index, 2, :], 1+rel_scale_unc[xs_index, 1, :], color=_order_color[order_item], alpha=0.6)
                ## denominator in ratio is lowest order in order_list
                ax2.fill_between(x_axis, (xs_all[xs_index]/xs_all[lowest_order])+rel_scale_unc[xs_index, 2, :],
                                (xs_all[xs_index]/xs_all[lowest_order])+rel_scale_unc[xs_index, 1, :], color=_order_color[order_item], alpha=0.46)


                #ax2.axis([xmin, xmax, 1+1.1*min(rel_scale_unc[xs_index, 2, :]), 1+1.1*max(rel_scale_unc[xs_index, 1, :])])
                ax2.set_ylabel('rel. scale_unc')

                ordernames += '_%s' %order_item

        ax2.axhline(y=1, xmin=0, xmax=xmax, color='k', linestyle='dotted')
        #now fixed axis, but should be flexible... (fixed it for testing)
        #ax2.axis([xmin, xmax, 0.3, 1.7]) #not final version
        fig.tight_layout()

        if given_filename is not None:
                filename = '%s.scaleunc-%s.%s' %(given_filename, vartype, ordernames[1:])
        else:
                filename = '%s.scaleunc-%s.%s' %(tablename, vartype, ordernames[1:])

        fig.savefig('%s.png' %filename)
        print '[fastnnlo_scaleunc]: saved as: %s.png' %filename



def main():
        #define arguments & options
        parser = argparse.ArgumentParser(epilog='')

        parser.add_argument('table', type=str, help='fastNLO table that shall be evaluated. To avoid mixing up the table argument with options taking multiple inputs like -o|--order, table must come first or must be separated from the -o option by another option!') #table is always required

        parser.add_argument('-p', '--pdfset', default='CT14nlo',
                                help='PDFset to evaluate fastNLO table.')
        parser.add_argument('-m', '--member', default=0, type=int,
                                help='Member of PDFset, default is 0.')
        parser.add_argument('-o', '--order', required=False, nargs='+', type=str, choices=['LO', 'NLO', 'NNLO'],
                            help='Blank-separated list of orders up to which the scale uncertainty is shown: LO, NLO, and/or NNLO. If nothing is chosen, successively produce plots for all orders that are available in table.')
        parser.add_argument('-s', '--scale', default=0, required=False, nargs='?', type=int,
                                help='For flex-scale tables define central scale choice for MuR and MuF by choosing 0 (kScale1) or 1 (kScale2).')
        parser.add_argument('-a', '--asymmetric', action="store_true",
                            help='If -a is chosen, use asymmetric (6P) scale variations; otherwise use symmetric ones (2P).')
        parser.add_argument('-f', '--filename', default=None, type=str,
                                help='Output filename (optional).')

        #parse arguments
        args = vars(parser.parse_args())

        #table name
        tablename = os.path.splitext(os.path.basename(args['table']))[0]
        tablename = os.path.splitext(tablename)[0] #to get rid of extension (.tab.gz or .tab)
        print '\n'
        print '[fastnnlo_scaleunc]: Analysing table', tablename

        #pdfset name
        pdfset = os.path.basename(args['pdfset'])
        print '[fastnnlo_scaleunc]: Using PDF set', pdfset

        #orders to be shown
        iorders = []
        iordmax = _text_to_order['NNLO']
        iordmin = 0
        if args['order'] is None:
                print '[fastnnlo_scaleunc]: Evaluate table up to highest available order.'
        else:
                for ord in args['order']:
                        iorders.append(_text_to_order[ord])
                iordmin = min(iorders)
                iordmax = max(iorders)
                print '[fastnnlo_scaleunc]: Evaluate table up to order(s)', args['order']

        #type of scale variation (symmetric vs asymmetric)
        if args['asymmetric']:
                scale_var_type = fastnlo.kAsymmetricSixPoint
                variation_type = 'scale uncertainty (6P)'
        else:
                scale_var_type = fastnlo.kSymmetricTwoPoint
                variation_type = 'scale uncertainty (2P)'

        #chosen scale
        scale_choice = args['scale']

        #given filename
        given_filename = args['filename']


        ###################### Start EVALUATION with fastNLO library #############################################
        fnlo = fastnlo.fastNLOLHAPDF(args['table'], args['pdfset'], args['member'])

        #Dictionary defining the fastNLO order settings
        #orders = {'LO':[True, False, False], 'NLO':[True, True, False], 'NNLO':[True, True, True]}

        #Get labeling for the x-axis
        #dimensionality of the table:
        ndim = fnlo.GetNumDiffBin()
        print '[fastnnlo_scaleunc]: Dimensions:', ndim

        #labels of all the dimensions:
        labels = fnlo.GetDimLabels()
        print '[fastnnlo_scaleunc]: Labels:', labels

        #label of first dimension:
        xlabel = fnlo.GetDimLabel(0)
        print '[fastnnlo_scaleunc]: x-label:', xlabel

        #creating x-axis
        bin_bounds = np.array(fnlo.GetObsBinsBounds(0))
        print '[fastnnlo_scaleunc]: bin_bounds.T: \n', bin_bounds.T, '\n'
        print '[fastnnlo_scaleunc]: bin_bounds.flatten()', bin_bounds.flatten(), '\n'

        x_axis = (bin_bounds.T[0]+bin_bounds.T[1])/2. #this is a list of bin centers
        xmin = 0.95*min(bin_bounds.ravel())
        xmax = 1.05*max(bin_bounds.ravel())
        print '[fastnnlo_scaleunc]: xmin =', xmin, ', xmax =', xmax

        #preparing x-errors (via bin_bounds) --> x_errors[0, :] are initially negative (positive via -1*), x_errors[1, :] positive
        x_errors = np.array([-1*(bin_bounds.T[0]-x_axis), bin_bounds.T[1]-x_axis])
        #print '[fastnnlo_scaleunc]: \n x_errors: ', x_errors, '\n'



        ##### Check: Which orders shall be investigated? Does table contain NNLO? ######
        #check existence of nnlo entry in table: ##
        o_existence = [False, False, False]
        cont_order = -1
        max_order  = 0
        for i in [0, 1, 2]:
                o_existence[i] = fnlo.SetContributionON(fastnlo.kFixedOrder, i, True)
                if o_existence[i]:
                        max_order   = i
                        if cont_order == i-1:
                                cont_order += 1
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

        #Set scale to user choice (default is 0)
        if (scale_choice in [0, 1]):
                fnlo.SetMuRFunctionalForm(scale_choice)
                fnlo.SetMuFFunctionalForm(scale_choice)
        else:
                sys.exit("Invalid choice of scale. Aborted!")
        #Now evaluate fastNLO table having a look at scale uncertainties
        xs_list = [] #will contain total cross section for LO, NLO, NNLO (note: could also be handled via GetScaleUnertaintyVec()[0] )
        rel_unc_list = [] #list for relative scale uncertainties (low, high) for LO, NLO, NNLO
        for n in order_list:
                for j in range(0, max_order+1):
                        if j <= _text_to_order[n]:
                                fnlo.SetContributionON(fastnlo.kFixedOrder, j, True)
                        else:
                                fnlo.SetContributionON(fastnlo.kFixedOrder, j, False)
                print '[fastnnlo_scaleunc]: \n'
                print '[fastnnlo_scaleunc]: Calculate XS for order: %s' %n, '\n'
                print '[fastnnlo_scaleunc]: ----  ----  ----  ----  ----  ----  ----  ----'
                fnlo.CalcCrossSection()
                xs_list.append(fnlo.GetCrossSection())

                ### Get scale uncertainties ###
                print '[fastnnlo_scaleunc]: Used scale factor MuF: ', fnlo.GetScaleFactorMuF()
                print '[fastnnlo_scaleunc]: Used scale factor MuR: ', fnlo.GetScaleFactorMuR(), '\n'
                print '[fastnnlo_scaleunc]: Calculate scale uncertainties \n'
                ## RELATIVE scale uncertainty with chosen type of scale variation (symmetric or asymmetric)
                rel_scale_unc_item = np.array(fnlo.GetScaleUncertaintyVec(scale_var_type)) #calculate this already for all accessible orders in any case
                rel_unc_list.append(rel_scale_unc_item)
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

        print '[fastnnlo_scaleunc]: Cross section xs_all uses ', order_list

        print xs_all, '\n \n'



        ## ABSOLUTE scale uncertainty
        num_orders = np.size(xs_all, 0) #length of axis 0 in xs_all equals number of orders
        abs_scale_unc = np.empty([num_orders, 2, len(x_axis)])
        for k in range(0, len(xs_all)):
                abs_scale_unc[k, 0, :] = rel_scale_unc[k, 2, :]*xs_all[k] #absolute uncertainties downwards (low)
                abs_scale_unc[k, 1, :] = rel_scale_unc[k, 1, :]*xs_all[k] #absolute uncertainties upwards (high)

        print '[fastnnlo_scaleunc]: Absolute Scale uncertainties downwards, upwards (order by order): \n'
        print abs_scale_unc, '\n'





        ############################## PLOTTING Scale Uncertainties ####################################################

        ### Plotting procedure if specific order has been chosen via -o ###
        ## Producing only one plot
        if (args['order'] is not None):
                if (len(order_list)==1):
                        order_index = _text_to_order[order_list[0]]
                        plotting_single(order_index, x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_choice)
                else: #plotting of combination plot which contains all requested orders (contains either (LO, NLO) or (LO, NNLO) or (NLO, NNLO) or (LO, NLO, NNLO))
                        # lowest_order will be used as denominator in ratio plot (rel_scale_unc)
                        plotting_multiple(iordmin, x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_choice)


        ### Plotting procedure to produce plots for all orders that are available ###
        elif (args['order'] is None):
                #for l in range(0, len(order_list)):
                for l in order_list:
                        order_index = _text_to_order[l]
                        plotting_single(order_index, x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_choice)
                        ###add also plot where all orders are plotted within the same figure

                lowest_order = 0 #because here all orders are considered --> 0 (=LO) is always the lowest.
                plotting_multiple(lowest_order, x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_choice)



if __name__ == '__main__':
        main()
