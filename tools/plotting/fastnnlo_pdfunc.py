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
import argparse, glob, os, re, sys, timeit
# Use matplotlib with Cairo offline backend for eps, pdf, png, or svg output
import matplotlib as mpl
#mpl.use('Agg')
mpl.use('Cairo')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import (FormatStrFormatter, ScalarFormatter, AutoMinorLocator, MultipleLocator)
from matplotlib import cm
# numpy
import numpy as np
# fastNLO for direct evaluation of interpolation grid
import fastnlo
from fastnlo import fastNLOLHAPDF
from fastnlo import SetGlobalVerbosity

# Redefine ScalarFormatter
class ScalarFormatterForceFormat(ScalarFormatter):
    def _set_format(self,vmin,vmax):  # Override function that finds format to use.
        self.format = "%1.2f"  # Give format here

# Action class to allow comma-separated list in options
class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

# Some global definitions
_formats        = {'eps':0, 'pdf':1, 'png':2, 'svg':3}
_text_to_order  = {'LO':0, 'NLO':1, 'NNLO':2}
_order_to_text  = {0:'LO', 1:'NLO', 2:'NNLO'}
_order_color    = {'LO':'g', 'NLO':'b', 'NNLO':'r'}
_colors         = ['tab:orange', 'tab:green', 'tab:purple', 'tab:blue', 'tab:brown']
_symbols        = ['s', 'X', 'o', '^', 'v']
_hatches        = ['-', '//', '\\', '|', '.']
_scale_to_text  = {0:'kScale1', 1:'kScale2', 2:'kQuadraticSum', 3:'kQuadraticMean', 4:'kQuadraticSumOver4',
                   5:'kLinearMean', 6:'kLinearSum', 7:'kScaleMax', 8:'kScaleMin', 9:'kProd',
                   10:'kS2plusS1half', 11: 'kPow4Sum', 12:'kWgtAvg', 13:'kS2plusS1fourth', 14:'kExpProd2', 15:'kExtern'}
_pdfbasenames   = ['ABM11', 'ABMP15', 'ABMP16', 'CJ12', 'CJ15', 'CT10', 'CT14', 'HERAPDF20', 'JR14',
                   'MMHT2014', 'MSTW2008', 'NNPDF23', 'NNPDF30', 'NNPDF31']
_debug          = False

#####################################################################################


# Plotting function to be elaborated
def plotting(x_axis, xmin, xmax, xs_chosen, rel_pdf_unc, abs_pdf_unc, xlabel, tablename, order_list, given_filename, scale_name, pdfsets, formats):

        pdfnicenames = []
        for pdf in pdfsets:
            nicename = 'Undefined'
            for pdfn in _pdfbasenames:
                if pdfn in pdf: nicename = pdfn
            pdfnicenames.append(nicename)

        gs = gridspec.GridSpec(3,3)
        fig = plt.figure(figsize=(7,7))
        ax1 = plt.subplot(gs[:-1,:])

        # Check whether only one PDF set is being looked at or several PDF sets.
        if (len(pdfsets)>1): #one plot contains all chosen pdf sets for one specific order + summary plot if 3 orders in order_list
                # Several PDFs will be plotted 'next to each other' --> need tiny shift to be displayed
                if len(pdfsets)==2:
                        shift_list=[0.98, 1.02]
                elif len(pdfsets)==3:
                        shift_list=[0.96, 1.00, 1.04]
                elif len(pdfsets)==4:
                        shift_list=[0.94, 0.98, 1.02, 1.06]
                elif len(pdfsets)==5:
                        shift_list=[0.92, 0.96, 1.00, 1.04, 1.08]
                else:
                        print '[fastnnlo_pdfunc]: Too many PDFs to produce combined plot. Aborted!'
                        print '[fastnnlo_pdfunc]: Max. possible: 5. Number of given PDF Sets: %s'%len(pdfsets)
                        sys.exit()

                order_index = 0 #necessary, because even though LO might not even be investigated, we need index 0, which could then correspond to NLO or NNLO
                ordernames = ''
                pdffilenames = ''
                for pdf in pdfsets:
                        pdffilenames += '_%s' %pdf

                for order_item in order_list: #producing one plot per order (comparing all the given pdfs)
                        fig = plt.figure(figsize=(7,7))
                        ax1 = plt.subplot(gs[:-1,:])

                        print '[fastnnlo_pdfunc]: Producing %s plot.'%order_item
#                        color = iter(cm.rainbow(np.linspace(0,1,len(pdfsets))))
                        patches = [] #later needed for legend

                        pdf_index = 0
#                        for pdf, shift in zip(pdfsets, shift_list): #go through all the pdfs
                        for pdf, shift in zip(pdfnicenames, shift_list): #go through all the pdfs
#                                c = next(color)
                                c = _colors[pdf_index]
                                ax1.errorbar(x_axis*shift, xs_chosen[pdf_index, order_index, :], yerr=abs(abs_pdf_unc[pdf_index, order_index, :, :]), elinewidth=1, linewidth=0.0, ms=6, marker=_symbols[pdf_index], color=c, fmt='.', label=pdf)
                                pdf_index += 1

                        ax1.set_xlim([xmin, xmax])
                        ax1.set_xscale('log', nonposx='clip')
                        ax1.set_yscale('log', nonposy='clip')
                        ax1.set_xlabel(r'%s' %xlabel, horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
                        ax1.set_ylabel(r'$\sigma \pm \Delta\sigma(\mathrm{PDF})$', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, rotation=90, labelpad=16)
                        ax1.legend(fontsize=10, numpoints=1)
                        ax1.text(0.03, 0.15, 'Reference PDF: %s' %pdfnicenames[0], horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)
                        ax1.text(0.03, 0.10, 'Scale: %s' %scale_name, horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)
                        ax1.text(0.03, 0.05, 'Order: %s' %order_item, horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)

                        ax1.set_title('%s' %tablename)



                        ax2 = plt.subplot(gs[2, :])
                        ax2.set_xlim([xmin, xmax])

                        #ratioplot normalised to first given pdf.
#                        color = iter(cm.rainbow(np.linspace(0,1,len(pdfsets))))
                        for p in range(0, len(pdfsets)): #for pdf in pdfsets
#                                c = next(color)
                                c = _colors[p]
                                # divide xs for each PDF by first given PDF (xs)
                                ax2.semilogx(x_axis, xs_chosen[p, order_index, :]/xs_chosen[0, order_index, :], '.',
                                             ms=6, marker=_symbols[p], ls='dashed', linewidth=0.2, color=c, label=pdfnicenames[p])
                                ax2.fill_between(x_axis, (xs_chosen[p, order_index, :]/xs_chosen[0, order_index, :])+rel_pdf_unc[p, order_index, 2, :],
                                                 (xs_chosen[p, order_index, :]/xs_chosen[0, order_index, :])+rel_pdf_unc[p, order_index, 1, :], color=c,
                                                 alpha=0.3, hatch=_hatches[p])
                                #patch for current pdf (needed for labeling/legend)
                                patches.append(mpl.patches.Rectangle((0, 0), 0, 0, color=c, label=pdfnicenames[p], alpha=0.4))
                                ax2.add_patch(patches[p])

                        ax2.set_ylabel(r'Ratio to ref. PDF')
                        ax2.axhline(y=1, xmin=0, xmax=1, linewidth=0.4, color='k', linestyle='dotted')

                        fig.tight_layout()

                        order_index += 1
                        ordernames += '_%s' %order_item

                        if given_filename is not None:
                                filename = '%s.pdfunc-%s.%s' %(given_filename, order_item, pdffilenames[1:])
                        else:
                                filename = '%s.pdfunc-%s.%s' %(tablename, order_item, pdffilenames[1:])

                        for fmt in formats:
                            figname = '%s.%s' %(filename, fmt)
                            fig.savefig(figname)
                            print '[fastnnlo_pdfunc]: Plot saved as:', figname
                        plt.close(fig)

                # Summary plot if 3 orders have been looked at
                if len(order_list)==3: #could adjust this to arbitrary length of order
                        gs2 = gridspec.GridSpec(2, 2)
                        gs2.update(wspace=0.4, hspace=0.4) #for more whitespace between the subplots
                        ##fig2 = plt.figure(figsize=(7,7))
                        fig2 = plt.figure(figsize=(12,10)) #different size compared to other plots that are produced
                        ax0 = plt.subplot(gs2[0,0])
                        ax1 = plt.subplot(gs2[0,1])
                        ax2 = plt.subplot(gs2[1,0])
                        ax3 = plt.subplot(gs2[1,1])

                        # subplots containing ratioplot for each order (comparing all pdfs, ratio to first given pdf)
                        # and one subplot for first pdf (comparing all orders, ratio to 1st given order)
                        ##for ax, order in zip([ax0, ax1, ax2], ['LO', 'NLO', 'NNLO']): #here ax0 will always contain LO, ax1 NLO and ax2 NNLO
                        order_index = 0 #to go through order_list item by item, keeping in mind that order_list is not sorted
                        patches = [] #later needed for legend
                        for ax, order in zip([ax0, ax1, ax2], order_list): #here it can happen, that NNLO is plotted in ax0, depending on content of order_list
                                patches=[] #currently necessary, not needed if one only makes one legend for all axes
#                                color = iter(cm.rainbow(np.linspace(0,1,len(pdfsets))))
                                for p in range(0, len(pdfsets)): #for pdf in pdfsets
#                                        c = next(color)
                                        c = _colors[pdf_index]
                                        ax.semilogx(x_axis, xs_chosen[p, order_index, :]/xs_chosen[0, order_index, :], '.', ms=6, marker=_symbols[p] ,ls='dashed', linewidth=0.2, color=c)#, label=pdfsets[p])
                                        ax.fill_between(x_axis, (xs_chosen[p, order_index, :]/xs_chosen[0, order_index, :])+rel_pdf_unc[p, order_index, 2, :],
                                                        (xs_chosen[p, order_index, :]/xs_chosen[0, order_index, :])+rel_pdf_unc[p, order_index, 1, :], color=c,
                                                        alpha=0.3)
                                        #Alternative: create the patches only once (for ax0) --> use the same for ax1 and ax2
                                        ###if ax==ax0:
                                        ###     patches.append(mpl.patches.Rectangle((0, 0), 0, 0, color=c, label=pdfsets[p], alpha=0.3))
                                        ###     ax.add_patch(patches[p])
                                        patches.append(mpl.patches.Rectangle((0, 0), 0, 0, color=c, label=pdfnicenames[p], alpha=0.3))
                                        ax.add_patch(patches[p])

                                ax.set_xlabel('%s'%xlabel)
                                ax.set_ylabel(r'$\sigma \pm \Delta\sigma(\text{PDF})$')
                                ax.text(0.03, 0.10, 'Reference PDF: %s'%pdfsets[0], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
                                ax.text(0.03, 0.05, 'Order: %s'%_order_to_text[order_index], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)

                                ax.axhline(y=1, xmin=0, xmax=1, linewidth=0.4, color='k', linestyle='dotted')

                                order_index += 1
                                ax.legend()
                        #plt.legend(handles=patches, mode='expand', bbox_to_anchor=(1.02, 1), loc='upper left') #'upper left' refers to the legend box ######

                        order_index = 0
                        for order_item in order_list:
                                ax3.semilogx(x_axis, xs_chosen[0, order_index, :]/xs_chosen[0, 0, :], '.', ms=6, ls='dashed', linewidth=0.2,
                                                color=_order_color[order_item], label=order_item)

                                ax3.fill_between(x_axis, (xs_chosen[0, order_index, :]/xs_chosen[0, 0, :])+rel_pdf_unc[0, order_index, 2, :],
                                                (xs_chosen[0, order_index, :]/xs_chosen[0, 0, :])+rel_pdf_unc[0, order_index, 1, :], color=_order_color[order_item],
                                                alpha=0.46)#, hatch=_order_hatch[order_item])

                                order_index += 1
                        ax3.set_title('%s'%pdfsets[0])
                        ax3.set_xlabel('%s'%xlabel)
                        ax3.set_ylabel('rel. pdf_unc')
                        ax3.set_ylim([0.82, None]) #to avoid problems with text
                        ax3.axhline(y=1, xmin=0, xmax=1, linewidth=0.4, color='k', linestyle='dotted')
                        ax3.text(0.03, 0.10, 'Reference order: %s'%order_list[0], horizontalalignment='left', verticalalignment='bottom', transform=ax3.transAxes)
                        ax3.text(0.02, 0.02, 'PDF Set: %s'%pdfnicenames[0], horizontalalignment='left', verticalalignment='bottom', transform=ax3.transAxes)
                        ax3.legend()

                        fig2.suptitle('%s'%tablename, fontsize=16)


                        ####gs2.tight_layout(fig2, rect=[0, 0.3, 1, 0.95]) #just testing
                        ####gs2.tight_layout(fig2) #just testing

                        if given_filename is not None:
                                filename = '%s.pdfunc-summary.%s' %(given_filename, pdffilenames[1:])
                        else:
                                filename = '%s.pdfunc-summary.%s' %(tablename, pdffilenames[1:])
                        for fmt in formats:
                            figname = '%s.%s' %(filename, fmt)
                            fig.savefig(figname)
                            print '[fastnnlo_pdfunc]: Plot saved as:', figname
                        plt.close(fig2)

                else: print '[fastnnlo_pdfunc]: For less than 3 orders --> no summary plot produced.'



        elif (len(pdfsets)==1): #plotting with given or default pdf, one plot contains all given/default orders
                # For plotting several orders 'next to each other', handling via shift from bincenter
                if len(order_list)==1:
                        shift_list=[1.00]
                elif len(order_list)==2:
                        shift_list=[0.98, 1.02]
                elif len(order_list)==3:
                        shift_list=[0.98, 1.00, 1.02]
                else:
                        print '[fastnnlo_pdfunc]: Too many orders to plot. Aborted!'
                        print '[fastnnlo_pdfunc]: Requested orders: %s'%order_list
                        sys.exit()

                order_index = 0 #this is necessary because the orders in order_list are unsorted (order_index is not tied to specific order)
                for order_item, shift in zip(order_list, shift_list):
                        ax1.errorbar(x_axis*shift, xs_chosen[0, order_index, :], yerr=abs(abs_pdf_unc[0, order_index, :, :]), elinewidth=1,
                                        linewidth=0.0, ms=6, color=_order_color[order_item], fmt='.', label=order_item)
                        order_index += 1

                ax1.set_xlim([xmin, xmax])
                ax1.set_xscale('log', nonposx='clip')
                ax1.set_yscale('log', nonposy='clip')
                ax1.set_xlabel('%s' %xlabel)
                ax1.set_ylabel('XS with abs_pdf_unc', rotation=90)
                ax1.legend(fontsize=10, numpoints=1)
                ax1.text(0.03, 0.05, 'PDF set: %s' %pdfnicenames[0], horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)

                ax1.set_title('%s' %tablename)

                # Ratio subplot with relative pdf uncertainties, denominator in ratio = first order in order_list (not necessarly LO!)
                ax2 = plt.subplot(gs[2, :])
                ax2.set_xlim([xmin, xmax])

                ordernames = ''
                order_index = 0
                for order_item in order_list:
                        ax2.semilogx(x_axis, xs_chosen[0, order_index, :]/xs_chosen[0, 0, :], '.', ms=6, ls='dashed', linewidth=0.2,
                                        color=_order_color[order_item], label=order_item)
                        # divide xs fo each order by xs of first given order! (--> ratio)
                        ax2.fill_between(x_axis, (xs_chosen[0, order_index, :]/xs_chosen[0, 0, :])+rel_pdf_unc[0, order_index, 2, :],
                                        (xs_chosen[0, order_index, :]/xs_chosen[0, 0, :])+rel_pdf_unc[0, order_index, 1, :], color=_order_color[order_item],
                                        alpha=0.46)#, hatch=_order_hatch[order_item])

                        order_index += 1
                        ordernames += '_%s' %order_item

                ax2.set_ylabel('rel. pdf_unc')
                ax2.text(0.02, 0.86, 'Reference order: %s'%order_list[0], horizontalalignment='left', verticalalignment='bottom', transform=ax2.transAxes)
                ax2.axhline(y=1, xmin=0, xmax=1, linewidth=0.4, color='k', linestyle='dotted')

                fig.tight_layout()

                if given_filename is not None:
                        filename = '%s.pdfunc-%s.%s' %(given_filename, pdfsets[0], ordernames[1:])
                else:
                        filename = '%s.pdfunc-%s.%s' %(tablename, pdfsets[0], ordernames[1:])

                fig.savefig('%s.png' %filename)
                print '[fastnnlo_pdfunc]: Saved as: %s.png' %filename


#####################################################################################

def main():
        start_time = timeit.default_timer() #just for measuring wall clock time - not necessary
        # Define arguments & options
        parser = argparse.ArgumentParser(epilog='',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        # Positional arguments
        parser.add_argument('table', type=str, nargs='+',
                            help='Filename glob of fastNLO tables to be evaluated. This must be specified!')
        # Optional arguments
        parser.add_argument('-f', '--filename', default=None, type=str,
                                help='Output filename (optional).')
        parser.add_argument('--format', required=False, nargs='?', type=str, action=SplitArgs,
                            help='Comma-separated list of plot formats to use: eps, pdf, png, svg. If nothing is chosen, png is used.')
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
        parser.add_argument('-v', '--verbose', action='store_true',
                            help='Increase output verbosity.')

        # Parse arguments
        args = vars(parser.parse_args())

        # List of table names
        files = args['table']
        print '\n'
        print '[fastnnlo_pdfunc]: Analysing table list: '
        for file in files:
                print '[fastnnlo_pdfunc]:   ', file

        # PDF set name
        pdfsets = []
        if args['pdfset'] is None:
                pdfsets = ['CT14nnlo']
                print '[fastnnlo_pdfunc]: Using default PDF set: ', pdfsets[0]
        else:
                for pdf in args['pdfset']:
                        pdfsets.append(pdf)
                print '[fastnnlo_pdfunc]: Using PDF set(s): ', pdfsets

        # Orders to be shown
        iorders = []
        iordmin = _text_to_order['LO']
        iordmax = _text_to_order['NNLO']
        if args['order'] is None:
                print '[fastnnlo_pdfunc]: Evaluate table up to highest available order.'
        else:
                for ord in args['order']:
                        if _text_to_order.has_key(ord):
                                iorders.append(_text_to_order[ord])
                        else:
                                print '[fastnnlo_pdfunc]: Illegal order specified, aborted!'
                                print '[fastnnlo_pdfunc]: Order list:', args['order']
                                exit(1)
                iordmin = min(iorders)
                iordmax = max(iorders)
                print '[fastnnlo_pdfunc]: Evaluate table up to order(s)', args['order']

        # Given filename
        given_filename = args['filename']

        # Plot formats to use
        formats = args['format']
        if formats is None: formats = ['png']
        for fmt in formats:
            if not _formats.has_key(fmt):
                print '[fastnnlo_pdfunc]: Illegal format specified, aborted!'
                print '[fastnnlo_pdfunc]: Format list:', args['format']
                exit(1)

        # Scale choice
        scale_choice = args['scale']

        # Verbosity
        verb = args['verbose']

        # Loop over table list
        for table in files:
                print '[fastnnlo_pdfunc]: Analysing table: ', table
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
                        print '\n'
                        print '[fastnnlo_pdfunc]: Table Dimensions: ', ndim

                # Labels of all the dimensions:
                dimlabels = fnlo.GetDimLabels()
                if verb:
                        print '[fastnnlo_pdfunc]: Dimension Labels: ', dimlabels

                # Label of first dimension:
                xlabel = fnlo.GetDimLabel(0)
                if verb:
                        print '[fastnnlo_pdfunc]: x-label: ', xlabel

                # Creating x-axis
                bin_bounds = np.array(fnlo.GetObsBinsBounds(0))
                if verb:
                        print '[fastnnlo_pdfunc]: bin_bounds.T: \n', bin_bounds.T, '\n'
                        print '[fastnnlo_pdfunc]: bin_bounds.flatten(): \n', bin_bounds.flatten(), '\n'

                x_axis = (bin_bounds.T[0]+bin_bounds.T[1])/2. #this is a list of bin centers
                xmin = 0.95*min(bin_bounds.ravel())
                xmax = 1.05*max(bin_bounds.ravel())
                if verb:
                        print '[fastnnlo_pdfunc]: xmin=%s, xmax=%s. \n' %(xmin, xmax)

                # Preparing x-errors (via bin_bounds) --> x_errors[0, :] are initially negative (positive via -1*), x_errors[1, :] positive
                x_errors = np.array([-1*(bin_bounds.T[0]-x_axis), bin_bounds.T[1]-x_axis])
                if verb:
                        print '[fastnnlo_pdfunc]: x_errors: \n', x_errors, '\n'

                # Check existence of orders in table
                lflex = fnlo.GetIsFlexibleScaleTable()
                scale_name = 'scale1'
                o_existence = [False, False, False]
                cnt_order = -1
                max_order  = 0
                for i in range(3):
                        o_existence[i] = fnlo.SetContributionON(fastnlo.kFixedOrder, i, True)
                        if o_existence[i]:
                                max_order = i
                                if not lflex:
                                        if scale_choice != 0:
                                                print '[fastnnlo_pdfunc]: Invalid choice of scale = ', scale_choice, ' Aborted!'
                                                print '[fastnnlo_pdfunc]: For fixed-scale tables only the default=0 is allowed.'
                                                exit(1)
                                        else:
                                                scale_name = fnlo.GetScaleDescription(i,0)
                                else:
                                        if scale_choice < 2:
                                                scale_name = fnlo.GetScaleDescription(i,scale_choice)
                                        else:
                                                scl0 = fnlo.GetScaleDescription(i,0)
                                                scl1 = fnlo.GetScaleDescription(i,1)
                                                scale_name = _scale_to_text[scale_choice]+'_'+scl0+'_'+scl1
                                if cnt_order == i-1:
                                        cnt_order += 1
                if verb:
                        print '[fastnnlo_pdfunc]: Table has continuous orders up to', cnt_order, 'and a maximal order of', max_order

                if iordmax > cnt_order:
                        print '[fastnnlo_pdfunc]: Invalid choice of orders. Aborted!'
                        print '[fastnnlo_pdfunc]: Highest order requested is', _order_to_text[iordmax], 'but orders are available only up to', cnt_order
                        exit(1)

                order_list = []
                if args['order'] is None:
                        for iord in range(cnt_order+1):
                                order_list.append(_order_to_text[iord])
                else:
                        order_list = args['order']

                print '[fastnnlo_pdfunc]: List of requested orders:', order_list

                # For flexible-scale tables set scale to user choice (default is 0)

                if lflex:
                        print '[fastnnlo_pdfunc]: Setting requested scale choice for flexible-scale table:', scale_choice
                        fnlo.SetMuRFunctionalForm(scale_choice)
                        fnlo.SetMuFFunctionalForm(scale_choice)
                else:
                        if scale_choice == 0:
                                print '[fastnnlo_pdfunc]: Evaluating fixed-scale table. Scale choice must be', scale_choice
                        else:
                                print '[fastnnlo_pdfunc]: No scale choice possible for fixed-scale table. Aborted!'
                                print '[fastnnlo_pdfunc]: scale_choice = ', scale_choice
                                exit(1)

                # Now evaluate fastNLO table for each of the given PDF sets
                xs_list = [] #will contain for each given pdf the total cross section for LO, NLO, NNLO (or requested orders)
                rel_unc_list = [] #list for relative pdf uncertainties (low, high) for each PDF in LO, NLO, NNLO

                # Outermost loop: go through all the given pdf sets
                # Evaluate fastNLO table focusing on PDF uncertainties
                for pdf in pdfsets:
                        xs_list_tmp = [] #temporary xs list for single pdf
                        rel_unc_list_tmp = [] #only temporary list. will be renewed for each pdf
                        print '-----------------------------------------------------------------------------------------------'
                        print '#############################  %s  ########################################'%pdf
                        print '-----------------------------------------------------------------------------------------------'
                        print '[fastnnlo_pdfunc]: Calculate XS and uncertainty for %s \n'%pdf
                        fnlo = fastnlo.fastNLOLHAPDF(table, pdf, args['member'])


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
                                        print '[fastnnlo_pdfunc]: Calculate XS for order: %s' %n
                                        print '[fastnnlo_pdfunc]: ---- ---- ---- ---- ---- ---- ---- \n'

                                fnlo.CalcCrossSection()
                                xs_list_tmp.append(fnlo.GetCrossSection())

                                ### Get PDF uncertainties ###
                                ## RELATIVE pdf uncertainty
                                rel_pdf_unc_item = np.array(fnlo.GetPDFUncertaintyVec(fastnlo.kLHAPDF6)) #now calculated for order n
                                ####rel_pdf_unc_item = np.ones((3,33))                  #####this is just for testing the plotting function without big calculation
                                ####rel_pdf_unc_item[2,:]=(-1)*rel_pdf_unc_item[2,:]    #####only for testing
                                rel_unc_list_tmp.append(rel_pdf_unc_item)
                                if verb:
                                        print '[fastnnlo_pdfunc]: Cross section in %s: \n'%n, np.array(xs_list_tmp)[-1], '\n' #print most recent xs
                                        print '[fastnnlo_pdfunc]: Relative PDF uncertainty in %s: \n'%n
                                        print rel_pdf_unc_item, '\n' #has 3 entries: central value (xs), unc_high, unc_low
                                        print '-----------------------------------------------------------------------------------------------'
                                        print '-----------------------------------------------------------------------------------------------'
                        xs_list.append(xs_list_tmp)
                        rel_unc_list.append(rel_unc_list_tmp)
                xs_chosen = np.array(xs_list)
                rel_pdf_unc = np.array(rel_unc_list) #Remember: both arrays here do only contain the CHOSEN orders! (or per default LO, NLO, (NNLO))
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
                        print '[fastnnlo_pdfunc]: For each pdf %s xs_chosen contains XS corresponding to %s.'%(pdfsets, order_list)
                        print xs_chosen, '\n \n'
                        print 'Size xs_chosen: ', np.shape(xs_chosen) ###test
                        print 'Size rel_pdf_unc: ', np.shape(rel_pdf_unc) ###test

                ## ABSOLUTE pdf uncertainty
                num_pdfsets = np.size(xs_chosen, 0) #length of axis 0 in xs_chosen equals number of (investigated) pdfsets
                num_orders = np.size(xs_chosen, 1) #length of axis 1 in xs_chosen equals number of (investigated) orders
                abs_pdf_unc = np.empty([num_pdfsets, num_orders, 2, len(x_axis)])
                for p in range(0, num_pdfsets):
                        for k in range(0, num_orders): # k is order index
                                abs_pdf_unc[p, k, 0, :] = rel_pdf_unc[p, k, 2, :]*xs_chosen[p, k, :] #absolute uncertainties downwards (low)
                                abs_pdf_unc[p, k, 1, :] = rel_pdf_unc[p, k, 1, :]*xs_chosen[p, k, :] #absolute uncertainties upwards (high)

                if verb:
                        print '[fastnnlo_pdfunc]: Absolute PDF uncertainties downwards, upwards for %s in %s: \n'%(pdfsets, order_list)
                        print abs_pdf_unc, '\n'




                                ################################## Start PLOTTING of the pdf uncertainties ######################################################################

                plotting(x_axis, xmin, xmax, xs_chosen, rel_pdf_unc, abs_pdf_unc, xlabel, tablename, order_list, given_filename, scale_name, pdfsets, formats)


        stop_time = timeit.default_timer()
        timediff = stop_time-start_time
        print 'Time: %s sec = %s min'%(timediff, round(timediff/60., 2))


if __name__ == '__main__':
        main()
