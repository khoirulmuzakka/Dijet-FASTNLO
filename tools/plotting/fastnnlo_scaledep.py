#!/usr/bin/env python2
#-*- coding:utf-8 -*-

##############################################
#
# Plot the scale dependence for each
# observable bin
#
# Created by K. Rabbertz, 16.02.2019
#
#############################################
#
import argparse, glob, os, re, sys
# Use matplotlib with Cairo offline backend for eps, pdf, png, or svg output
import matplotlib as mpl
#mpl.use('Agg')
mpl.use('Cairo')
import matplotlib.gridspec as gridspec
import matplotlib.lines as mpllines
import matplotlib.patches as mplpatches
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

#params = {'xtick.labelsize':'large',
#          'xtick.major.size': 5,
#          'xtick.major.width': 2,
#          'xtick.minor.size': 3.5,
#          'xtick.minor.width': 1,
#          'ytick.labelsize':'large',
#          'ytick.major.size': 5,
#          'mathtext.fontset':'cm'}
#pylab.rcParams.update(params)

#print 'a',mpl.get_configdir()
#exit(1)

# Redefine ScalarFormatter
class ScalarFormatterForceFormat(ScalarFormatter):
    def _set_format(self,vmin,vmax):  # Override function that finds format to use.
        self.format = "%1.2f"  # Give format here

# Action class to allow comma-separated list in options
class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

# Some global definitions for orders to show
_formats        = {'eps':0, 'pdf':1, 'png':2, 'svg':3}
_text_to_order  = {'LO':0, 'NLO':1, 'NNLO':2}
_order_to_text  = {0:'LO', 1:'NLO', 2:'NNLO'}
_order_to_color = {'LO':'g', 'NLO':'b', 'NNLO':'r'}
_scale_to_text  = {0:'kScale1', 1:'kScale2', 2:'kQuadraticSum', 3:'kQuadraticMean', 4:'kQuadraticSumOver4',
                   5:'kLinearMean', 6:'kLinearSum', 7:'kScaleMax', 8:'kScaleMin', 9:'kProd',
                   10:'kS2plusS1half', 11: 'kPow4Sum', 12:'kWgtAvg', 13:'kS2plusS1fourth', 14:'kExpProd2', 15:'kExtern'}
_debug          = False

#####################################################################################


# Function for plotting the mur scale dependence of a list of orders into one figure
# Optionally, ratios are shown as subplot with respect to the first order given in the order list.
def plotting(x_axis, xmin, xmax, iobs, xs_cn, xs_fl, xs_fu, dxsr_cn, xind, tablename, order_list, filename, scale_name, pdfset, labels, ylabel, borders, ratio, formats, verb):
        if ratio:
                gs = gridspec.GridSpec(3,3)
                fig = plt.figure(figsize=(7,7))
                ax1 = plt.subplot(gs[:-1,:])
        else:
                gs = gridspec.GridSpec(2,3)
                fig = plt.figure(figsize=(7,5))
                ax1 = plt.subplot(gs[:,:])
        ax1.set_xlim([xmin, xmax])
        ax1.set_xscale('log', nonposx='clip')
        ax1.set_yscale('linear')
        # Texify labels
#        tlabels = []
#        for lb in labels:
#            lb = lb.replace('_', '\_')
#            lb = lb.replace('[', '\[')
#            lb = lb.replace(']', '\]')
#            tlabels.append(lb)
#        if verb: print '[fastnnlo_scaledep]: Texified labels:', tlabels
        # Label of last dimension
        last = labels[-1].split('_',1)[0]
#        xlabel = '$\mu_R\,/\,%s$' %last
        xlabel = '$\mu_R\,/\,\mu_0$'
        if verb:
                print '[fastnnlo_scaledep]: x-label:', xlabel
        ax1.set_xlabel(r'%s' %xlabel, horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
        myxticks = [0.2, 0.5, 1.0, 2.0, 4.0, 8.0]
        ax1.set_xticks(myxticks)
        ax1.set_xticklabels(myxticks)
        ax1.set_ylabel(ylabel, horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, rotation=90, labelpad=16)
        ax1.grid(None)
        ax1.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        ax1.tick_params(direction='in', which='both')
        yfmt = ScalarFormatterForceFormat()
        yfmt.set_powerlimits((0,0))
        ax1.yaxis.set_major_formatter(yfmt)
        plt.setp(ax1.get_xticklabels(), fontsize=12)
        ax1.text(0.03, 0.10, 'PDF set: %s' %pdfset, horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)
        ax1.text(0.03, 0.05, 'Scale: %s' %scale_name, horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)
        yt = 1.00
        for id in range(len(labels)):
            yt -= 0.05
            text = '%1.1f < %s < %1.1f' % (borders[iobs][id][0], labels[id], borders[iobs][id][1])
            ax1.text(0.5, yt, text, horizontalalignment='center', verticalalignment='top', transform=ax1.transAxes)
        ax1.set_title('%s' %tablename)

        # Plot all x section results in order_list with muf variation as uncertainty bands
        iorder = -1
        xs_min = np.minimum.reduce([xs_cn,xs_fl,xs_fu])
        xs_max = np.maximum.reduce([xs_cn,xs_fl,xs_fu])
        dxu    = xs_max - xs_cn
        dxl    = xs_cn - xs_min
        dxur   = np.divide(dxu, xs_cn, out=np.ones_like(dxu), where=xs_cn!=0)
        dxlr   = np.divide(dxl, xs_cn, out=np.ones_like(dxl), where=xs_cn!=0)

        handles = []
        labels  = []
        for order in order_list:
                iorder += 1
                if xind[0] > -1:
                    dxstfl = np.multiply(dxsr_cn,xs_fl[:,xind[0]])
                    ax1.errorbar(x_axis[xind[0]], xs_fl[iorder,xind[0]], yerr=dxstfl[iorder], elinewidth=1, linewidth=0.0, ms=12, color=_order_to_color[order], fmt='.', label='_')
                if xind[1] > -1:
                    dxstcn = np.multiply(dxsr_cn,xs_cn[:,xind[1]])

                    points = ax1.errorbar(x_axis[xind[1]], xs_cn[iorder,xind[1]], yerr=dxstcn[iorder], elinewidth=1, linewidth=0.0, ms=12, color=_order_to_color[order], fmt='.', label=order)
                    handles.append(points)
                    labels.append(order)
                if xind[2] > -1:
                    dxstfu = np.multiply(dxsr_cn,xs_fu[:,xind[2]])
                    ax1.errorbar(x_axis[xind[2]], xs_fu[iorder,xind[2]], yerr=dxstfu[iorder], elinewidth=1, linewidth=0.0, ms=12, color=_order_to_color[order], fmt='.', label='_')
                ax1.fill_between(x_axis, xs_min[iorder,:], xs_max[iorder,:], color=_order_to_color[order], alpha=0.5)
                linel, = ax1.plot(x_axis, xs_fl[iorder,:], '--', color=_order_to_color[order], label=r'$\mu_F\,/\,\mu_R=0.5$')
                linec, = ax1.plot(x_axis, xs_cn[iorder,:], '-', color=_order_to_color[order], label=r'$\mu_F\,/\,\mu_R=1.0$')
                lineu, = ax1.plot(x_axis, xs_fu[iorder,:], ':', color=_order_to_color[order], label=r'$\mu_F\,/\,\mu_R=2.0$')
                if iorder==len(order_list)-1:
                    handles.append(linel)
                    handles.append(linec)
                    handles.append(lineu)
                    labels.append(r'$\mu_F\,/\,\mu_R=0.5$')
                    labels.append(r'$\mu_F\,/\,\mu_R=1.0$')
                    labels.append(r'$\mu_F\,/\,\mu_R=2.0$')

        leg = ax1.legend(handles, labels, fontsize=10, numpoints=1)
        # Set line colors in legend to black
        ih = 0
        for h in leg.legendHandles:
            if (ih > 2): h.set_color('black')
            ih+=1


        # Do subplot with ratios and relative muf uncertainties; denominator in ratio = first order in order_list
        if ratio:
                ax2 = plt.subplot(gs[2, :])
                #        ax2 = plt.subplot(gs[2, :], sharex=ax1)
                ax2.set_xlim([xmin, xmax])
                ax2.set_xscale('log', nonposx='clip')
                ax2.set_yscale('linear')
                ax2.set_xlabel(r'%s' %xlabel, horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
                myxticks = [0.2, 0.5, 1.0, 2.0, 4.0, 8.0]
                ax2.set_xticks(myxticks)
                ax2.set_xticklabels(myxticks)
                ax2.set_ylabel(r'Ratios to %s' %order_list[0])

                # Plot all x section ratios with respect to first in order_list
                iorder = -1
                ordernames = ''
                for order in order_list:
                        iorder += 1
                        ordernames += '_%s' %order
                        ax2.semilogx(x_axis, xs_cn[iorder,:]/xs_cn[0,:], ls='-', lw=1.0, color=_order_to_color[order], label=order)
                        ax2.fill_between(x_axis, xs_min[iorder,:]/xs_cn[0,:], xs_max[iorder,:]/xs_cn[0,:], color=_order_to_color[order], alpha=0.50)

        fig.tight_layout()

        for fmt in formats:
            figname = '%s.%s' %(filename, fmt)
            fig.savefig(figname)
            print '[fastnnlo_scaledep]: Plot saved as:', figname
#        print(plt.get_fignums())
        plt.close(fig)


#####################################################################################

def main():
        # Define arguments & options
        parser = argparse.ArgumentParser(epilog='',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        # Positional arguments
        parser.add_argument('table', type=str, nargs='+',
                            help='Filename glob of fastNLO tables to be evaluated. This must be specified!')
        # Optional arguments
        parser.add_argument('-b', '--binpower', default=3, type=int,
                            choices=range(-5,6), metavar='[-5,5]',
                            help='Power b for binary log range 2^(-b) to 2^b of mur scale factor.')
        parser.add_argument('-d', '--datfiles', required=False, nargs='?', type=str, action=SplitArgs,
                            help='Comma-separated list of NNLOJET datfile basenames with statistical uncertainties for each order to show. If nothing is chosen, statistical uncertainties are ignored.')
        parser.add_argument('-f', '--filename', default=None, type=str,
                            help='Set desired basename for output filenames instead of tablename.')
        parser.add_argument('--format', required=False, nargs='?', type=str, action=SplitArgs,
                            help='Comma-separated list of plot formats to use: eps, pdf, png, svg. If nothing is chosen, png is used.')
        parser.add_argument('-l', '--logpoints', default=7, type=int,
                            choices=range(3,31), metavar='[3-30]',
                            help='Number of equidistant points in log_2(xmur) from 2^(-b) to 2^b.')
        parser.add_argument('-m', '--member', default=0, type=int,
                            help='Member of PDFset to use.')
        parser.add_argument('-o', '--order', required=False, nargs='?', type=str, action=SplitArgs,
                            help='Comma-separated list of orders to show: LO, NLO, and/or NNLO. If nothing is chosen, show all orders available in table.')
        parser.add_argument('-p', '--pdfset', default='CT14nlo',
                            help='PDFset to use with fastNLO table.')
        parser.add_argument('-r', '--ratio', action="store_true",
                            help="Include ratio subplot in figure.")
        parser.add_argument('-s', '--scale', default=0, required=False, nargs='?', type=int,
                            choices=range(16), metavar='[0-15]',
                            help='For flexible-scale tables define central scale choice for MuR and MuF by selection enum fastNLO::ScaleFunctionalForm ("0"=kScale1, "1"=kScale2, "2"=kQuadraticSum), ...')
        parser.add_argument('-v', '--verbose', action="store_true",
                            help="Increase output verbosity.")
        parser.add_argument('-x', '--xbin', default=-1, type=int,
                            help='Observable bin no. to plot. Note: There are always ALL bins evaluated in one go.')

        # Parse arguments
        args = vars(parser.parse_args())

        # List of table names
        files = args['table']
        print '\n'
        print '[fastnnlo_scaledep]: Analysing table list: '
        for file in files:
                print '[fastnnlo_scaledep]:   ', file

        # PDF set name
        pdfset = os.path.basename(args['pdfset'])
        print '[fastnnlo_scaledep]: Using PDF set', pdfset

        # Orders to be shown
        iorders = []
        iordmin = _text_to_order['LO']
        iordmax = _text_to_order['NNLO']
        if args['order'] is None:
                print '[fastnnlo_scaledep]: Evaluate table up to highest available order.'
        else:
                for ord in args['order']:
                        if _text_to_order.has_key(ord):
                                iorders.append(_text_to_order[ord])
                        else:
                                print '[fastnnlo_scaledep]: Illegal order specified, aborted!'
                                print '[fastnnlo_scaledep]: Order list:', args['order']
                                exit(1)
                iordmin = min(iorders)
                iordmax = max(iorders)
                print '[fastnnlo_scaledep]: Evaluate table for order(s)', args['order']

        # Check existence of NNLOJET dat file for each order if desired
        datfilenames = []
        if args['datfiles'] is None:
                print '[fastnnlo_scaledep]: No statistical uncertainties requested.'
        else:
            for dat in args['datfiles']:
                fname = dat+'.dat'
                lstat = os.path.isfile(fname)
                if lstat:
                    datfilenames.append(fname)
                else:
                    print '[fastnnlo_scaledep]: Given file ', fname, 'for statistical uncertainties not found, aborted!'
                    exit(1)

        # Scale choice
        scale_choice = args['scale']

        # Given filename
        given_filename = args['filename']

        # Plot formats to use
        formats = args['format']
        if formats is None: formats = ['png']
        for fmt in formats:
            if not _formats.has_key(fmt):
                print '[fastnnlo_scaledep]: Illegal format specified, aborted!'
                print '[fastnnlo_scaledep]: Format list:', args['format']
                exit(1)

        # x axis range specification
        binpower  = args['binpower']
        logpoints = args['logpoints']

        # Include ratio plots
        ratio = args['ratio']

        # Verbosity
        verb = args['verbose']

        # x bin
        xbin = args['xbin']
        if xbin==0:
            print '[fastnnlo_scaledep]: Illegal observable bin no. specified, aborted!'
            print '[fastnnlo_scaledep]: The first bin has no. 1 here!'
            print '[fastnnlo_scaledep]: xbin:', args['xbin']
            exit(1)

        # Loop over table list
        for table in files:
                print '[fastnnlo_scaledep]: Analysing table: ', table
                # Get rid of extensions (.tab.gz or .tab)
                tablename = os.path.splitext(os.path.basename(table))[0]
                tablename = os.path.splitext(tablename)[0]

                ###################### Start EVALUATION with fastNLO library ###################################################
                # SetGlobalVerbosity(0) # Does not work since changed to default in the following call
                # Take the general information (bin_bounds, labels, order_existence, etc.) from given pdfset.
                fnlo = fastnlo.fastNLOLHAPDF(table, args['pdfset'], args['member'])

                # Get no. of observable bins to plot
                nobs = fnlo.GetNObsBin()
                if verb: print '[fastnnlo_scaledep]: NObsBins:', nobs

                # Get x section description
                ylabel = fnlo.GetXSDescr()
                if verb: print '[fastnnlo_scaledep]: X Section:', ylabel

                # Get plot labeling
                # Dimensionality of the table:
                ndim = fnlo.GetNumDiffBin()
                if verb:
                    print '[fastnnlo_scaledep]: Dimensions:', ndim

                # Labels of all the dimensions:
                labels = fnlo.GetDimLabels()
                if verb: print '[fastnnlo_scaledep]: Labels:', labels

                # Get multidimensional binning
                borders = []
                Dim0BinBounds = fnlo.GetDim0BinBounds()
                if verb: print '[fastnnlo_scaledep]: Dim0BinBounds', Dim0BinBounds
                for i0 in range(len(Dim0BinBounds)):
                    if ndim == 1:
                        borders.append( [Dim0BinBounds[i0]] )
                    else:
                        Dim1BinBounds = fnlo.GetDim1BinBounds(i0)
                        if verb: print '[fastnnlo_scaledep]: Dim1BinBounds',Dim1BinBounds
                        for i1 in range(len(Dim1BinBounds)):
                            if ndim == 2:
                                borders.append( [Dim0BinBounds[i0], Dim1BinBounds[i1]] )
                            else:
                                Dim2BinBounds = fnlo.GetDim2BinBounds(i0,i1)
                                if verb: print '[fastnnlo_scaledep]: Dim2BinBounds',Dim2BinBounds
                                for i2 in range(len(Dim2BinBounds)):
                                    if ndim == 3:
                                        borders.append( [Dim0BinBounds[i0], Dim1BinBounds[i1], Dim2BinBounds[i2]] )
                                    else:
                                        print '[fastnnlo_scaledep]: Invalid no. of dimensions. Aborted!'
                                        print '[fastnnlo_scaledep]: ndim = ', ndim
                                        exit(1)

                # Creating x-axis
                # Range in mur scale factors to show with step size
                lbxmin = -float(binpower)
                lbxmax = +float(binpower)
                stepsize = (lbxmax-lbxmin)/(logpoints-1)
                lbxmur = np.arange(lbxmin,lbxmax+1e-6,stepsize)
                xmin = np.power(2,lbxmin)
                xmax = np.power(2,lbxmax)
                xmur = np.power(2,lbxmur)
                i = -1
                xind = [-1,-1,-1]
                for xmu in xmur:
                    i+=1
                    if np.abs(xmu-0.5)/0.5 < 1.e-6: xind[0] = i
                    if np.abs(xmu-1.0)/1.0 < 1.e-6: xind[1] = i
                    if np.abs(xmu-2.0)/2.0 < 1.e-6: xind[2] = i

                # muf scale factors to use for muf uncertainty band
                xmuf = np.array([0.5, 1.0, 2.0])
                if verb:
                    print '[fastnnlo_scaledep]: mur scale factors:', xmur
                    print '[fastnnlo_scaledep]: muf scale factors:', xmur

                x_axis = xmur
                if verb:
                    print '[fastnnlo_scaledep]: xmin =', xmin, ', xmax =', xmax

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
                                print '[fastnnlo_scaledep]: Invalid choice of scale = ', scale_choice, ' Aborted!'
                                print '[fastnnlo_scaledep]: For fixed-scale tables only the default=0 is allowed.'
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
                    print '[fastnnlo_scaledep]: Table has continuous orders up to', cnt_order, 'and a maximal order of', max_order

                # If nothing requested, set to cnt_order for fixed-scale tables and max_order for flex-scale tables
                if args['order'] is None:
                    if lflex:
                        iordmax = max_order
                    else:
                        iordmax = cnt_order
                        if iordmax < 0:
                            print '[fastnnlo_scaledep]: Noncontinuous order availability for fixed-scale table. Aborted!'
                            print '[fastnnlo_scaledep]: Fixed-scale tables require presence of all orders to allow scale variations.'
                            print '[fastnnlo_scaledep]: Available orders are:', o_existence
                            exit(1)
                order_list = []
                if args['order'] is None:
                    for iord in range(iordmax+1):
                        if o_existence[iord]:
                            order_list.append(_order_to_text[iord])
                else:
                    order_list = args['order']

                print '[fastnnlo_scaledep]: List of requested orders:', order_list

                # Read in statistical uncertainty for each order if requested
                dxsr = []
                lstat = (len(datfilenames) > 0)
                if lstat and len(datfilenames) != len(order_list):
                    print '[fastnnlo_scaledep]: Mismatch between no. of requested orders and no. of filenames for statistical uncertainties, aborted!'
                    exit(1)
                for fname in datfilenames:
                    print '[fastnnlo_scaledep]: Taking statistical uncertainties from', fname
                    cols = np.loadtxt(fname,usecols=range(3,5))
                    xs_dat   = np.array(cols[:,0])
                    dxs_dat  = np.array(cols[:,1])
                    dxsr_dat = np.divide(dxs_dat, xs_dat, out=np.ones_like(dxs_dat), where=xs_dat!=0)
                    dxsr.append(dxsr_dat)
                dxsr_cn = abs(np.array(dxsr))

                # For flexible-scale tables set scale to user choice (default is 0)
                if lflex:
                    print '[fastnnlo_scaledep]: Setting requested scale choice for flexible-scale table:', scale_choice
                    fnlo.SetMuRFunctionalForm(scale_choice)
                    fnlo.SetMuFFunctionalForm(scale_choice)
                else:
                    if scale_choice == 0:
                        print '[fastnnlo_scaledep]: Evaluating fixed-scale table. Scale choice must be', scale_choice
                    else:
                        print '[fastnnlo_scaledep]: No scale choice possible for fixed-scale table. Aborted!'
                        print '[fastnnlo_scaledep]: scale_choice = ', scale_choice
                        exit(1)

                # Get bin bounds in 1st dimension
                bounds = []
                bounds.append(fnlo.GetDim0BinBounds())
                dim0_bounds = np.squeeze(np.array(bounds),axis=0)

                # Now evaluate fastNLO table having a look at scale uncertainties
                xsfls = []
                xss   = []
                xsfus = []
                for n in order_list:
                    for j in range(0, max_order+1):
                        if j <= _text_to_order[n]:
                            fnlo.SetContributionON(fastnlo.kFixedOrder, j, True)
                        else:
                            fnlo.SetContributionON(fastnlo.kFixedOrder, j, False)
                        if verb:
                            print '[fastnnlo_scaledep]: \n'
                            print '[fastnnlo_scaledep]: Calculate XS for order: %s' %n, '\n'
                            print '[fastnnlo_scaledep]: ----  ----  ----  ----  ----  ----  ----  ----'

                    xsfl = []
                    xs   = []
                    xsfu = []
                    for xmr in xmur:
                        # If xmuf not available for fixed-scale table, try to use HOPPET
                        lsclok = fnlo.SetScaleFactorsMuRMuF(xmr,xmuf[0])
                        if not lsclok: fnlo.UseHoppetScaleVariations(True)
                        fnlo.CalcCrossSection()
                        xsfl.append(fnlo.GetCrossSection())
                        lsclok = fnlo.SetScaleFactorsMuRMuF(xmr,xmuf[1])
                        if not lsclok: fnlo.UseHoppetScaleVariations(True)
                        fnlo.CalcCrossSection()
                        xs.append(fnlo.GetCrossSection())
                        lsclok = fnlo.SetScaleFactorsMuRMuF(xmr,xmuf[2])
                        if not lsclok: fnlo.UseHoppetScaleVariations(True)
                        fnlo.CalcCrossSection()
                        xsfu.append(fnlo.GetCrossSection())

                        #                if debug:
                        #                        print 'xsfl',xsfl
                        #                        print 'xs',xs
                        #                        print 'xsfu',xsfu

                    xsfls.append(xsfl)
                    xss.append(xs)
                    xsfus.append(xsfu)

                    if verb:
                        print '[fastnnlo_scaledep]: \n'
                        print '[fastnnlo_scaledep]: Relative scale uncertainty in %s: \n'%n
                        print '---------------------------------------------------------------------------------------'
                        print '---------------------------------------------------------------------------------------'

                xs_cn = np.array(xss)
                xs_fl = np.array(xsfls)
                xs_fu = np.array(xsfus)

                if verb:
                    print '[fastnnlo_scaledep]: Cross section xs_cn uses ', order_list

                if verb:
                    print xs_cn, '\n \n'

                ############################## Do the plotting ####################################################

                basename   = ''
                ordernames = ''
                for order in order_list:
                    ordernames += '_%s' %order
                if given_filename is not None:
                    basename = '%s_scaledep_%s.%s.%s' %(given_filename, pdfset, ordernames[1:], scale_name)
                else:
                    basename = '%s-scaledep-%s.%s.%s' %(tablename, pdfset, ordernames[1:], scale_name)

                    # Without statistical uncertainties create zero array of proper dimensions here
                if len(dxsr_cn)==0:
                    dxsr_cn = np.zeros((xs_cn.shape[0],xs_cn.shape[2]))
                if _debug: print 'dxsr_cn',dxsr_cn

                for iobs in range(nobs):
                    if xbin==-1 or xbin==iobs+1:
                        filename = '%s.bin%s' %(basename, str(iobs+1).zfill(2))
                        plotting(x_axis, xmin, xmax, iobs, xs_cn[:,:,iobs], xs_fl[:,:,iobs], xs_fu[:,:,iobs], dxsr_cn[:,iobs], xind, tablename, order_list, filename, scale_name, pdfset, labels, ylabel, borders, ratio, formats, verb)

if __name__ == '__main__':
        main()
