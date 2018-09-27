#!/usr/bin/env python2
#-*- coding:utf-8 -*-

###################################################
#
# Script for creating plots of the kfactor.
# Processing one table at a time.
# Using the fastNLO library functionalities.
# Possible update: several tables within one plot?
#
# Created by B.Schillinger, 10.09.2018
#
###################################################

import argparse
import numpy as np
import os, sys
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.gridspec as grispec
from matplotlib.pyplot import cm
import matplotlib.pylab as pylab
import matplotlib.ticker #to fix x-axis ticks

import fastnlo

def main():
	parser = argparse.ArgumentParser()

	parser.add_argument('-t', '--table', default='table.tab', required=True, type=str,
				help='FastNLO table that shall be evaluated.')
	parser.add_argument('-p', '--pdfset', default='CT14nlo',
				help='PDFset to evaluate fastNLO table.')
	parser.add_argument('-m', '--member', default=0, type=int,
				help='Member of PDFset, default is 0.')

	#this option is not working yet, will be updated
	parser.add_argument('-n', '--numerator', default=None, type=int,
				help='The higher order for kfactor.')
	parser.add_argument('-d', '--denominator', default=None, type=int,
				help='The lower order for kfactor.')
	parser.add_argument('-f', '--filename', default=None, type=str,
				help='Output filename (optional).')


	#parse arguemnts
	args = vars(parser.parse_args())

	#table name
	tablename = os.path.basename(args['table'])
	#tablename = os.path.splitext(tablename)[0] #to get rid of extension (.tab.gz or .tab)
	print 'Table: ', tablename, '\n'

	#pdfset name
	pdfset = os.path.basename(args['pdfset'])
	pdfname = os.path.splitext(pdfset)[0]
	print 'PDF Set: ', pdfname, '\n'

	#chosen higher order and lower order for kfactor
	print 'Higher order: ', args['numerator'] #was earlier called args['order']
	print 'Lower order: ', args['denominator'] #used to be args['normoder']

	if args['numerator'] is not None:
		if (args['numerator'] in [1, 'nlo', 'NLO']):
			high = 'NLO'
		elif (args['numerator'] in [2, 'nnlo', 'NNLO']):
			high = 'NNLO'
		else:
			print 'ERROR: Invalid choice of order. Aborted.'
			sys.exit('Higher order can be NLO (=1) or NNLO (=2).')
		print 'Chosen highest order: %s.' %high
	else:
		print 'No specific order chosen.'
	
	if args['denominator'] is not None:
		if (args['denominator'] in [0, 'lo', 'LO']):
			low = 'LO'
		elif (args['denominator'] in [1, 'nlo', 'NLO']):
			low = 'NLO'
		else:
			print 'ERROR: Invalid choice of lower order. Aborted.'
			sys.exit('Lower (=denominator) order can be LO (=0) or NLO (=1).')
		print 'Chosen lower (=denominator) order: %s.' %low
	else:
		print 'No specific lower order chosen.'


	if (args['numerator'] is None) and (args['denominator'] is None):
		print 'No specific orders chosen: Look at NLO/LO, NNLO/NLO and NNLO/LO.'
		all_orders=['LO', 'NLO', 'NNLO']
	
	

	############### Start EVALUATION with fastNLO library ###################################
	fnlo = fastnlo.fastNLOLHAPDF(args['table'], args['pdfset'], args['member'])

	#Dictionary containing the fastNLO settings for certain orders
	orders = {'LO': [True, False, False], 'NLO': [True, True, False], 'NNLO': [True, True, True]}

	#Get labeling for the x-axis
	#dimensionality of table:
	ndim = fnlo.GetNumDiffBin()
	print '\n', 'Dimensions: ', ndim

	#labels of all the dimensions:
	labels = fnlo.GetDimLabels()
	print 'Labels: ', labels

	#label of first dimension:
	xlabel = fnlo.GetDimLabel(0)
	print 'x-label: ', xlabel


	# Now evaluate fastNLO table for creating 3 plots
	## if we allow options -o and -n, we will need an if-condition here!
	print 'Start table evaluation for creating three plots. \n'

	#cross section for all 3 orders
	xs_list = []
	for n in orders:
		for j in range(0,3):
			fnlo.SetContributionON(fastnlo.kFixedOrder, j, orders[n][j])
		print '\n'
		print 'Calc XS for order: %s' %n, '\n'
		fnlo.CalcCrossSection()
		xs_list.append(fnlo.GetCrossSection())
	xs_all = np.array(xs_list)
	print 'Cross section with all subprocesses xs_all: \n'
	print xs_all, '\n \n'

	#calculate k-facators
	fractions_three_plots = []
	for i in range(0, 2): #take care of NLO/LO and NNLO/NLO
		fractions_three_plots.append(np.divide(xs_all[i+1], xs_all[i]))
	fractions_three_plots.append(np.divide(xs_all[2], xs_all[0])) #here NNLO/LO
	
	#kfactors: NLO/LO, NNLO/NLO, NNLO/LO
	fractions_three = np.array(fractions_three_plots)
	print 'k-factors (line by line): NLO/LO, NNLO/NLO, NNLO/LO: \n', fractions_three


	################ Start PLOTTING of kfactor #############################################
	#set parameters
	params = {'xtick.labelsize':'large',
			'xtick.major.size': 5,
			'xtick.major.width': 2,
			'xtick.minor.size': 3.5,
			'xtick.minor.width': 1,
			'ytick.labelsize':'large',
			'ytick.major.size': 5}
	pylab.rcParams.update(params)

	plt.close()
	fig0 = plt.figure(figsize=(8,7))
	ax0 = fig0.add_subplot(111)

	color = iter(cm.rainbow(np.linspace(0,1,len(fractions_three))))
	colors_orders = {'NLO/LO':'r', 'NNLO/NLO':'b', 'NNLO/LO':'g'}

	## axis settings? color settings?
	bin_bounds = np.array(fnlo.GetObsBinsBounds(0))
	print 'bin bounds: ', bin_bounds.flatten()
	xstart = bin_bounds.flatten()[0]
	xstop = bin_bounds.flatten()[-1]
	xaxis_ticks = np.linspace(xstart, xstop, num=5)

	# labels
	labels = ['NLO/LO', 'NNLO/NLO', 'NNLO/LO']

	## plot all the kfactors into one plot. can be changed to 3 plots by introducing some for-loop
	print 'Start plotting. \n'
	for k in range(0, 3):
		print 'Current index in fractions array: ', k
		#c = next(color)
		c = colors_orders[labels[k]]
		ax0.plot(bin_bounds.flatten(), steppify_bin(fractions_three[k]),
			label=labels[k], color=c, alpha=1.0)
	
	#get limits for axes:
	xlim = ax0.get_xlim()
	ylim = ax0.get_ylim()

	print 'xlim: ', xlim
	print 'ylim: ', ylim, '\n'

	#settings for the whole plot
	ax0.set_xscale('log', nonposx='clip')
	ax0.axis([xlim[0], xlim[1], ylim[0]-0.05, ylim[1]+0.05]) #flexible axis
	ax0.set_xlim(xlim[0], xlim[1]) #adjust x-axis
	#plt.xticks(xaxis_ticks, xaxis_ticks
	ax0.autoscale(enable=True, axis='x', tight=True) #so that only xrange with data is plotted
	## The following line is only working with matplotlib version 2.0.2 or higher (just some x-axis settings, might be removed eventually)
	##matplotlib.ticker.LogLocator(base=10.0, subs='all') #to show all xticks, also the lowest

	plt.axhline(y=1, xmin=0, xmax=1, color='k', linestyle='dotted')
	plt.axhline(y=0, xmin=0, xmax=1, color='k', linestyle='dotted')

	plt.title(x=0.5, y=1.06, s='%s' %tablename, fontsize=16) #, tight=True)

	ax0.set_xlabel('%s' %xlabel, fontsize=12, horizontalalignment='right')
	#ax0.xaxis.set_label_coords(0.98, -0.06)
	ax0.xaxis.set_label_coords(1.00, -0.06)
	ax0.set_ylabel('k-factor', rotation=90, fontsize=12)
	#ax0.yaxis.set_label_coords()

	ax0.text(0.96, 0.03, args['pdfset'], alpha=0.6, transform=ax0.transAxes, fontsize=14, verticalalignment='bottom', horizontalalignment='right')
	plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=12)

	plt.tight_layout()
	fig0.tight_layout()

	#naming of the plot		##include more information on process etc here!!
	if args['filename'] is not None:
		plotname = '%s_kfactor.png' %args['filename']
	else:
		plotname = '%s.kfactor.png' %tablename
	
	fig0.savefig(plotname, bbox_inches='tight')
	print 'k-factor plot saved as: %s' %plotname,'\n \n'



########################################################
## function for producing stepped arrays for plotting
def steppify_bin(arr, isx=False):
	"""
	Produce stepped array of arr, needed for example for stepped fill_betweens.
	Pass all x bin edges to produce stepped x arr and all y bincontents to produce
	stepped bincontents representation
	steppify_bin([1,2,3], True)
	-> [1,2,2,3]
	steppify_bin([5,6])
	-> [5,5,6,6]
	"""
	if isx:
		newarr = np.array(zip(arr[:-1], arr[1:])).ravel()
	else:
		newarr = np.array(zip(arr, arr)).ravel()
	return newarr
#########################################################


if __name__ == '__main__':
	main()
