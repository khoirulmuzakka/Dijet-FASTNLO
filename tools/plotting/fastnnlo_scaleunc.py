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

#function for plotting uncertainties for one single order
def plotting_single(order_index, x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_choice):
	order_color={'LO':'g', 'NLO':'b', 'NNLO':'r'}
	scale_name={0:'kScale1', 1:'kScale2'}
	if variation_type=='scale uncertainty (2P)':
		vartype='2P'
	elif variation_type=='scale uncertainty (6P)':
		vartype='6P'

	gs = gridspec.GridSpec(3,3)
	fig = plt.figure(figsize=(7,7))
	ax1 = plt.subplot(gs[:-1,:])

	order_name = order_list[order_index]
	#ax1.loglog(x_axis, xs_all[order_index], 'gd', alpha=0.1)
	ax1.errorbar(x_axis, xs_all[order_index], #xerr=x_errors, 
			yerr=abs(abs_scale_unc[order_index]), elinewidth=1, linewidth=0.0, ms=4, color=order_color[order_name], fmt='.', alpha=0.9, label=order_name)
	
	#ax1.axis()
	ax1.set_xscale('log', nonposx='clip')
	ax1.set_yscale('log', nonposy='clip')
	ax1.set_xlabel('%s' %xlabel)
	ax1.set_ylabel('XS with abs_scale_unc', rotation=90)
	ax1.legend(fontsize=10, numpoints=1)
	ax1.text(0.02, 0.08, '%s' %variation_type, horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)
	ax1.text(0.02, 0.04, '%s' %scale_name[scale_choice], horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)	

	ax1.set_title('%s' %tablename)
		
	#subplot with relative scale uncertainties
	ax2 = plt.subplot(gs[2, :])
	
	ax2.semilogx(x_axis, xs_all[order_index]/xs_all[order_index], ls='solid', linewidth=0.02, color='red', alpha=0.8) #is this necessary?

	ax2.fill_between(x_axis, 1+rel_scale_unc[order_index, 2, :], 1+rel_scale_unc[order_index, 1, :], color=order_color[order_name], alpha=0.6, label=order_name)
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
	print 'saved as: %s.png' %filename

#function for plotting multiple orders into one figure
def plotting_multiple(lowest_order, x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_choice):
	order_color={'LO':'g', 'NLO':'b', 'NNLO':'r'}
	index_dict = {'LO':0, 'NLO':1, 'NNLO':2}
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
	for order_item, shift in zip(order_list, shift_list):
		order_index = index_dict[order_item]
		ax1.errorbar(x_axis*shift, xs_all[order_index], yerr=abs(abs_scale_unc[order_index]), elinewidth=1, linewidth=0.0, ms=4, color=order_color[order_item], fmt='.', label=order_item)

	ax1.set_xscale('log', nonposx='clip')
	ax1.set_yscale('log', nonposy='clip')
	ax1.set_xlabel('%s' %xlabel)
	ax1.set_ylabel('XS with abs_scale_unc', rotation=90)
	ax1.legend(fontsize=10, numpoints=1)
	ax1.text(0.02, 0.08, '%s' %variation_type, horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)
	ax1.text(0.02, 0.04, '%s' %scale_name[scale_choice], horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes)	


	ax1.set_title('%s' %tablename)
	

	#subplot with relative scale uncertainties, denominator in ratio = lowest given order
	ax2 = plt.subplot(gs[2, :])
	
	ordernames = ''
	for order_item in order_list:
		order_index = index_dict[order_item]
		ax2.semilogx(x_axis, xs_all[order_index]/xs_all[lowest_order], '.', ms=4, ls='solid', linewidth=0.2, color=order_color[order_item], label=order_item)
		#ax2.fill_between(xs_axis, 1+rel_scale_unc[order_index, 2, :], 1+rel_scale_unc[order_index, 1, :], color=order_color[order_item], alpha=0.6)
		## denominator in ratio is lowest order in order_list
		ax2.fill_between(x_axis, (xs_all[order_index]/xs_all[lowest_order])+rel_scale_unc[order_index, 2, :], 
				(xs_all[order_index]/xs_all[lowest_order])+rel_scale_unc[order_index, 1, :], color=order_color[order_item], alpha=0.46)


		#ax2.axis([xmin, xmax, 1+1.1*min(rel_scale_unc[order_index, 2, :]), 1+1.1*max(rel_scale_unc[order_index, 1, :])])
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
	print 'saved as: %s.png' %filename

def main():
	parser = argparse.ArgumentParser()

	parser.add_argument('table', type=str, help='FastNLO table that shall be evaluated.') #table is always required

	parser.add_argument('-p', '--pdfset', default='CT14nlo', 
				help='PDFset to evaluate fastNLO table.')
	parser.add_argument('-m', '--member', default=0, type=int,
				help='Member of PDFset, default is 0.')
	parser.add_argument('-o', '--order', required=False, nargs='*', type=int, 
				help='Order(s) for which scale uncertainty is investigated. LO=0, NLO=1, NNLO=2.'
				'If nothing is chosen, produce plots for all orders that are available in table.')

	#parser.add_argument('-v', '--variation', default=0, type=int,
	#			help='Chose between symmetric or asymmetric variations, default is symmetric.')

	parser.add_argument('-a', '--asymmetric', default=False, const=True, type=bool, nargs='?',
				help='If -a is chosen, use asymmetric scale variations, default is symmetric.')
	parser.add_argument('-s', '--scale', default=0, required=False, nargs='?', type=int, 
				help='Set scales MuR and MuF by choosing 0 (kScale1) or 1 (kScale2).')
	parser.add_argument('-f', '--filename') 


	#parse arguments
	args = vars(parser.parse_args())

	#table name
	tablename = os.path.splitext(os.path.basename(args['table']))[0]
	tablename = os.path.splitext(tablename)[0] #to get rid of extension (.tab.gz or .tab)
	print '\n'
	print 'Table: ', tablename, '\n'

	#pdfset name
	pdfset = os.path.basename(args['pdfset'])
	print 'PDF Set: ', pdfset, '\n'
	

	#order --> further evaluation, see evaluation part
	order = args['order']
	print 'Specific order(s)?: ', order
	

	#type of scale variation (symmetric vs asymmetric)
	if args['asymmetric']==True:
		scale_var_type = fastnlo.kAsymmetricSixPoint
		variation_type = 'scale uncertainty (6P)'
		#variation_type = 'asymmetric variation (6 point)'
	else:
		scale_var_type = fastnlo.kSymmetricTwoPoint
		variation_type = 'scale uncertainty (2P)'
		#variation_type = 'symmetric variation (2 point)'

	#chosen scale
	scale_choice = args['scale']

	#given filename
	given_filename = args['filename']


	###################### Start EVALUATION with fastNLO library #############################################
	fnlo = fastnlo.fastNLOLHAPDF(args['table'], args['pdfset'], args['member'])

	#Dictionary containing the fastNLO settings for certain orders
	orders = {'LO':[True, False, False], 'NLO':[True, True, False], 'NNLO':[True, True, True]}

	#Get labeling for the x-axis
	#dimensionality of the table:
	ndim = fnlo.GetNumDiffBin()
	print '\n', 'Dimensions: ', ndim

	#labels of all the dimensions:
	labels = fnlo.GetDimLabels()
	print 'Labels: ', labels

	#label of first dimension:
	xlabel = fnlo.GetDimLabel(0)
	print 'x-label: ', xlabel

	#creating x-axis
	bin_bounds = np.array(fnlo.GetObsBinsBounds(0))
	print 'bin_bounds.T: \n', bin_bounds.T, '\n'
	print 'bin_bounds.flatten()', bin_bounds.flatten(), '\n'
	
	x_axis = (bin_bounds.T[0]+bin_bounds.T[1])/2. #this is a list of bin centers
	xmin = min(bin_bounds.ravel())-0.1*min(bin_bounds.ravel())
	xmax = max(bin_bounds.ravel())+0.1*max(bin_bounds.ravel())
	
	#preparing x-errors (via bin_bounds) --> x_errors[0, :] are initially negative (positive via -1*), x_errors[1, :] positive
	x_errors = np.array([-1*(bin_bounds.T[0]-x_axis), bin_bounds.T[1]-x_axis])
	#print '\n x_errors: ', x_errors, '\n'



	##### Check: Which orders shall be investigated? Does table contain NNLO? ######
	#check existence of nnlo entry in table: ##
	nnlo_existence = fnlo.SetContributionON(fastnlo.kFixedOrder, 2, True) #returns True or False depending on whether NNLO exists
	
	index_dict = {'LO':0, 'NLO':1, 'NNLO':2}
	order_list = []
	if (args['order'] is not None):
		for given_order in args['order']:
			if given_order not in [0, 1, 2]:
				sys.exit('Invalid choice of orders. Aborted! \n'
					 'Possible choices are 0, 1 and 2. Order which caused problems: %s' %given_order)
			elif given_order==0:
				order_list.append('LO')
			elif given_order==1:
				order_list.append('NLO')
			elif given_order==2:
				if nnlo_existence==True:
					order_list.append('NNLO')
				else: 	
					print "No NNLO entry in given table."
					sys.exit("Aborted! Choose different orders (0, 1) or None for plotting all available orders.")

		##if 0 in args['order']:
		##	order_list.append('LO')
		##if 1 in args['order']:
		##	order_list.append('NLO')


	elif (args['order'] is None):
		if (nnlo_existence==False):
			print "No NNLO entry in given table."
			order_list = ['LO', 'NLO']
		elif (nnlo_existence==True):
			print "Table contains NNLO entry."
			order_list = ['LO', 'NLO', 'NNLO']


	print "order_list: ", order_list


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
		for j in range(0, 3):
			fnlo.SetContributionON(fastnlo.kFixedOrder, j, orders[n][j])
		print '\n'
		print 'Calculate XS for order: %s' %n, '\n'
		print '----  ----  ----  ----  ----  ----  ----  ----'
		fnlo.CalcCrossSection()
		xs_list.append(fnlo.GetCrossSection())

		### Get scale uncertainties ###
		print 'Used scale factor MuF: ', fnlo.GetScaleFactorMuF()
		print 'Used scale factor MuR: ', fnlo.GetScaleFactorMuR(), '\n'
		print 'Calculate scale uncertainties \n'
		## RELATIVE scale uncertainty with chosen type of scale variation (symmetric or asymmetric)
		rel_scale_unc_item = np.array(fnlo.GetScaleUncertaintyVec(scale_var_type)) #calculate this already for all accessible orders in any case 
		rel_unc_list.append(rel_scale_unc_item)
		print '\n'
		print 'Relative scale uncertainty in %s: \n'%n
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

	if nnlo_existence==False:
		## xs_all = np.delete(xs_all, -1, 0)  #because the last entry ('NNLO xs') would just be a duplicate of the NLO xs
		print 'Cross section xs_all in LO, NLO: \n'
	elif nnlo_existence==True:
		print 'Cross section xs_all in LO, NLO, NNLO: \n'
	
	print xs_all, '\n \n'




	## ABSOLUTE scale uncertainty
	num_orders = np.size(xs_all, 0) #length of axis 0 in xs_all equals number of orders
	abs_scale_unc = np.empty([num_orders, 2, len(x_axis)])
	for k in range(0, len(xs_all)):
		abs_scale_unc[k, 0, :] = rel_scale_unc[k, 2, :]*xs_all[k] #absolute uncertainties downwards (low)
		abs_scale_unc[k, 1, :] = rel_scale_unc[k, 1, :]*xs_all[k] #absolute uncertainties upwards (high)

	print 'Absolute Scale uncertainties downwards, upwards (order by order): \n'
	print abs_scale_unc, '\n'





	############################## PLOTTING Scale Uncertainties ####################################################

	### Plotting procedure if specific order has been chosen via -o ###
	## Producing only one plot
	if (args['order'] is not None):
		if (len(order_list)==1):
			order_index = args['order']
			plotting_single(order_index, x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_choice) 
			#should some of these just be global variables?
		else: #plotting of combination plot which contains all requested orders (contains either (LO, NLO) or (LO, NNLO) or (NLO, NNLO) or (LO, NLO, NNLO))
			lowest_order = min(args['order']) #will be used as denominator in ratio plot (rel_scale_unc)
			plotting_multiple(lowest_order, x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_choice)

	
	### Plotting procedure to produce plots for all orders that are available ###
	elif (args['order'] is None):
		#for l in range(0, len(order_list)):
		for l in order_list:
			order_index = index_dict[l]
			plotting_single(order_index, x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_choice)
			###add also plot where all orders are plotted within the same figure

		lowest_order = 0 #because here all orders are considered --> 0 (=LO) is always the lowest.
		plotting_multiple(lowest_order, x_axis, xmin, xmax, xs_all, rel_scale_unc, abs_scale_unc, xlabel, tablename, order_list, variation_type, given_filename, scale_choice)



if __name__ == '__main__':
	main()
