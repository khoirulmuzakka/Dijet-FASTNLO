#! /usr/bin/env python

###############################################################
#Script for comparing the contribution of several subprocesses
#to the total cross section. 
#Either (per default) as stackplot or for single process (-e).
#
#Choice of specific orders for subprocess (-o) and normalisation
#(-n) is possible here.
#If no specific order is chosen --> all 5 plots are produced.
#
#General script, not customized for any certain tables.
###############################################################


import argparse
import shutil
import numpy as np
import os
import time, datetime
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.pyplot import cm

import matplotlib.style
matplotlib.style.use('classic')			#to produce plots that look like they used to in matplotlib 1.x 

import fastnlo



def main():
	
	parser = argparse.ArgumentParser()
	
	#make it possible to use different options when running
	parser.add_argument('-t','--table', default='table.tab', required=True,
						help='FastNLO tables that shall be evaluated.')
	parser.add_argument('-p', '--pdfset', default='CT14nlo',
						help='PDFset(s) to evaluate fastNLO tables.')
	parser.add_argument('-m', '--member', default=0, type=int,
						help='Member of PDFset, default is 0.')
						
	parser.add_argument('-s', '--subproc', default=['all'], type=str, nargs='*',
						help='Subprocesses that shall be compared. \n'
							 'Options: gg, gq, gaq, qq, q=q, q!q, uu, auu, auau, etc. \n'
							 'See SelectProcesses() in fastNLOReader.cc for more information.')
	
	parser.add_argument('-e', '--extraplot', default=False, const=True, type=bool, nargs='?',
						help='Produces extra single-process-plots if option -e is chosen. \n'
							'Otherwise stackplot (containing all chosen processes) per default.')
						
	parser.add_argument('-f', '--filename', type=str, #default='stackplot.png', 
						help='Optional name for output file. \n'
						'If nothing is chosen: string containing tablename and orders.')
						
	parser.add_argument('-x', '--xaxis', default='x quantity', type=str,
						help='Set label for x-axis. Default is \'x quantity\'.')
	
	parser.add_argument('-c', '--constyaxis', default=False, const=True, type=bool, nargs='?',
						help='Set constant limits for y-axis from -0.4 to 1.2. \n'
						'Otherwise (per default) y-axis will be flexibly adjusted to the data.')
	
	
	#flexible choice of order for the subprocess contribution (-o) and normalisation (-n)
	#(NLO to NLO+LO for instance)
	parser.add_argument('-o', '--order', type=int,
						help='Order in which the (single) subprocess contributions to XS are calculated. \n'
								'Type 0 for LO, 1 for NLO, 2 for NNLO, etc. '
								'If nothing is chosen order=normorder. '
								'If normorder is not chosen either, all five plots are produced.')
	
	parser.add_argument('-n', '--normorder', type=int,
						help='Normalisation order in which the XS with all subprocesses included is calculated. \n'
								'Type 0 for LO, 1 for NLO, 2 for NNLO, etc. '
								'If nothing is chosen normorder=order. '
								'If order is not chosen either, all five plots are produced.')


################# TAKING CARE OF THE INPUT #################################
	#parse arguments
	args = vars(parser.parse_args())
	
	#table name
	table_ = os.path.basename(args['table'])
	tablename = os.path.splitext(table_)[0]
	tablename = os.path.splitext(tablename)[0] #because yb tables have .tab.gz ending (getting rid of both here)
	print '\n'
	print 'Table: ', tablename, '\n'
	
	#pdfset name
	pdfset_ = os.path.basename(args['pdfset'])
	pdfname = os.path.splitext(pdfset_)[0]
	print 'PDF Set: ', pdfname, '\n'
	
	#Which subprocesses are evaluated?
	print 'Subprocesses that are investigated: '
	print args['subproc']


	#handling of subprocess order and normalisation order:

	#for labeling purposes
	print 'subprocess order: ', args['order']
	print 'normalisation order: ', args['normorder']

	#subprocess order (either chosen or per default all combinations)
	sp_order_all = ['LO', 'NLO', 'NNLO']
	if (args['order'] is not None) and (args['order'] in range(0,3)):
		sp_order = sp_order_all[args['order']]
	elif (args['order'] is None) and (args['normorder'] is not None):
		sp_order = sp_order_all[args['normorder']]
		print "Set subprocess order to chosen normalisation order."
	else:
		print "No certain (or valid) subprocess order given."
		sp_order = sp_order_all #in case of no specification --> sp_order is array, will produce all 5 plots
	
	#normalisation order (as above)
	norm_order_all = ['LO', 'NLO', 'NNLO']
	if (args['normorder'] is not None) and (args['normorder'] in range(0,3)):
		norm_order = norm_order_all[args['normorder']]
	elif (args['normorder'] is None) and (args['order'] is not None):
		norm_order = norm_order_all[args['order']]
		print "Set normalisation order to chosen subprocess order."
	else:
		print "No certain (or valid) normalisation order given."
		norm_order = norm_order_all #watch out! here norm_order becomes an array



############################ EVALUATING ####################################
	fnlo = fastnlo.fastNLOLHAPDF(args['table'], args['pdfset'], args['member'])
	
	#dictionary with settings for the different orders (when switching them on/off)
	b_order = { "LO": [True, False, False], "NLO": [True, True, False], "NNLO": [True, True, True] }

	'''
	if (norm_order < sp_order):
		print "Choice of orders not reasonable." #depends on purpose of the plot. (maybe useful for kfactors if args['subproc']='all')
	else:
		#usual plotting
	'''
	
	print "\n"
	#check whether all 5 plots are needed or just one certain --> evaluate accordingly
	if isinstance(norm_order, list) and isinstance(sp_order, list): #actually one condition to check should be enough in both cases
		#produce all 5 plots (this also happens if both input orders were invalid...)
		print 'Start table evaluation for creating 5 plots. \n'
		
		#cross section with every process included
		xs_all_list = []
		for n in norm_order:
			print "\n", "n: ", n
			for j in range(0,3):
				fnlo.SetContributionON(fastnlo.kFixedOrder, j, b_order[n][j]);
		#for order, setting in sorted(b_order.iteritems()):
		#	index = 0
		#	for j in setting:
		#		fnlo.SetContributionON(fastnlo.kFixedOrder, index, j)
		#		index+=1
			
			fnlo.CalcCrossSection()
			xs_all_list.append(fnlo.GetCrossSection())
		xs_all = np.array(xs_all_list)
		print '\n'
		print 'Cross section with all subprocesses xs_all: \n'
		print xs_all, '\n \n'
		
		#cross section of single processes
		xs_sub_list = []
		for s in sp_order:			#go through lo, nlo, nnlo
			print "\n", "Subprocess order s: ", s
			xs_subproc_list = [] #will be emptied after each iteration
			for l in range(0,3):	#go through the settings for this order s
				fnlo.SetContributionON(fastnlo.kFixedOrder, l, b_order[s][l])
		
			for i in range(0, len(args['subproc'])):
				subproc_name = args['subproc'][i]
				print 'Selected subprocess: %s' %subproc_name
				fnlo.SelectProcesses(subproc_name, True) #optional second argument: symmetric==True per default
				fnlo.CalcCrossSection()
				xs_subproc = np.array(fnlo.GetCrossSection())
				xs_subproc_list.append(xs_subproc) #append xs of chosen subprocess in that certain order
				print 'XS for subprocess %s: \n' %args['subproc'][i]
				print xs_subproc, '\n \n' #should be one line in xs_sub[k,i,:]
			#all selected processes in lo or nlo or nnlo
			xs_sub_list.append(np.array(xs_subproc_list))
		xs_sub = np.array(xs_sub_list)
		print "XS of the selected processes in different pTZ-regions, in LO, NLO and NNLO: "
		print "From outmost level to innermost: order, subprocess, pTZ-bin. \n"
		print xs_sub, '\n \n'


		#### Investigating the fractions ####
		#Now, check all five combinations of lo, nlo, nnlo that shall be plotted
		#How much do the selected subprocesses contribute to the total xs?

		#Three plots where subprocess is of same order as normalisation is
		fractions_eq_order = []
		for i in range(0,3):
			fractions_eq_order.append(np.divide(xs_sub[i,:],xs_all[i]))
		fractions_eq = np.array(fractions_eq_order)
		print 'fractions, equal order of subprocess and normalisation: \n', fractions_eq, '\n'

		'''
		#output for checking what has been calculated
		for n in range(0, 3):
			for i in range(0, len(fractions_eq[0,:])):
				print 'fraction of %s on total XS, both in %s: ' %(args['subproc'][i], sp_order_all[n])
				print fractions_eq[n,i,:], '\n \n'
		'''

		#make the two additional plots for (LO to NNLO) and (NLO to NNLO)
		fractions_diff_order = []
		for s in range(0,2):			#s is current subprocess order (LO, NLO)
			fractions_diff_order.append(np.divide(xs_sub[s,:,:], xs_all[2,:]))
		fractions_diff = np.array(fractions_diff_order)
		print 'fractions, different order for subprocess/normalisation: LO/NNLO, NLO/NNLO: '
		print fractions_diff, '\n'
		
		#combine these arrays to 'fractions'-array (5 outmost indices, one for each plot)
		fractions = np.append(fractions_eq, fractions_diff, axis=0)
		print "------------------------------------------------------------------"
		print 'Fractions for all the plots to be produced, first index indicates:'
		print 'LO/LO, NLO/NLO, NNLO/NNLO, LO/NNLO, NLO/NNLO'
		print fractions, '\n'
		print "------------------------------------------------------------------ \n"


		#checking whether contribution is positive or negative 
		#(process all 5 "lines" of fractions-array)
		pos_contr = np.array(fractions)
		neg_contr = np.array(fractions)
		bin_bounds = np.array(fnlo.GetObsBinsBounds(0))

		for p in range(0, len(fractions[:,0,0])): #should be range(0,6), as there are 5 plots
			for i in range(0, len(args['subproc'])):	#go through 'fractions' (second index line by line) = subprocess after subprocess
				for j in range(0, len(bin_bounds)):					#check the different pTZ bins
					if (fractions[p,i,j] < 0):
						pos_contr[p,i,j] = 0.0
					else:
						neg_contr[p,i,j] = 0.0
				#print 'Positive contributions of %s on total XS in plot %s: \n' %(args['subproc'][i], p)
				#print pos_contr[p,i,:], '\n'
				#print 'Negative contributions of %s on total XS in plot %s: \n' %(args['subproc'][i], p)
				#print neg_contr[p,i,:], '\n \n'

		print 'Summary of positive and negative contributions of the subprocesses for the 5 plots: \n'
		print 'pos_contr: \n', pos_contr, '\n'
		print 'neg_contr: \n', neg_contr, '\n'
		
	############################ Five-Plots-PLOTTING ############################################
		#Bins for x-axis (according to first pdf (0))
		bin_bounds = np.array(fnlo.GetObsBinsBounds(0))
		#print bin_bounds.T		#bin_bounds.T[0] = lower bounds, bin_bounds.T[1] upper bounds
		print 'bin bounds:', bin_bounds.flatten()
	
		#for naming the plots afterwards and give correct label to y-axis:
		plots_order = np.array([["LO", "LO"], ["NLO", "NLO"], ["NNLO", "NNLO"], ["LO", "NNLO"], ["NLO", "NNLO"]], dtype=str)
		print 'The five plots will contain: '
		print plots_order, '\n'

		#Check whether comparison-plot or single-plot is required:
		print 'Single plots for each process?: ', args['extraplot']
		if (args['extraplot']==False):		#stackplot that compares the contribution of the subprocesses to the total XS
			plt.close()

			for p in range(0, 5): #p is index that says which plot 0 to 5 is produced
				#plt.close('all')
				fig0 = plt.figure(figsize=(8,7))
				ax0 = fig0.add_subplot(111)
				#number of bins via bin_bounds variable
				y0_pos = np.zeros([len(bin_bounds)]) #for tracking the current 'height' of stackplot, in order to staple and not overlap the values for each subprocess
				y0_neg = np.zeros([len(bin_bounds)]) #same as above for the negative contributions, plotted below x-axis
				print '-	-	-	-'
				print 'Start plotting of plot %s \n' %p
				color = iter(cm.rainbow(np.linspace(0,1,len(args['subproc']))))
				patches = []  #later needed for legend

				for i in range(0, len(args['subproc'])):	#go through fractions, subproc by subproc
					c = next(color)
					#positive contributions
					ax0.fill_between(bin_bounds.flatten(), steppify_bin(y0_pos), steppify_bin(y0_pos+pos_contr[p, i]), #label=labeling[args['subproc'][i]],
									facecolor=c, alpha=0.4)
					y0_pos += pos_contr[p, i] #calculate new height of stackplot

					#negative contributions
					ax0.fill_between(bin_bounds.flatten(), steppify_bin(y0_neg), steppify_bin(y0_neg+neg_contr[p, i]), #label=labeling[args['subproc'][i]],
									facecolor=c, hatch='X', alpha=0.4)
					y0_neg += neg_contr[p, i] #new 'height' below x-axis

					#patch for subprocess i (for the legend) ##labeling happens here
					patches.append(matplotlib.patches.Rectangle((0, 0), 0, 0, color=c,
									label=args['subproc'][i], alpha=0.4))
					ax0.add_patch(patches[i])
			
				#get limits for axes:
				xlim = ax0.get_xlim()
				ylim = ax0.get_ylim()
			
				print 'xlim:', xlim
				print 'ylim:', ylim, '\n'
			
				xmin_ = xlim[0]
				print 'xmin', xmin_
				xmax_ = xlim[1]
				#ymin = ylim[0]
				#ymax = ylim[1]
				ymin_ = np.amin(y0_neg)+0.3*np.amin(y0_neg)
				ymax_ = np.amax(y0_pos)+0.1*np.amax(y0_pos)
				
				#settings for the whole stackplot
				ax0.set_xscale('log', nonposx='clip')
				#ax0.axis([30, 1000, -0.20, 1.20])
				#ax0.axis([bin_bounds.flatten()[0], xlim[1], ylim[0]-0.005, ylim[1]+0.005]) #flexible axes
				#ax0.axis([bin_bounds.flatten()[0], xlim[1], ylim[0]-0.01*ylim[0], ylim[1]+0.01*ylim[1]])
				
				if args['constyaxis']==True:
					ax0.axis([0, xmax, -0.20, 1.20])
				elif args['constyaxis']==False:
					ax0.set_xlim(xmin_, xmax_) ##
					ax0.set_ylim(ymin_, ymax_) ##set axis limit according to data

				##ax0.axis([0, xmax, -0.20, 1.20]) #for test with fixed y-axis
				ax0.autoscale(enable=True, axis='x', tight=True)
				plt.axhline(y=1, xmin=0, xmax=1, color='k', linestyle='dotted') #dotted line at xs_sub/xs=1=100%
				plt.axhline(y=0, xmin=0, xmax=1, color='k', linestyle='-')		#line at xs_sub/xs_tot=0
				
				plt.title(x=0.5, y=1.06, s='%s' %(table_), fontsize=22)

				ax0.set_xlabel('x quantity', fontsize=20)
				ax0.xaxis.set_label_coords(0.9, -0.05)
				ax0.set_ylabel(r'$\frac{\mathrm{XS_{sp, %s}}}{\mathrm{XS_{tot, %s}}}$' %(plots_order[p,0], plots_order[p,1]), rotation=0, fontsize=24)
				ax0.yaxis.set_label_coords(-0.16, 0.9)

				ax0.text(0.96, 0.03, args['pdfset']+', %s to %s' %(plots_order[p,0], plots_order[p,1]), alpha=0.6, transform=ax0.transAxes, fontsize=14, verticalalignment='bottom', horizontalalignment='right')
				plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left') #'upper left' refers to the legend box

				plt.tight_layout()
				fig0.tight_layout()

				#naming of the plots
				if args['filename'] is not None:
					stackplotname = '%s_%s.png' %(args['filename'], plots_order[p,0]+'_'+plots_order[p,1])
				else:
					stackplotname = '%s.stackplot_%s.png' %(tablename, plots_order[p,0]+'_'+plots_order[p,1])
				fig0.savefig(stackplotname, bbox_inches='tight')

				print 'Stackplot for subprocess comparison saved as: %s' %(stackplotname), '\n \n'


		else:
			print "----------------------------------------------------------"
			print "Create individual plot for each single subprocess \n"

			#individual plot for each subproc (NumberOfPlots = 5*NumberOfSubprocesses)
			#create folder for saving the plots, will be overwritten in case it already exists (avoid countless folders while testing)
			current_dir = os.getcwd()
			
			if args['filename'] is not None:
				prename = args['filename']
			else: prename = tablename

			final_dir = os.path.join(current_dir, r'%s_subproc' %(prename))

			if os.path.exists(final_dir):
				shutil.rmtree(final_dir) #remove final_dir if it already exists
			os.makedirs(final_dir)

			#loop where single plot for each of the selected subprocesses is created

			#plot the subprocess contributions for each chosen process individually
			#do this for all the five combinations: LO/LO, NLO/NLO, NNLO/NNLO, LO/NNLO, NLO/NNLO
			for p in range(0, 5): #p is index that says which plot 0 to 5 is produced
				#plt.close('all')
				color = iter(cm.rainbow(np.linspace(0,1,len(args['subproc']))))
				fig0 = plt.figure(figsize=(8,7))
				ax0 = fig0.add_subplot(111)
				
				for i in range(0, len(args['subproc'])):
					c = next(color)
					plt.close()
					fig0 = plt.figure(figsize=(8, 7))
					ax0 = fig0.add_subplot(111)

					#plot the fractions
					y0 = np.zeros([len(bin_bounds)])

					#ax0.fill_between(bin_bounds.flatten(), steppify_bin(y0), steppify_bin(fractions_oneplot[i]), 
					#			facecolor=c, alpha=0.4)
				
					ax0.fill_between(bin_bounds.flatten(), steppify_bin(y0), steppify_bin(pos_contr[p,i]), 
								facecolor=c, alpha=0.4)
					ax0.fill_between(bin_bounds.flatten(), steppify_bin(y0), steppify_bin(neg_contr[p,i]), 
								facecolor=c, hatch='X', alpha=0.4)

					#get limits for axes:
					xlim = ax0.get_xlim()
					ylim = ax0.get_ylim()
					
					
					ax0.set_xscale('log', nonposx='clip')
					#ax0.axis([30, 1000, -0.4, 1.2])
					
					if args['constyaxis']==True:
						ax0.axis([xlim[0], xlim[1], -0.2, 1.2])
					elif args['constyaxis']==False:
						ax0.axis([xlim[0], xlim[1], ylim[0]+0.3*ylim[0], ylim[1]+0.1*ylim[1]])
					
					#ax0.axis([bin_bounds.flatten()[0], xlim[1], ylim[0]-0.1*ylim[0], ylim[1]+0.1*ylim[1]]) #leave y-limits fixed for better comparison between the subprocess-plots
					ax0.autoscale(enable=True, axis='x', tight=True)
					plt.axhline(y=1, xmin=0, xmax=1, color='k', linestyle='dotted')
					plt.axhline(y=0, xmin=0, xmax=1, color='k', linestyle='-')


					titlename = args['subproc'][i]

					plt.title(x=0.5, y=1.01, s=titlename, fontsize=20)
					ax0.set_xlabel('x quantity', fontsize=20)
					ax0.xaxis.set_label_coords(0.9, -0.05)
					ax0.set_ylabel(r'$\frac{\mathrm{XS_{sp, %s}}}{\mathrm{XS_{tot, %s}}}$' %(plots_order[p,0], plots_order[p,1]), rotation=0, fontsize=24)
					ax0.yaxis.set_label_coords(-0.16, 0.9)
					ax0.text(0.96, 0.03, args['pdfset']+', %s to %s' %(plots_order[p,0], plots_order[p,1]), transform=ax0.transAxes, fontsize=14, verticalalignment='bottom', horizontalalignment='right')
					ax0.text(0.96, 0.08, '%s' %table_, transform=ax0.transAxes, fontsize=14, verticalalignment='bottom', horizontalalignment='right')

					#plt.legend() #location parameter loc='upper right'

					plt.tight_layout()
					fig0.tight_layout()
					plotname = '%s_%s_singlefrac_%s_%s.png' %(tablename, args['subproc'][i], plots_order[p,0], plots_order[p,1])

					#save the plots
					fig0.savefig(os.path.join(final_dir,plotname))
					print 'Plot saved in %s' %(final_dir)
					print 'fraction plot saved as: %s' %(plotname), "\n \n"



	elif isinstance(norm_order, str) and isinstance(sp_order, str):
		#produce only demanded plot for chosen orders
		print 'Start table evaluation for creating 1 plot. \n'

		print "Normalisation Order: ", norm_order
		for j in range(0,3):
			fnlo.SetContributionON(fastnlo.kFixedOrder, j, b_order[norm_order][j])
		fnlo.CalcCrossSection()
		xs_all_oneplot = np.array(fnlo.GetCrossSection())
		print "\n","Cross Section with all subprocesses included in %s: \n" %norm_order
		print xs_all_oneplot, "\n \n"
		
		print "Subprocess Order: ", sp_order
		for l in range(0,3):
			fnlo.SetContributionON(fastnlo.kFixedOrder, l, b_order[sp_order][l])
		xs_subproc_list = []
		for i in range(0, len(args['subproc'])):	#go through subprocesses
			subproc_name = args['subproc'][i]
			print 'Selected subprocess: %s' %subproc_name
			fnlo.SelectProcesses(subproc_name, True) #optional second argument: symmetric==True per default
			fnlo.CalcCrossSection()
			xs_subproc = np.array(fnlo.GetCrossSection())
			xs_subproc_list.append(xs_subproc)
			print "XS for subprocess %s in %s: \n" %(subproc_name, sp_order)
			print xs_subproc, '\n \n'
		xs_sub_oneplot = np.array(xs_subproc_list)
		print "XS of the selected processes in different pTZ-regions in %s: " %sp_order
		print xs_sub_oneplot, "\n \n"
		
		#calculate fraction of subprocess in sp_order on total XS in norm_order
		fractions_sub_oneplot = []
		for k in range(0, len(args['subproc'])):
			fractions_sub_oneplot.append(np.divide(xs_sub_oneplot[k,:], xs_all_oneplot[:]))
		fractions_oneplot = np.array(fractions_sub_oneplot)
		print "--------------------------------------------------------------"
		print "Fractions for subprocesses in %s to total XS in %s" %(sp_order, norm_order)
		print fractions_oneplot, "\n"
		print "--------------------------------------------------------------"

		#checking whether contribution is positive or negative 
		pos_contr = np.array(fractions_oneplot)
		neg_contr = np.array(fractions_oneplot)
		bin_bounds = np.array(fnlo.GetObsBinsBounds(0))

		for i in range(0, len(args['subproc'])):	#go through 'fractions' = subprocess after subprocess
			for j in range(0, len(bin_bounds)):					#check the different pTZ bins
				if (fractions_oneplot[i,j] < 0):
					pos_contr[i,j] = 0.0
				else:
					neg_contr[i,j] = 0.0
			print 'Positive contributions of %s in %s on total XS in %s: \n' %(args['subproc'][i], sp_order, norm_order)
			print pos_contr[i,:], '\n'
			print 'Negative contributions of %s in %s on total XS in %s: \n' %(args['subproc'][i], sp_order, norm_order)
			print neg_contr[i,:], '\n \n'

		print 'Summary of positive and negative contributions of the subprocesses: \n'
		print "positive: \n", pos_contr, '\n'
		print "negative: \n", neg_contr, '\n'


		######################## One-Plot-PLOTTING #########################################
		#Bins for x-axis (according to first pdf (0))
		bin_bounds = np.array(fnlo.GetObsBinsBounds(0))
		print 'bin bounds:', bin_bounds.flatten()

		#Check whether comparison-plot or single-plot is required:
		print 'Single plots for each process?: ', args['extraplot']
		if (args['extraplot']==False):		#stackplot that compares the contribution of the subprocesses to the total XS
			plt.close()

			fig0 = plt.figure(figsize=(8,7))
			ax0 = fig0.add_subplot(111)
			#number of bins via bin_bounds variable
			y0_pos = np.zeros([len(bin_bounds)]) #for tracking the current 'height' of stackplot, in order to staple and not overlap the values for each subprocess
			y0_neg = np.zeros([len(bin_bounds)]) #same as above for the negative contributions, plotted below x-axis

			color = iter(cm.rainbow(np.linspace(0,1,len(args['subproc']))))
			patches = []  #later needed for legend

			for i in range(0, len(args['subproc'])):	#go through fractions, subproc by subproc
				c = next(color)
				#positive contributions
				ax0.fill_between(bin_bounds.flatten(), steppify_bin(y0_pos), steppify_bin(y0_pos+pos_contr[i,:]), #label=labeling[args['subproc'][i]],
								facecolor=c, alpha=0.4)
				y0_pos += pos_contr[i, :] #calculate new height of stackplot

				#negative contributions
				ax0.fill_between(bin_bounds.flatten(), steppify_bin(y0_neg), steppify_bin(y0_neg+neg_contr[i, :]), #label=labeling[args['subproc'][i]],
								facecolor=c, hatch='X', alpha=0.4)
				y0_neg += neg_contr[i, :] #new 'height' below x-axis

				#patch for subprocess i (for the legend) ##labeling happens here
				patches.append(matplotlib.patches.Rectangle((0, 0), 0, 0, color=c,
								label=args['subproc'][i], alpha=0.4))
				ax0.add_patch(patches[i])
			
			#get limits for axes:
			xlim = ax0.get_xlim()
			ylim = ax0.get_ylim()
		
			print 'xlim:', xlim
			print 'ylim:', ylim, '\n'
		
			xmin_ = xlim[0]
			print 'xmin', xmin_
			xmax_ = xlim[1]
			#ymin = ylim[0]
			#ymax = ylim[1]
			ymin_ = np.amin(y0_neg)+0.3*np.amin(y0_neg)
			ymax_ = np.amax(y0_pos)+0.1*np.amax(y0_pos)
			
			#settings for the whole stackplot
			if args['constyaxis']==True:
				ax0.axis([0, xmax_, -0.20, 1.20])
			elif args['constyaxis']==False:
				ax0.set_xlim(xmin_, xmax_) ##
				ax0.set_ylim(ymin_, ymax_) ##set axis limit according to data
			
			#ax0.axis([30, 1000, -0.20, 1.20])
			ax0.set_xscale('log', nonposx='clip')
			ax0.autoscale(enable=True, axis='x', tight=True)

			##ax0.axis([0, xmax_, -0.20, 1.20]) #for test with fixed y-axis
			plt.axhline(y=1, xmin=0, xmax=1, color='k', linestyle='dotted') #dotted line at xs_sub/xs=1=100%
			plt.axhline(y=0, xmin=0, xmax=1, color='k', linestyle='-')		#line at xs_sub/xs_tot=0
			
			plt.title(x=0.5, y=1.06, s='%s' %(table_), fontsize=22)

			ax0.set_xlabel('x quantity', fontsize=20)
			ax0.xaxis.set_label_coords(0.9, -0.05)
			ax0.set_ylabel(r'$\frac{\mathrm{XS_{sp, %s}}}{\mathrm{XS_{tot, %s}}}$' %(sp_order, norm_order), rotation=0, fontsize=24)
			ax0.yaxis.set_label_coords(-0.16, 0.9)

			ax0.text(0.96, 0.03, args['pdfset']+', %s to %s' %(sp_order, norm_order), alpha=0.6, transform=ax0.transAxes, fontsize=14, verticalalignment='bottom', horizontalalignment='right')
			plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left') #'upper left' refers to the legend box

			plt.tight_layout()
			fig0.tight_layout()

			#naming of the plots
			if args['filename'] is not None:
				stackplotname = '%s_%s.png' %(args['filename'], sp_order+'_'+norm_order)
			else:
				stackplotname = '%s.stackplot_%s.png' %(tablename, sp_order+'_'+norm_order)
			fig0.savefig(stackplotname, bbox_inches='tight')

			print 'Stackplot for subprocess comparison saved as: %s' %(stackplotname)



		else:
			print "\n"
			print "----------------------------------------------------------"
			print "Create individual plot for each single subprocess \n"

			#individual plot for each subproc
			#create folder for saving the plots, will be overwritten in case it already exists (avoid countless folders while testing)
			current_dir = os.getcwd()
			
			if args['filename'] is not None:
				prename = args['filename']
			else: prename = tablename
			
			final_dir = os.path.join(current_dir, r'%s_sp_%s_%s' %(prename, sp_order, norm_order))

			if os.path.exists(final_dir):
				shutil.rmtree(final_dir) #remove final_dir if it already exists
			os.makedirs(final_dir)

			#loop where single plot for each of the selected subprocesses is created
			#start the plot
			color = iter(cm.rainbow(np.linspace(0,1,len(args['subproc']))))

			for i in range(0, len(args['subproc'])):
				c = next(color)
				plt.close()
				fig0 = plt.figure(figsize=(8, 7))
				ax0 = fig0.add_subplot(111)

				#plot the fractions
				y0 = np.zeros([len(bin_bounds)])
				
				#ax0.fill_between(bin_bounds.flatten(), steppify_bin(y0), steppify_bin(fractions_oneplot[i]), 
				#			facecolor=c, alpha=0.4)
				
				ax0.fill_between(bin_bounds.flatten(), steppify_bin(y0), steppify_bin(pos_contr[i]), 
							facecolor=c, alpha=0.4)
				ax0.fill_between(bin_bounds.flatten(), steppify_bin(y0), steppify_bin(neg_contr[i]), 
							facecolor=c, hatch='X', alpha=0.4)

				#get limits for axes:
				xlim = ax0.get_xlim()
				ylim = ax0.get_ylim()

				ax0.set_xscale('log', nonposx='clip')
				if args['constyaxis']==True:
					ax0.axis([xlim[0], xlim[1], -0.2, 1.2])
				elif args['constyaxis']==False:
					ax0.axis([xlim[0], xlim[1], ylim[0]+0.3*ylim[0], ylim[1]+0.1*ylim[1]])

				ax0.autoscale(enable=True, axis='x', tight=True)
				#ax0.axis([bin_bounds.flatten()[0], xlim[1], ylim[0]-0.1*ylim[0], ylim[1]+0.1*ylim[1]]) #leave y-limits fixed for better comparison between the subprocess-plots
				plt.axhline(y=1, xmin=0, xmax=1, color='k', linestyle='dotted')
				plt.axhline(y=0, xmin=0, xmax=1, color='k', linestyle='-')


				titlename = args['subproc'][i]

				plt.title(x=0.5, y=1.01, s=titlename, fontsize=20)
				ax0.set_xlabel('x quantity', fontsize=20)
				ax0.xaxis.set_label_coords(0.9, -0.05)
				ax0.set_ylabel(r'$\frac{\mathrm{XS_{sp, %s}}}{\mathrm{XS_{tot, %s}}}$' %(sp_order, norm_order), rotation=0, fontsize=24)
				ax0.yaxis.set_label_coords(-0.16, 0.9)
				ax0.text(0.96, 0.03, args['pdfset']+', %s to %s' %(sp_order, norm_order), transform=ax0.transAxes, fontsize=14, verticalalignment='bottom', horizontalalignment='right')
				ax0.text(0.96, 0.08, '%s' %table_, transform=ax0.transAxes, fontsize=14, verticalalignment='bottom', horizontalalignment='right')

				#plt.legend() #location parameter loc='upper right'

				plt.tight_layout()
				fig0.tight_layout()
				plotname = '%s_%s_singlefrac_%s_%s.png' %(tablename, args['subproc'][i], sp_order, norm_order)

				#save the plots
				fig0.savefig(os.path.join(final_dir,plotname))
				print 'Plot saved in %s' %(final_dir)
				print 'fraction plot saved as: %s' %(plotname), "\n \n"



	else:
		print "Something went wrong. Check choice of order(s). \n"







""" #other alternative possibility to switch orders on/off
		for order, setting in sorted(b_order.iteritems()): #go through tuples in b_order, evaluate
			print "\n", "order: ", order
			print "setting: ", setting
			index=0
			for j in setting:
				print "index", index
				print "j", j
				fnlo.SetContributionON(fastnlo.kFixedOrder,index,j);
				index+=1
"""












##################################################
##function to make better uncertainty plot (shaded)
##could maybe be replaced by .flatten() ??
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


#################################################

if __name__ == '__main__':
	main()
