#! /usr/bin/env python

####################################################
#
#Here, the contributions of subprocesses in LO are 
#compared to the total cross section in LO+NLO.
#
####################################################

import argparse
import numpy as np
import os
import time, datetime
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
						
	parser.add_argument('-s', '--subproc', default=['gg'], type=str, nargs='*',
						help='Subprocesses that shall be compared. \n'
							'Options: gg, gq, qiqi, qiai, qiqj, qiaj')
	parser.add_argument('-c', '--compare', default=False, const=True, type=bool, nargs='?',
						help='Stackplot if option -c is chosen. Otherwise single-process-plots.')
						#help='Stackplot if True. Single-Process-Plot if False.')
						


	#Dictionary for how the subprocesses are called
	processes = {"gg":"Gluon-Gluon",
				"gq":"Gluon-Quark/Antiquark",
				"qiqi":"Quark-Quark and $\mathrm{\overline{q_i} \; \overline{q_i}}$ (same flavor)",
				"qiai":"Quark-Antiquark (same flavor)",
				"qiqj":"Quark-Quark and $\mathrm{\overline{q_i} \; \overline{q_j}}}$ (different flavors)",
				"qiaj":"Quark-Antiquark (different flavors)"}

	#shorter for labeling the stackplot
	labeling = {"gg":r"$\mathrm{gg}$",
				"gq":r"$\mathrm{gq} \ & \ \mathrm{g\overline{q}}$",
				"qiqi":r"$\mathrm{q_i q_i} \ & \ \mathrm{\overline{q_i} \; \overline{q_i}}$",
				"qiai":r"$\mathrm{q_i} \; \mathrm{\overline{q_i}}$",
				"qiqj":r"$\mathrm{q_i q_j} \ & \ \mathrm{\overline{q_i} \; \overline{q_j}}$",
				"qiaj":r"$\mathrm{q_i} \; \mathrm{\overline{q_j}}$"}


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
	
	
	#checks for investigated phase space region --> mainly for plot-labeling
	
	#range of y_boost
	def yb():
		if 'yb0' in args['table']:
			yb_range = r'$\mathrm{0.0 \leq \/ y_b < 0.5}$'
			info_yb = r'0 <= yb < 0.5'
		elif 'yb1' in args['table']:
			yb_range = r'$\mathrm{0.5 \leq \/ y_b < 1.0}$'
			info_yb = r'0.5 <= yb < 1.0'
		elif 'yb2' in args['table']:
			yb_range = r'$\mathrm{1.0 \leq \/ y_b < 1.5}$'
			info_yb = r'1.0 <= yb < 1.5'
		elif 'yb3' in args['table']:
			yb_range = r'$\mathrm{1.5 \leq \/ y_b < 2.0}$'
			info_yb = r'1.5 <= yb < 2.0'
		elif 'yb4' in args['table']:
			yb_range = r'$\mathrm{2.0 \leq \/ y_b < 2.5}$'
			info_yb = r'2.0 <= yb < 2.5'
		else:
			print 'No info about range of yb.'
			yb_range = ''
			info_yb = ''
		return yb_range, info_yb		#yb_range is in correct format for labeling
	
	
	#range of ystar
	def ystar():
		if 'ystar0' in args['table']:
			ystar_range = r'$\mathrm{0.0 \leq \/ y^{\ast} < 0.5}$'
			info_ystar = r'0 <= ystar < 0.5'
		elif 'ystar1' in args['table']:
			ystar_range = r'$\mathrm{0.5 \leq \/ y^{\ast} < 1.0}$'
			info_ystar = r'0.5 <= ystar < 1.0'
		elif 'ystar2' in args['table']:
			ystar_range = r'$\mathrm{1.0 \leq \/ y^{\ast} < 1.5}$'
			info_ystar = r'1.0 <= ystar < 1.5'
		elif 'ystar3' in args['table']:
			ystar_range = r'$\mathrm{1.5 \leq \/ y^{\ast} < 2.0}$'
			info_ystar = r'1.5 <= ystar < 2.0'
		elif 'ystar4' in args['table']:
			ystar_range = r'$\mathrm{2.0 \leq \/ y^{\ast} < 2.5}$'
			info_ystar = r'2.0 <= ystar < 2.5'
		else:
			print 'No info about range of ystar.'
			ystar_range = ''
			info_ystar = ''
		return ystar_range, info_ystar



	#show information about phase space region
	yb_range = yb()[0]
	ystar_range = ystar()[0]
	
	info_yb = yb()[1]
	info_ystar = ystar()[1]

	print 'Phase Space Region:'
	print '%s' %info_yb
	print '%s' %info_ystar

	
	
	#Which subprocesses are evaluated?
	print 'Subprocesses that are investigated: \n'
	print args['subproc']
	
	
	######## EVALUATING ########
	fnlo = fastnlo.fastNLOLHAPDF(args['table'], args['pdfset'], args['member'])
	fnlo.SetContributionON(fastnlo.kFixedOrder,0,True);     
	fnlo.SetContributionON(fastnlo.kFixedOrder,1,True);         
	
	
	#cross section with every process included (for LO+NLO)
	fnlo.CalcCrossSection()
	xs_all = np.array(fnlo.GetCrossSection())
	print 'Cross section with all subprocesses xs_all: \n'
	print xs_all, '\n \n'
	
	#cross section of single processes in LO (only LO-contribution considered, no NLO!)
	fnlo.SetContributionON(fastnlo.kFixedOrder, 0, True); #switch LO on
	fnlo.SetContributionON(fastnlo.kFixedOrder, 1, False); #switch NLO off
	xs_subproc_list = []
	for i in range(0, len(args['subproc'])):
		subproc_name = args['subproc'][i]
		print 'Selected subprocess: %s' %subproc_name
		fnlo.SelectProcesses(subproc_name)
		fnlo.CalcCrossSection()
		xs_subproc = np.array(fnlo.GetCrossSection())
		xs_subproc_list.append(xs_subproc)
		print 'XS for subprocess %s: \n' %args['subproc'][i]
		print xs_subproc, '\n \n'
	xs_sub = np.array(xs_subproc_list) #Array with xs for selected processes in certain pTZ-bins
	print 'XS of the selected processes in different pTZ-regions: \n'
	print xs_sub, '\n \n'
	
	#how much do the selected subprocesses contribute to the total xs?
	fractions = np.divide(xs_sub,xs_all)
	print 'fractions: \n', fractions, '\n'
	
	for i in range(0, len(fractions)):
		print 'fraction of %s in LO on total XS: ' %args['subproc'][i]
		print fractions[i,:], '\n \n'
	
	
	#checking whether contribution is positive or negative
	pos_contr = np.array(fractions)
	neg_contr = np.array(fractions)
	bin_bounds = np.array(fnlo.GetObsBinsBounds(0))
	
	for i in range(0, len(args['subproc'])):	#go through 'fractions' line by line = subprocess after subprocess
		for j in range(0, len(bin_bounds)):					#check the different pTZ bins
			if (fractions[i,j] < 0):
				pos_contr[i,j] = 0.0
			else:
				neg_contr[i,j] = 0.0
		print 'Positive LO contributions of %s on total XS: \n' %args['subproc'][i]
		print pos_contr[i,:], '\n'
		print 'Negative LO contributions of %s on total XS: \n' %args['subproc'][i]
		print neg_contr[i,:], '\n \n'
	
	print 'Summary of positive and negative LO contributions of the subprocesses: \n'
	print pos_contr, '\n'
	print neg_contr, '\n'
	
	
	
	######## PLOTTING ########
	#Bins for x-axis (according to first pdf (0))
	bin_bounds = np.array(fnlo.GetObsBinsBounds(0))
	#print bin_bounds.T		#bin_bounds.T[0] = lower bounds, bin_bounds.T[1] upper bounds
	#print bin_bounds.flatten()
	
	
	#dictionary for using specific color per process
	process_color = {"gg":"g",
					"gq":"r",
					"qiqi":"y",
					"qiai":"c",
					"qiqj":"m",
					"qiaj":"b"}

	#Check whether comparison-plot or single-plot is required:
	print args['compare']
	if (args['compare']==True):		#stackplot that compares the contribution of the subprocesses to the total XS
		plt.close()
		fig0 = plt.figure(figsize=(8,7))
		ax0 = fig0.add_subplot(111)
		
		patches = []  #later needed for legend
		#number of bins via bin_bounds variable
		y0_pos = np.zeros([len(bin_bounds)]) #for tracking the current 'height' of stackplot, in order to staple and not overlap the values for each subprocess
		y0_neg = np.zeros([len(bin_bounds)]) #same as above for the negative contributions, plotted below x-axis
	
		for i in range(0, len(args['subproc'])):	#go through fractions, subproc by subproc
			
			#positive contributions
			ax0.fill_between(bin_bounds.flatten(), steppify_bin(y0_pos), steppify_bin(y0_pos+pos_contr[i]), #label=labeling[args['subproc'][i]],
							facecolor=process_color[args['subproc'][i]], alpha=0.4)
			y0_pos += pos_contr[i] #calculate new height of stackplot
							
			#negative contributions
			ax0.fill_between(bin_bounds.flatten(), steppify_bin(y0_neg), steppify_bin(y0_neg+neg_contr[i]), #label=labeling[args['subproc'][i]],
							facecolor=process_color[args['subproc'][i]], hatch='X', alpha=0.4)
			y0_neg += neg_contr[i] #new 'height' below x-axis

			#patch for subprocess i (for the legend) ##labeling happens here
			patches.append(matplotlib.patches.Rectangle((0, 0), 0, 0, color=process_color[args['subproc'][i]],
							label=labeling[args['subproc'][i]], alpha=0.4))
			ax0.add_patch(patches[i])
			
		#settings for the whole stackplot
		ax0.set_xscale('log', nonposx='clip')
		ax0.axis([30, 1000, -0.20, 1.20])
		plt.axhline(y=1, xmin=0, xmax=1, color='k', linestyle='dotted') #dotted line at xs_sub/xs=1=100%
		plt.axhline(y=0, xmin=0, xmax=1, color='k', linestyle='-')		#line at xs_sub/xs_tot=0
		
		#use phasespace region as title
		title_phasespace = yb_range + ",     " + ystar_range
		plt.title(x=0.5, y=1.01, s='%s' %(title_phasespace), fontsize=22)

		ax0.set_xlabel('$\mathrm{p_{T,Z}} \ \mathrm{[GeV]}$', fontsize=20)
		ax0.xaxis.set_label_coords(0.9, -0.05)
		ax0.set_ylabel(r'$\frac{\mathrm{XS_{sp,LO}}}{\mathrm{XS_{tot,LO+NLO}}}$', rotation=0, fontsize=22)
		ax0.yaxis.set_label_coords(-0.16, 0.9)
		ax0.text(0.96, 0.03, args['pdfset'], transform=ax0.transAxes, fontsize=14, verticalalignment='bottom', horizontalalignment='right')
		plt.legend()

		plt.tight_layout()
		fig0.tight_layout()
		
		#make string with the subprocess names 
		processnames = ""
		for i in range(0, len(args['subproc'])):
			processnames += ("_"+args['subproc'][i])
		stackplotname = 'compare_spLO_xstot%s_%s_%s.png' %(processnames, tablename, args['pdfset'])
		fig0.savefig(stackplotname)
    
		print 'Stackplot for subprocess comparison saved as: %s' %(stackplotname)
		
		
	
	
	else:					#individual plot for each subproc
		#create folder for saving the plots
		current_dir = os.getcwd()
		
		final_dir = os.path.join(current_dir, r'single_spLO_plots_%s' %tablename) 
		now = datetime.datetime.now()
		nowstr = 'single_spLO_xstot_plots_%s_' %(tablename)
		nowstr += now.strftime('%Y-%m-%d_%H-%M-%S')
		final_dir_date = os.path.join(current_dir, nowstr)
		
		#check if directory already exists, if so: include datetime to name
		if not os.path.exists(final_dir):
			os.makedirs(final_dir)
		else:
			os.makedirs(final_dir_date)
		
		#loop where single plot for each of the selected subprocesses is created
		#start the plot
		for i in range(0, len(args['subproc'])):
			plt.close()
			fig0 = plt.figure(figsize=(8, 7))
			ax0 = fig0.add_subplot(111)

			#plot the fractions
			y0 = np.zeros([len(bin_bounds)])
			ax0.fill_between(bin_bounds.flatten(), steppify_bin(y0), steppify_bin(fractions[i]), 
							facecolor=process_color[args['subproc'][i]], alpha=0.4) #improve labeling?
	
	
	
			#info about phase space region included as text (see below)
			ax0.set_xscale('log', nonposx='clip')
			ax0.axis([30, 1000, -0.4, 1.2])
			plt.axhline(y=1, xmin=0, xmax=1, color='k', linestyle='dotted')
			plt.axhline(y=0, xmin=0, xmax=1, color='k', linestyle='-')
	
			titlename = processes[args['subproc'][i]]

			plt.title(x=0.5, y=1.01, s=titlename, fontsize=20)
			ax0.set_xlabel('$\mathrm{p_{T,Z}} \ \mathrm{[GeV]}$', fontsize=20)
			ax0.xaxis.set_label_coords(0.9, -0.05)
			ax0.set_ylabel(r'$\frac{\mathrm{XS_{sp,LO}}}{\mathrm{XS_{tot,LO+NLO}}}$', rotation=0, fontsize=22)
			ax0.yaxis.set_label_coords(-0.16, 0.9)
			ax0.text(0.96, 0.03, args['pdfset'], transform=ax0.transAxes, fontsize=14, verticalalignment='bottom', horizontalalignment='right')
			ax0.text(0.96, 0.94, yb_range, transform=ax0.transAxes, fontsize=16, verticalalignment='bottom', horizontalalignment='right')
			ax0.text(0.96, 0.88, ystar_range, transform=ax0.transAxes, fontsize=16, verticalalignment='bottom', horizontalalignment='right')
			plt.legend() #location parameter loc='upper right'

			plt.tight_layout()
			fig0.tight_layout()
			plotname = '%s_LOfrac_xstot_%s_%s.png' %(args['subproc'][i], tablename, args['pdfset'])

		
		#save the plots
			if os.path.exists(final_dir_date):
				fig0.savefig(os.path.join(final_dir_date,plotname))
				print 'Plot saved in %s' %(final_dir_date)
			else:
				fig0.savefig(os.path.join(final_dir,plotname))
				print 'Plot saved in %s' %(final_dir)

			print '\n'
			print 'fraction plot saved as: %s' %(plotname)





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





