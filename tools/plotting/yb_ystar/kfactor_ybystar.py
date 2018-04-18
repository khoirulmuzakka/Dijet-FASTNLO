#! /usr/bin/env python

######################################################################
#Looking at the k-factor
#(total cross section in NLO to total xs in LO)
#
#Comparison of multiple tables that have been chosen.
#No distinction between different subprocesses (always include ALL).
######################################################################

import argparse
import numpy as np
import os
import time, datetime
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.pyplot import cm

import fastnlo


def main():
	
	parser = argparse.ArgumentParser()
	
	#make it possible to use different options when running
	parser.add_argument('-t','--table', default=['table.tab'], required=True, type=str, nargs='*',
						help='FastNLO tables that shall be evaluated.')
	parser.add_argument('-p', '--pdfset', default='CT14nlo',
						help='PDFset(s) to evaluate fastNLO tables.')
	parser.add_argument('-m', '--member', default=0, type=int,
						help='Member of PDFset, default is 0.')
						


	#parse arguments
	args = vars(parser.parse_args())
	
	#table name
	table_names = []
	for i in range(0, len(args['table'])):
		table_ = os.path.basename(args['table'][i])
		tablename = os.path.splitext(table_)[0]
		tablename = os.path.splitext(tablename)[0] #because yb tables have .tab.gz ending (getting rid of both here)
		table_names.append(table_)
	tables = np.array(table_names)
	print '\n'
	print 'Tables: ', tables, '\n'
	
	#pdfset name
	pdfset_ = os.path.basename(args['pdfset'])
	pdfname = os.path.splitext(pdfset_)[0]
	print 'PDF Set: ', pdfname, '\n'
	
	
	#checks for investigated phase space region --> mainly for plot-labeling
	#range of y_boost
	def yb(table):
		if 'yb0' in table:
			yb_range = r'$\mathrm{0.0 \leq \/ y_b < 0.5}$'
			info_yb = r'0 <= yb < 0.5'
		elif 'yb1' in table:
			yb_range = r'$\mathrm{0.5 \leq \/ y_b < 1.0}$'
			info_yb = r'0.5 <= yb < 1.0'
		elif 'yb2' in table:
			yb_range = r'$\mathrm{1.0 \leq \/ y_b < 1.5}$'
			info_yb = r'1.0 <= yb < 1.5'
		elif 'yb3' in table:
			yb_range = r'$\mathrm{1.5 \leq \/ y_b < 2.0}$'
			info_yb = r'1.5 <= yb < 2.0'
		elif 'yb4' in table:
			yb_range = r'$\mathrm{2.0 \leq \/ y_b < 2.5}$'
			info_yb = r'2.0 <= yb < 2.5'
		else:
			print 'No info about range of yb.'
			yb_range = ''
			info_yb = ''
		return yb_range, info_yb		#yb_range is in correct format for labeling
			
	
	#range of ystar 
	def ystar(table):
		if 'ystar0' in table:
			ystar_range = r'$\mathrm{0.0 \leq \/ y^{\ast} < 0.5}$'
			info_ystar = r'0 <= ystar < 0.5'
		elif 'ystar1' in table:
			ystar_range = r'$\mathrm{0.5 \leq \/ y^{\ast} < 1.0}$'
			info_ystar = r'0.5 <= ystar < 1.0'
		elif 'ystar2' in table:
			ystar_range = r'$\mathrm{1.0 \leq \/ y^{\ast} < 1.5}$'
			info_ystar = r'1.0 <= ystar < 1.5'
		elif 'ystar3' in table:
			ystar_range = r'$\mathrm{1.5 \leq \/ y^{\ast} < 2.0}$'
			info_ystar = r'1.5 <= ystar < 2.0'
		elif 'ystar4' in table:
			ystar_range = r'$\mathrm{2.0 \leq \/ y^{\ast} < 2.5}$'
			info_ystar = r'2.0 <= ystar < 2.5'
		else:
			print 'No info about range of ystar.'
			ystar_range = ''
			info_ystar = ''
		return ystar_range, info_ystar

	

	
	######## EVALUATING ########

	fractions_list = []
	binbounds_list = [[] for i in range(len(args['table']))]
	#binbounds_array = np.array(
	for i in range(0, len(args['table'])):
		tab = tables[i]
		print 'Selected table: %s' %tab
		#show information about phase space region
		yb_range = yb(tab)[0]
		ystar_range = ystar(tab)[0]
		info_yb = yb(tab)[1]
		info_ystar = ystar(tab)[1]
		
		print 'Phase Space Region:'
		print '%s' %info_yb
		print '%s' %info_ystar, '\n'
		
		#Evaluate selected table:
		fnlo = fastnlo.fastNLOLHAPDF(args['table'][i], args['pdfset'], args['member'])
		
		####leading order cross section##################################################
		fnlo.SetContributionON(fastnlo.kFixedOrder,0,True);  #Switch LO on     
		fnlo.SetContributionON(fastnlo.kFixedOrder,1,False); #Switch NLO off
	
		#cross section with every process included (only Leading Order)
		fnlo.CalcCrossSection()
		xs_lo = np.array(fnlo.GetCrossSection())
		print '\n'
		print 'Cross section in Leading Order with all subprocesses xs_lo: \n'
		print xs_lo, '\n \n'

		###next-to-leading order cross section (LO+NLO):#################################
		fnlo.SetContributionON(fastnlo.kFixedOrder, 0, True); #switch LO on!
		fnlo.SetContributionON(fastnlo.kFixedOrder, 1, True); #switch NLO on
		
		#cross section with every process included (Next-to-Leading Order)
		fnlo.CalcCrossSection()
		xs_nlo = np.array(fnlo.GetCrossSection())
		print '\n'
		print 'Cross section in Next-to-Leading Order with all subprocesses xs_nlo: \n'
		print xs_nlo, '\n \n'

		#XS in NLO compared to XS in LO
		#alternative could be: first calculate all the nlo and lo XS for each table, then divide
		#here it is the other way round: first k-factors for each table, then put them into one list
		fraction = np.divide(xs_nlo,xs_lo) #nlo/lo for one table (in pTZ bins)
		fractions_list.append(fraction) #collect k-factor for all the tables in this list
		obs_binbounds = np.array(fnlo.GetObsBinsBounds(0)) #must be done for each table, as binbounds may differ
		#obs_binbounds = fnlo.GetObsBinBounds(0) #leave it as a list
		binbounds_list[i]= obs_binbounds
	fractions = np.array(fractions_list)
	binbounds = np.array(binbounds_list)
	print 'k-factors for the selected tables: \n'
	print 'fractions: \n', fractions, '\n'
	
	#print 'binbounds:', binbounds

	#NOTE: fractions and binbounds sometimes are arrays of lists that differ in length
	
	
	
	######## PLOTTING ########
	#taking care of the directory
	current_dir = os.getcwd()
	
	final_dir = os.path.join(current_dir, r'kfactor_plots_%s' %pdfname)
	now = datetime.datetime.now()
	nowstr = 'kfactor_plots_%s' %(pdfname)
	nowstr += now.strftime('_%Y-%m-%d_%H-%M-%S')
	final_dir_date = os.path.join(current_dir, nowstr)
	
	#check if final directory already exists, if so: include datetime to name
	if not os.path.exists(final_dir):
		os.makedirs(final_dir)
	else:
		os.makedirs(final_dir_date)
		
	
	plt.close()
	fig0 = plt.figure(figsize=(8,7))
	ax0 = fig0.add_subplot(111)
	
	
	color = iter(cm.rainbow(np.linspace(0,1,len(args['table']))))
	#print 'color:', color


	ymin = -0.20
	ymax = 1.20
	xmin = 30
	xmax = 1000
	for i in range(0, len(args['table'])):	#go through fractions, line by line = table by table
		bin_bounds = binbounds[i][:]
#		print 'x',bin_bounds
		fraction = fractions[i][:]
#		print 'y',fraction

		#plot the k-factors in each pTZ bin for table i:
		c = next(color) #use the iterator for choosing the linecolor fractions[i,:]
		ax0.plot(bin_bounds.flatten(), steppify_bin(fraction),
				 label=yb(args['table'][i])[0] + ",   " + ystar(args['table'][i])[0],
				 color=c, alpha=1.0)

		#check current limits of y-axis:
		ylim = ax0.get_ylim()
#		if (ylim[0]-0.05 < ymin):
#			ymin = ylim[0]-0.05
#		if (ylim[1]+0.05 > ymax):
#			ymax = ylim[1]+0.05

		#check current limits of x-axis:
		xlim = ax0.get_xlim()
#		if (xlim[0]-0.05 < xmin):
#			xmin = ylim[0]-0.05
#		if (xlim[1]+0.05 > xmax):
#			xmax = xlim[1]+0.05


	#settings for the whole plot
	ax0.set_xscale('log', nonposx='clip')
#	ax0.axis([xmin, xmax, ymin, ymax])
	ax0.axis([xlim[0], xlim[1], ylim[0]-0.05, ylim[1]+0.05]) #flexible axes
	plt.axhline(y=1, xmin=0, xmax=1, color='k', linestyle='dotted') 
	plt.axhline(y=0, xmin=0, xmax=1, color='k', linestyle='dotted')


	#plt.title(x=0.5, y=1.01, s='k-factor', fontsize=22)
	ax0.set_xlabel('$\mathrm{p_{T,Z}} \ \mathrm{[GeV]}$', fontsize=18)
	ax0.xaxis.set_label_coords(0.9, -0.05)
	ax0.set_ylabel(r'$\frac{\mathrm{XS_{tot,NLO}}}{\mathrm{XS_{tot,LO}}}$', rotation=0, fontsize=24)
	ax0.yaxis.set_label_coords(-0.14, 0.9)

	ax0.text(0.96, 0.03, args['pdfset'], transform=ax0.transAxes, fontsize=14, verticalalignment='bottom', horizontalalignment='right')
	#plt.legend(loc="best", fancybox=True, framealpha=0.5)
	plt.legend(bbox_to_anchor=(1.02,1), loc='upper left') #'upper left' refers to the legend box
	
	plt.tight_layout()
	fig0.tight_layout()

	figname = 'kfactors_%s.png' %(pdfname)
	#save the plot
	if os.path.exists(final_dir_date):
		fig0.savefig(os.path.join(final_dir_date,figname), bbox_inches='tight')
		print 'Plot saved in %s' %(final_dir_date)
	else:
		fig0.savefig(os.path.join(final_dir,figname), bbox_inches='tight')
		print 'Plot saved in %s' %(final_dir)
		print '\n'

	print 'k-factor plot saved as: %s' %(figname)
		






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





