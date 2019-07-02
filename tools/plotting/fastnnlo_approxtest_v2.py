#!/usr/bin/env python2
#-*- coding:utf-8 -*-
#
# Make statistical evaluation plots of ensemble of one-to-one comparisons
# between fastNLO interpolation tables and NNLOJET original results
#
# Version:
#
# created by K. Rabbertz: 13.07.2017
# modified by B. Schillinger: 04.09.2018
#
#-----------------------------------------------------------------------
#
import argparse
import glob, os, pylab, sys
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator) #for xticks
import numpy as np
# from copy import deepcopy
from matplotlib import cm
from fastnlo import fastNLOLHAPDF
from fastnlo import SetGlobalVerbosity
import re
from StringIO import StringIO
#import warnings
#warnings.filterwarnings("error")

def main():
    parser = argparse.ArgumentParser()

    #add arguments
    parser.add_argument('-d','--datfile', default='file.dat', required=True, nargs='?',
                        help='.dat file for evaluation. All other matching .dat files that '
                                'only differ in seed no. will be evaluated as well.')
    parser.add_argument('-w','--weightfile', default='weight.txt', required=True, nargs='?',
                        help='.txt file containing weights.')
    parser.add_argument('-v', '--variation', action='store_true', #boolean value, True if -v chosen
                        help='If option is chosen: Plotting for all scale variations.'
                            'Otherwise (per default): Produce only central scale plots.'
                            'If specific scale is chosen via -f <int>, this option (-v) is ignored.')

    parser.add_argument('-f','--fscl', nargs='?', type=int, required=False,
                        help='Possible choices: int from fscl=1 to fscl=7.')

    parser.add_argument('-o', '--outputfilename', required=False, nargs='?', type=str,
                        help='Customise the first part of the output filename.'
                                'Default: Same structure as datfile name.')


    #parse arguments
    args = vars(parser.parse_args())
    namesp = parser.parse_args()
    print "Plots for all scale variations chosen? -->", args['variation']

    ################################################################################
    # Style settings
    # Location of matplotlibrc
    # print mpl.matplotlib_fname()
    # List active style settings
    # print mpl.rcParams
    # Needs matplotlib > 1.3
    # plt.style.use('ggplot')
    # plt.style.use('presentation')
    params = {'legend.fontsize': 'x-large',
              'figure.figsize': (16, 12),
              'axes.labelsize':  'x-large',
              'axes.titlesize':  'x-large',
              #'axes.linewidth':  2, #increase default linewidth
              'xtick.labelsize': 'x-large',
              'ytick.labelsize': 'x-large'}
    pylab.rcParams.update(params)


    #
    # Start
    #
    print "\n###########################################################################################"
    print "# fastnnlo_approxtest.py: Plot statistical evaluation of fastNLO interpolation quality"
    print "###########################################################################################\n"

    #
    # Default arguments             ##probably not anymore necessary
    #
    proc = '1jet'
    jobn = 'LO-CMS7'
    kinn = 'vBa'
    obsv = 'fnl2332d_xptji_y1'
    nmax = 99999
    fscl = 1

    #
    # Parse arguments
    #
    datfile = args['datfile']
    wgtfile = args['weightfile']
    fscl = args['fscl']
    print "fscl chosen by user: ", fscl
    #arguments
    print "Given datfile: ", datfile #, '\n'
    datbase = os.path.basename(datfile)
    datargs = datbase.split(".") #array containing filename-parts
    print "Given weightfile: ", wgtfile #, '\n'


    #arguments from datfile
    print 'arguments in datfile name: ', len(datargs)
    print "datargs: ", datargs, '\n'
    if len(datargs)==4:
        proc, jobn, obsv, ext = datargs
    elif len(datargs)==5:
        proc, jobn, kinn, obsv, ext = datargs
    elif len(datargs)==6:
        proc, jobn, kinn, obsv, seed, ext = datargs

    print 'proc: ', proc
    print 'jobn: ', jobn
    print 'kinn: ', kinn
    print 'obsv: ', obsv
    print 'seed: ', seed
    print '\n'

    # Extract order/contribution from job type (substring before first '-')
    order  = jobn.split('-')[0]

    # Prepare result arrays which are not used by Read_XS() function, but later
    xm = []      # bin "center"
    ntab = 0
    xs_fnlt = [] # fastNLO results
    nlin = 0
    weights = [] # Weight factors


    # Read binning and cross sections from NNLOJET dat file
    datglob = datfile[:-(len(seed)+len(ext))]+'*.dat'
    #print 'datglob ', datglob, '\n'

    datfiles = glob.glob(datglob)
    if not datfiles:
        print >> sys.stderr, 'No NNLOJET dat files matching', datglob ,'found, aborted!'
        sys.exit(1)
    datfiles.sort()

    xs_nnlo = None #will be updated in case only one scale is being looked at
                    #otherwise xs_nnlo_list will be used (and all scales investigated)
    xs_nnlo_list = None #see comment above (just the other way round)
    xs_fnll = None #same logic as above
    xs_fnll_list = None #same logic as above


    # New functionality! Check chosen options for scale setting.
    if fscl in np.arange(1,8,1): #if valid fscl is chosen via -f <int>
        xl, xu, ndat, xs_nnlo, nobs, seeds = Read_XS(datfiles, fscl)
    elif fscl is None: #option -f not used
        if (args['variation']==True): #evaluate all 7 scale variations (central + 6 others)
            xs_nnlo_list = [] #will become a list of arrays
            for fscl in np.arange(1,8,1):
                #take bin bounds from first iteration (fscl=1), they stay the same
                if fscl==1:
                    xl, xu, ndat, xs_nnlo_item, nobs, seeds = Read_XS(datfiles, fscl)
                    xs_nnlo_list.append(xs_nnlo_item)
                else:
                    xs_nnlo_item = Read_XS(datfiles, fscl)[3]
                    xs_nnlo_list.append(xs_nnlo_item)
            #xs_nnlo_list contains 7 arrays now
        elif (args['variation']==False): #evaluate central scale only
            fscl=1
            xl, xu, ndat, xs_nnlo, nobs, seeds = Read_XS(datfiles, fscl)
    else: #in case -f is chosen with invalid value
        sys.exit("No valid value for fscl chosen. Possible integers from fscl=1 up to fscl=7 .")

    if xs_nnlo is not None: print "shape of xs_nnlo: ", np.shape(xs_nnlo)
    if xs_nnlo_list is not None: print "shape of xs_nnlo_list", np.shape(xs_nnlo_list)


    # Read weights per file per bin from Alex
    wgtnams = np.genfromtxt(wgtfile, dtype=None, usecols=0)
    wgttmps = np.loadtxt(wgtfile, usecols=(list(range(1,nobs+1))))
    ntmp = len(wgtnams)
    # Combine to key-value tuple ( name, weights[] ) and sort according to key=name
    wgttup = zip(wgtnams,wgttmps)
    wgttup.sort(key = lambda row: (row[0]))
    # Unzip again
    allnames,allweights = zip(*wgttup)
    print 'datfile names in weightfile: \n', np.array(allnames), '\n'

    #print 'datfiles array: \n', np.array(datfiles), '\n'
    #print "weights..", weights
    for dfile in datfiles: #does not have to be changed, as choice of fscl makes no difference here
        newna = './'+order+'/'+os.path.basename(dfile)
        indexlin = allnames.index(newna)
        print 'Weight file line no. ', indexlin, ' is for ' , allnames[indexlin]
        weights.append(allweights[indexlin])
        nlin+=1
    weights = np.array(weights)

    print '\n', 'Using %s weight lines.' %nlin
    #print '\n','weights-array: \n', weights,'\n'


    # Evaluate cross sections from fastNLO tables
    # INFO=0, WARNING=1
    #SetGlobalVerbosity(1)
    #fnlotabs = glob.glob(pord+'/'+scen+'.'+proc+'.'+pord+'-'+ecms+'.???.'+obsv+'.*.tab.gz')
    #fnlotabs.sort()
    #for fnlotab in fnlotabs:
    #    print 'fastNLO table no. ', ntab, ' is ', fnlotab
    #    fnlo = fastNLOLHAPDF(fnlotab)
    #    fnlo.SetLHAPDFFilename('CT14nnlo')
    #    fnlo.SetLHAPDFMember(0)
    #    fnlo.SetMuFFunctionalForm(1); # kScale1=0
    #    fnlo.SetMuRFunctionalForm(1); # kScale2=1
    #    fnlo.CalcCrossSection()
    #    xs_fnla.append(fnlo.GetCrossSection())
    #    ntab += 1
    #    if ntab==nmax: break
    #print 'Using ', ntab, 'table files.'
    #xs_fnla = np.array(xs_fnla)

    # Evaluate cross sections from pre-evaluated fastNLO tables
    # assume that corresponding logfiles are located in the same directory as the datfiles
    log_list=[]

    for dfile in datfiles:
        log_list.append(dfile[:-(len(ext)+1)]+'.log')
    fnlologs = np.array(log_list)
    print 'fnlologs: \n', fnlologs, '\n'


    ############### adjusted to new logfiles (check if -v true or false) #########
    fnlologs.sort()
    if (xs_nnlo is not None): #only one certain scale choice fscl is considered
        xs_fnll, nlog = Read_logfile(fnlologs, fscl, seeds, nobs)
    elif (xs_nnlo_list is not None): #have list with values for all 7 fscl variations
        xs_fnll_list = []
        for fscl in range(1, 8, 1):
            xs_fnll_item, nlog = Read_logfile(fnlologs, fscl, seeds, nobs)
            xs_fnll_list.append(xs_fnll_item)

    #print 'xs_fnll', xs_fnll ####just to have a look at it

    # Check on identical file numbers, either for table or log files and weight lines
    if ndat == nlog == nlin and ndat*nlog*nlin != 0:
        print 'OK: Have equal number of dat files, fastNLO results, and weight lines: ndat = ', ndat, ', nlog = ', nlog, ', nlin = ', nlin
    else:
    #    sys.exit('ERROR: No matching file found or file number mismatch! ndat = '+str(ndat)+', ntab = '+str(ntab))
        sys.exit('ERROR: No matching file found or file number mismatch! ndat = '+str(ndat)+', nlog = '+str(nlog)+', nlin = '+str(nlin))

    # prepare info_values tuple for plotting function
    info_values = (weights, nobs, ndat, nlin, proc, jobn, order, obsv, args)

    # Use plotting function:
    print "\n"
    if xs_fnll is not None: #should be calculated already, in case only one fscl is investigated
        Statistics_And_Plotting(fscl, xs_fnll, xs_nnlo, info_values)
        #plot here and save the two plots
        print "Plotting done. (fscl=%s)" %fscl
    elif xs_fnll_list is not None: #that list exists if all fscl from 1 to 7 are considered, otherwise it is None
        r_f2nl_list=[]
        a_f2nl_list=[]
        for fscl in range(1,8,1):
            xs_fnll = xs_fnll_list[fscl-1]
            xs_nnlo = xs_nnlo_list[fscl-1]

            #plot here and save plots
            Statistics_And_Plotting(fscl, xs_fnll, xs_nnlo, info_values)
            print "Plotting done for fscl= %s" %fscl


    exit(0)
    #what are these plots like?
    for i_bin, (nnlo, fnlo, ratios, asyms, wgts) in enumerate(zip(xs_nnlo.T, xs_fnll.T, r_f2nl.T, a_f2nl.T, weights.T)):

        if i_bin > 0:
            exit(2222)
        print 'NNLOJET cross sections'
        bin_dists(nnlo, 'NNLOJET', r'$\sigma$ [pb]', 'blue', 'violet')
        print 'fastNLO cross sections'
        bin_dists(fnlo, 'fastNLO', r'$\sigma$ [pb]', 'blue', 'violet')
        print 'Ratios'
        bin_dists(ratios, 'Ratio', 'fastNLO/NNLOJET', 'blue', 'violet', True, rmed_f2nl[i_bin]-10.*riqd_f2nl[i_bin], rmed_f2nl[i_bin]+10.*riqd_f2nl[i_bin])
        print 'Asymmetries'
        bin_dists(asyms, 'Asymmetry', '(fastNLO-NNLOJET)/(fastNLO+NNLOJET)', 'blue', 'violet', True, amed_f2nl[i_bin]-10.*aiqd_f2nl[i_bin], amed_f2nl[i_bin]+10.*aiqd_f2nl[i_bin])

################################################################################################
##################################### END OF MAIN ##############################################
################################################################################################



### FUNCTIONS to be used in main() ###

################################################################################
# Define method for histogram statistics
################################################################################
def hist_stats(x, weights, axis=None):

# Fill proper np.array in case only 1. is given as weight
    try:
        iter(weights)
    except TypeError:
        weights = np.ones_like(x)*weights

# Unweighted mean and sample variance (not population variance)
    _ave = np.mean(x, axis=axis)
    if x.shape[axis] == 1:
        _var = np.zeros_like(_ave)
    else:
        _var = np.var(x, ddof=1, axis=axis)
    _std = np.sqrt(_var)

# Useful sums for weighted quantities
    sumw   = np.sum(weights, axis=axis)
    sumw2  = np.sum(weights * weights, axis=axis)
    sumwx  = np.sum(weights * x, axis=axis)
    sumwx2 = np.sum(weights * x * x, axis=axis)

# Weighted mean and reliability weighted sample variance
    _wave = []
    _wvar = []
    for i in range(len(sumw)):
        if sumw[i] == 0:
            _wave.append(np.nan)
        else:
            _wave.append(sumwx[i]/sumw[i])

    if x.shape[axis] == 1:
        _wvar = np.zeros_like(_ave)
    else:
        for i in range(len(sumw)):
            if sumw[i] == 0:
                _wvar.append(np.nan)
            else:
                _wvar.append((sumwx2[i] - sumwx[i]*sumwx[i]/sumw[i]) / (sumw[i] - sumw2[i]/sumw[i]))

    _wstd = np.sqrt(_wvar)

# Median and half interquartile distance
    _med = np.median(x, axis=axis)
    _med_err = np.subtract(*np.percentile(x, [75, 25], axis=axis, interpolation='linear'))/2.

    return dict(mean=_ave,
                stddev=_std,
                weighted_mean=_wave,
                weighted_stddev=_wstd,
                median=_med,
                iqd2=_med_err)
################################################################################
# End of method for histogram statistics
################################################################################

################################################################################
##### Function to read bin bounds and xs_nnlo from datfiles
################################################################################
def Read_XS(datfiles, fscl): #takes list of datfiles (different seeds) and fscl
    # Prepare result arrays
    xl = []      # left bin border
    xu = []      # right bin border
    ndat = 0
    xs_nnlo = [] # NNLOJET results
    seeds = []   # Seed numbers for matching
    # default maximum number of datfiles
    nmax = 99999

    print "Reading datfiles for fscl=%s ." %fscl
    ixscol = 3 + 2 * (fscl-1)
    for datfile in datfiles:
        print 'Datfile no. ', ndat, ' is ', datfile
        if ndat == 0:
            xl.append(np.loadtxt(datfile,usecols=(0,)))
            xu.append(np.loadtxt(datfile,usecols=(2,)))
        xs_nnlo.append(np.loadtxt(datfile,usecols=(ixscol,)))
        #print "xs_nnlo \n", xs_nnlo ## just to see what's happening
        parts = datfile.split(".")
        seed  = parts[len(parts)-2]
        seeds.append(seed)
        ndat += 1
        if ndat == nmax:
            break
    print 'Using ', ndat, 'dat files.'
    xl = np.array(xl)
    xu = np.array(xu)
    xs_nnlo = np.array(xs_nnlo)/1000. # Conversion of fb to pb

    # Determine no. of observable bins
    nobs = xl.size
    print 'Number of observable bins: ', nobs, '\n'
    return xl, xu, ndat, xs_nnlo, nobs, seeds
################### End of function ###########################################
###############################################################################

###############################################################################
############### Function to read XS from logfile
###############################################################################
def Read_logfile(fnlologs, fscl, seeds, nobs): #takes list of logfiles and fscl,
                                                #as well as list of seeds and nobs (number of observable bins)
    nlog = 0
    # default maximum number of logfiles
    nmax = 99999
    xs_fnll = []
    for fnlolog in fnlologs:
        print 'fastNLO file no. ', nlog, ' is ', fnlolog
        if not seeds[nlog] in fnlolog:
            print 'Mismatch in result sort order between NNLOJET and fastNLO. Aborted!'
            print 'seeds[',nlog,'] = ', seeds[nlog],', fastNLO file = ', fnlolog
        else:
            print 'NNLOJET and fastNLO result correctly matched. seed is ', seeds[nlog]

        xs_tmp = np.loadtxt(fnlolog,usecols=(6,),comments=['#',' #','C','L','N'])
        #print "xs_tmp \n", xs_tmp

        indi = (fscl-1)*nobs #skip lines of lower fscl
        indf = indi + nobs #last line in logfile which belongs to certain fscl
        xs_sub = xs_tmp[indi:indf]
        xs_fnll.append(xs_sub)
        nlog += 1
        if nlog==nmax: break
    print 'Using ', nlog, 'pre-evaluated table files.'
    xs_fnll = np.array(xs_fnll)
    return xs_fnll, nlog
############### End of function to read logfile ###############################
###############################################################################

###############################################################################
################ Function for Statistics and Plotting
###############################################################################
def Statistics_And_Plotting(fscl, xs_fnll, xs_nnlo, info_values): #takes pre-calculated values for specific fscl (process variables as tuple)
    weights, nobs, ndat, nlin, proc, jobn, order, obsv, args = info_values #get required variables from tuple

    r_f2nl = xs_fnll/xs_nnlo
    a_f2nl = (xs_fnll - xs_nnlo)/(xs_fnll + xs_nnlo)

    #Ratio statistics
    _ratio_stats = hist_stats(r_f2nl, weights=weights, axis=0)
    rave_f2nl    = _ratio_stats['mean']
    rstd_f2nl    = _ratio_stats['stddev']
    rwav_f2nl    = _ratio_stats['weighted_mean']
    rwst_f2nl    = _ratio_stats['weighted_stddev']
    rmed_f2nl    = _ratio_stats['median']
    riqd_f2nl    = _ratio_stats['iqd2']

    # Asymmetry statistics
    _asymm_stats = hist_stats(a_f2nl, weights=weights, axis=0)
    aave_f2nl = _asymm_stats['mean']
    astd_f2nl = _asymm_stats['stddev']
    awav_f2nl = _asymm_stats['weighted_mean']
    awst_f2nl = _asymm_stats['weighted_stddev']
    amed_f2nl = _asymm_stats['median']
    aiqd_f2nl = _asymm_stats['iqd2']

    #Settings for x-axis of plots
    #For setting small xticks for each bin with only certain labels
    #Set minimum > unity in MultipleLocator
    majorLocator = MultipleLocator(max((nobs//6),1.1))
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(1) #set minor tick for each observable bin


    # Ratio plot
    limfs = 'x-large'
    fig1 = plt.figure(figsize=(16,12)) #renamed fig to fig1, so function can return it
    ax1  = fig1.gca() #
    plt.title(r'Ratio: {} {} for {} using scale {}'.format(proc, order, obsv, fscl), fontsize='x-large', fontweight='bold', loc='left')
    plt.xlabel('Observable bin index', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
    plt.ylabel('fastNLO/NNLOJET', horizontalalignment='right', x=1.0, verticalalignment='bottom', y=1.0)
    plt.axhline(y=1.001, linestyle='--', linewidth=1.0, color='black')
    plt.axhline(y=0.999, linestyle='--', linewidth=1.0, color='black')
    plt.fill_between([0.0,34.0],0.999,1.001, color='black', alpha=0.1)
    #plt.text(33.6,1.00085,u'+1‰')
    #plt.text(33.7,0.99885,u'–1‰')
    plt.text(nobs+1.1,+1.00085,u'+1‰')
    plt.text(nobs+1.1,+0.99885,u'–1‰') #location flexible adjusted to plotting range

    dx = 1./3.
    xa = np.arange(1-dx, nobs-dx+1.e-6)
    xb = np.arange(1   , nobs   +1.e-6)
    xc = np.arange(1+dx, nobs+dx+1.e-6)

    rmean = plt.errorbar(xa, rave_f2nl, yerr=0.*rstd_f2nl, marker='v', linestyle='none', label=r'$\mu$ $\pm$ $\Delta\mu$', color='blue')
    plt.errorbar(xa, rave_f2nl, yerr=rstd_f2nl/np.sqrt(ndat), marker='.', linestyle='none', color='blue')
    if nlin>0:
        rmwgt = plt.errorbar(xb, rwav_f2nl, yerr=0.*rwst_f2nl, marker='^', linestyle='none', label=r'$\mu_w$ $\pm$ $\Delta\mu_w$', color='violet')
        plt.errorbar(xb, rwav_f2nl, yerr=rwst_f2nl/np.sqrt(nlin), marker='.', linestyle='none', color='violet')
    rmedian = plt.errorbar(xc, rmed_f2nl, yerr=riqd_f2nl, marker='s', linestyle='none', label=r'median $\pm$ IQD/2', color='red')

    #plt.xlim(0.0,34.0)
    plt.xlim(0, nobs+1)
    #plt.xticks(range(0, nobs+2, nobs//6))
    ax1.xaxis.set_major_locator(majorLocator)
    ax1.xaxis.set_major_formatter(majorFormatter)
    ax1.xaxis.set_minor_locator(minorLocator)
    ax1.tick_params(which='both', width=1)
    ax1.tick_params(which='major', length=7)
    ax1.tick_params(which='minor', length=5)
    plt.ylim(0.99,1.01)
    #plt.ylim(-50, 12000) #to display points for fscl=0 (log6file), was just a test
    if nlin>0:
        handles = [rmean,rmwgt,rmedian]
    else:
        handles = [rmean,rmedian]
    labels  = [h.get_label() for h in handles]
    legend = ax1.legend(handles, labels, title=r'No. of entries/bin = {}'.format(ndat), loc='upper left', numpoints=1, handlelength=0)
    legend.get_title().set_fontsize(limfs)

    if args['outputfilename'] is None:
        fignam = proc+'.'+jobn+'.'+obsv+'.'+'fscl'+str(fscl)+'.'+'approx_ratio_'+str(ndat)+'.png'
    else:
        fignam = args['outputfilename']+'.fscl'+str(fscl)+'.approx_ratio_'+str(ndat)+'.png'
    plt.savefig(fignam)
    #print "Ratio plot saved as: %s" %fignam
    #plt.show()

    # Asymmetry plot
    limfs = 'x-large'
    fig2 = plt.figure(figsize=(16,12)) #renamed fig to fig2
    ax2  = fig2.gca() #
    plt.title(r'Asymmetry: {} {} for {} using scale {}'.format(proc, order, obsv, fscl), fontsize='x-large', fontweight='bold', loc='left')
    plt.xlabel('Observable bin index', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
    plt.ylabel('(fastNLO-NNLOJET)/(fastNLO+NNLOJET)', horizontalalignment='right', x=1.0, verticalalignment='bottom', y=1.0)
    plt.axhline(y=-0.001, linestyle='--', linewidth=1.0, color='black')
    plt.axhline(y=+0.001, linestyle='--', linewidth=1.0, color='black')
    plt.fill_between([0.0,34.0],-0.001,0.001, color='black', alpha=0.1)
    #plt.text(33.6,+0.00085,u'+1‰')
    #plt.text(33.7,-0.00115,u'–1‰')
    plt.text(nobs+1.1,+0.00085,u'+1‰')
    plt.text(nobs+1.1,-0.00115,u'–1‰') #location flexible adjusted to plotting range

    dx = 1./3.
    xa = np.arange(1-dx, nobs-dx+1.e-6)
    xb = np.arange(1   , nobs   +1.e-6)
    xc = np.arange(1+dx, nobs+dx+1.e-6)

    amean = plt.errorbar(xa, aave_f2nl, yerr=0.*astd_f2nl, marker='<', linestyle='none', label='$\mu$ $\pm$ $\Delta\mu$', color='orange')
    plt.errorbar(xa, aave_f2nl, yerr=astd_f2nl/np.sqrt(ndat), marker='.', linestyle='none', color='orange')
    if nlin>0:
        amwgt = plt.errorbar(xb, awav_f2nl, yerr=0.*awst_f2nl, marker='>', linestyle='none', label=r'$\mu_w$ $\pm$ $\Delta\mu_w$', color='brown')
        plt.errorbar(xb, awav_f2nl, yerr=awst_f2nl/np.sqrt(ndat), marker='.', linestyle='none', color='brown')
    amedian = plt.errorbar(xc, amed_f2nl, yerr=aiqd_f2nl, marker='8', linestyle='none', label=r'median $\pm$ IQD/2', color='green')

    #plt.xlim(0.0,34.0)
    plt.xlim(0, nobs+1)
    ax2.xaxis.set_major_locator(majorLocator)
    ax2.xaxis.set_major_formatter(majorFormatter)
    ax2.xaxis.set_minor_locator(minorLocator)
    ax2.tick_params(which='both', width=1)
    ax2.tick_params(which='major', length=7)
    ax2.tick_params(which='minor', length=5)
    plt.ylim(-0.01,0.01)
    #plt.ylim(-1.02, 1.02) ## to display points for fscl=0 (log6file), just a test
    if nlin>0:
        handles = [amean,amwgt,amedian]
    else:
        handles = [amean,amedian]
    labels  = [h.get_label() for h in handles]
    legend = ax2.legend(handles, labels, title=r'No. of entries/bin = {}'.format(ndat), loc='upper left', numpoints=1, handlelength=0)
    legend.get_title().set_fontsize(limfs)

    if args['outputfilename'] is None:
        fignam = proc+'.'+jobn+'.'+obsv+'.'+'fscl'+str(fscl)+'.'+'approx_asymm_'+str(ndat)+'.png'
    else:
        fignam = args['outputfilename']+'.fscl'+str(fscl)+'.approx_asymm_'+str(ndat)+'.png' #should fscl necessarily be added, too?
    plt.savefig(fignam)
    #print "Asymmetry plot saved as: %s" %fignam
    #plt.show()

    #produces 1 plot for asymmetry and one for ratio (so it is unnecessary to return all intern variables)
    return fig1, fig2
################### End of function for statistics & plotting ################
##############################################################################





#the following method was not used in default settings of original script

################################################################################
# Define method for binwise distribution plots
################################################################################
def bin_dists(q, label, xlabel, col1, col2, indicators=False, lval=0., rval=0. ):
    _min, _max = np.min(q), np.max(q)
    _range = _max - _min
    _nbins = min(int(ndat/5), 50)
    _bins = np.linspace(_min, _max, _nbins + 1)

    fig_dist = plt.figure(figsize=(12,12))

    _h1, _c1, _ = plt.hist(q, bins=_bins, color=col1)
    _h2, _c2, _ = plt.hist(q, weights=wgts/np.sum(wgts)*ndat, bins=_bins, color=col2)

    plt.gcf().clf()
    plt.bar(_bins[:-1], _h1, color=col1, width=0.5*_range/_nbins)
    plt.bar(_bins[:-1]+0.5*_range/_nbins, _h2, color=col2, width=0.5*_range/_nbins)

    plt.title(r'{}: {} {} for {}'.format(label, proc, order, obsv), fontsize='x-large', fontweight='bold', loc='left')
    plt.xlabel(xlabel, horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
    plt.ylabel('Frequency', horizontalalignment='right', x=1.0, verticalalignment='bottom', y=1.0)
    raw = mpatches.Patch(color=col1, label='raw')
    wgt = mpatches.Patch(color=col2, label='weighted')
    if indicators:
        vl1 = plt.axvline(x=lval, linestyle='--', linewidth=2.0, color='red', label=r'$\pm10$ IQD/2')
        vl2 = plt.axvline(x=rval, linestyle='--', linewidth=2.0, color='red')
        handles = [vl1, raw, wgt]
    else:
        handles = [raw, wgt]
    labels  = [h.get_label() for h in handles]
#    plt.fill_between([-0.5,33.5],0.999,1.001, color='black', alpha=0.1)
#    plt.text(33.6,1.00085,u'+1‰')
#    plt.text(33.7,0.99885,u'–1‰')

    plt.legend(handles, labels, title=r'Observable bin no. {}'.format(i_bin), loc='upper left', numpoints=1)

    plt.yscale('log')
    plt.ylim((0.1, None))

    fig_dist_name = "{}.{}.{}.approx_{:d}_{}_bin{:d}.png".format(proc, jobn, obsv, ndat, label, i_bin)
    fig_dist.savefig(fig_dist_name)
    plt.close(fig_dist)

###################################################################################




if __name__=="__main__":
    main()
