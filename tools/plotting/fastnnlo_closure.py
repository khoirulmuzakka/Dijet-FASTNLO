#!/usr/bin/env python2
#-*- coding:utf-8 -*-
#
# Make statistical evaluation plots of ensemble of one-to-one comparisons
# between fastNLO interpolation tables and original NNLOJET results
#
# Version:
#
# created by K. Rabbertz: 13.07.2017
# modified by B. Schillinger: 04.09.2018
# modified & renamed by K. Rabbertz: 09.04.2019
#
#-----------------------------------------------------------------------
#
# Use matplotlib with Cairo offline backend for png, eps, or svg output
import matplotlib as mpl
mpl.use('Cairo')
import argparse, glob, os, pylab, re, sys
from StringIO import StringIO
import matplotlib.lines as mpllines
import matplotlib.gridspec as gridspec
import matplotlib.patches as mplpatches
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import (FormatStrFormatter, ScalarFormatter, AutoMinorLocator, MultipleLocator)
from matplotlib import cm
# numpy
import numpy as np

# Redefine ScalarFormatter
class ScalarFormatterForceFormat(ScalarFormatter):
    def _set_format(self,vmin,vmax):  # Override function that finds format to use.
        self.format = "%1.2f"  # Give format here

# Action class to allow comma-separated list in options
class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

# Some global definitions
_formats        = {'eps':0, 'png':1, 'svg':2}
_debug          = False

########################################################################################################################

def main():
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

    # Define arguments & options
    parser = argparse.ArgumentParser(epilog='',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Positional arguments
    parser.add_argument('datfile', type=str,
                        help='NNLOJET result (with .dat extension) for closure test. All other matching .dat files that only differ in seed no. will be evaluated as well up to a maximum of nmax files in total, see below.')
    # Optional arguments
    parser.add_argument('-d', '--debug', required=False, default=False, type=bool,
                        help='Switch on debug/verbose mode.')
    parser.add_argument('-f', '--filename', required=False, default=None, type=str,
                        help='Customise the basename for output filenames instead of datfile name.')
    parser.add_argument('--format', required=False, nargs='?', type=str, action=SplitArgs,
                        help='Comma-separated list of plot formats to use: eps, png, svg. If nothing is chosen, png is used.')
    parser.add_argument('-n', '--nmax', default=1, type=int, choices=range(1,10000), metavar='[1-9999]',
                        help='Number of NNLOJET results to compare one-by-one to fastNLO.')
    parser.add_argument('-s','--scale', default=0, type=int, choices=range(0,8),
                        help='Preset NNLOJET scale setting to compare. By default (0) compare all seven, else choose 1 to 7.')
    parser.add_argument('-w','--weightfile', required=False, default=None, type=str,
                        help='NNLOJET combination APPLfast.txt file containing merging weights.')

    # Parse arguments
    args = vars(parser.parse_args())

    #
    # Start
    #
    print ''
    print '###########################################################################################'
    print '# fastnnlo_closure.py: Plot statistical evaluation of fastNLO interpolation quality'
    print '###########################################################################################'
    print ''

    #
    # Parse arguments
    #
    # NNLOJET dat file
    datfile = args['datfile']
    print '[fastnnlo_closure]: Given datfile: ', datfile
    # Get details from datfile basename
    datbase  = os.path.basename(datfile)
    datargs  = datbase.split(".")
    ndatargs = len(datargs)
    print '[fastnnlo_closure]: No. of arguments derived from datfile name: ', ndatargs
    kinn = ''
    if ndatargs == 4:
        proc, jobn, obsv, ext = datargs
    elif ndatargs == 5:
        proc, jobn, kinn, obsv, ext = datargs
    elif ndatargs == 6:
        proc, jobn, kinn, obsv, seed, ext = datargs
    else:
        print >> sys.stderr, '[fastnnlo_closure]: ERROR! Incompatible datfile name: ', datfile
        sys.exit(1)
    if ext != 'dat':
        print >> sys.stderr, '[fastnnlo_closure]: ERROR! Incompatible datfile extension: ', ext
        sys.exit(1)
    # Extract order/contribution from job type (substring before first '-')
    order  = jobn.split('-')[0]
    print '[fastnnlo_closure]: NNLOJET process: ', proc
    print '[fastnnlo_closure]: NNLOJET job name: ', jobn
    print '[fastnnlo_closure]: NNLOJET order/channel: ', order
    print '[fastnnlo_closure]: NNLOJET kinematics (if any): ', kinn
    print '[fastnnlo_closure]: Observable: ', obsv
    print '[fastnnlo_closure]: Seed (if any): ', seed

    # Debug mode
    _debug = args['debug']
    print '[fastnnlo_closure]: Debug mode is set to: ', _debug

    # Given output filename base
    if args['filename']:
        print '[fastnnlo_closure]: Desired output filename base: ', args['filename']

    # Plot formats to use
    formats = args['format']
    if formats is None: formats = ['png']
    for fmt in formats:
        if not _formats.has_key(fmt):
            print '[fastnnlo_closure]: Illegal format specified, aborted!'
            print '[fastnnlo_closure]: Format list: ', args['format']
            exit(1)

    # No. of dat files to use for statistical evaluation
    nmax = args['nmax']
    print '[fastnnlo_closure]: No. of dat files to check for closure: ', nmax

    # Scale choice
    fscl = args['scale']
    print '[fastnnlo_closure]: Scale settings to check for closure: ', nmax

    # Weight file
    wgtfile = args['weightfile']
    if wgtfile:
        print '[fastnnlo_closure]: Weight file to use: ', wgtfile



    # Prepare result arrays which are not used by Read_XS() function, but later
    xm = []      # bin "center"
    ntab = 0
    xs_fnlt = [] # fastNLO results



    # Read binning and cross sections from NNLOJET dat files
    datglob = datfile[:-(len(seed)+len(ext))]+'*.dat'
    if _debug: print '[fastnnlo_closure]: datfile glob: ', datglob

    datfiles = glob.glob(datglob)
    if not datfiles:
        print >> sys.stderr, '[fastnnlo_closure]: ERROR! No NNLOJET dat files found matching ', datglob ,', aborted!'
        sys.exit(1)
    datfiles.sort()
    xl, xu, seeds, xs_ref = Read_XS_from_datfiles(datfiles, nmax)

    # Determine no. of observable bins
    nobs = xl.size
    print '[fastnnlo_closure]: Number of observable bins: ', nobs

    # Print no. of dat files
    ndat = seeds.size
    print '[fastnnlo_closure]: Using ', ndat, 'dat files.'

    # Read weights per file per bin from Alex' APPLfast.txt file
    wgtnams = np.genfromtxt(wgtfile, dtype=None, usecols=0)
    wgttmps = np.loadtxt(wgtfile, usecols=(list(range(1,nobs+1))))
    ntmp = len(wgtnams)
    # Combine to key-value tuple ( name, weights[] ) and sort according to key=name
    wgttup = zip(wgtnams,wgttmps)
    wgttup.sort(key = lambda row: (row[0]))
    # Unzip again
    allnames, allweights = zip(*wgttup)
    # Loop over datfiles up to nmax
    nlin = 0
    weights = [] # Weight factors
    for dfile in datfiles:
        # Adapt datfile name to match APPLfast.txt entry
        afname = './'+order+'/'+os.path.basename(dfile)
        ilin = allnames.index(afname)
        weights.append(allweights[ilin])
        nlin+=1
        if nlin == nmax: # Stop with statistics if nmax reached
            break
    weights = np.array(weights)

    # Print no. of weight lines
    print '[fastnnlo_closure]: Using ', nlin, 'weight lines.'

    # Read cross sections from pre-evaluated fastNLO tables
    # Assume that corresponding logfiles are located in the same directory as the datfiles
    logfiles = []
    for dfile in datfiles:
        logfiles.append(dfile[:-(len(ext)+1)]+'.log')
    logfiles = np.array(logfiles)
    logfiles.sort()
    nlog, xs_grid = Read_XS_from_logfiles(logfiles, nmax, nobs, seeds)

    # Print no. of log files
    print '[fastnnlo_closure]: Using ', nlog, 'pre-evaluated table files.'

    # Check on identical file numbers for log files and weight lines
    if ndat == nlog == nlin and ndat*nlog*nlin != 0:
        if _debug: print '[fastnnlo_closure]: OK: Have equal number of dat files, fastNLO results, and weight lines: ndat = ', ndat, ', nlog = ', nlog, ', nlin = ', nlin
    else:
        sys.exit('[fastnnlo_closure]: ERROR: No matching file found or file number mismatch! ndat = '+str(ndat)+', nlog = '+str(nlog)+', nlin = '+str(nlin))

    # Prepare info_values tuple for plotting function
    info_values = (weights, nobs, ndat, nlin, proc, jobn, order, obsv, args)



    # DEBUG
    print 'xl',xl,xl.shape
    print 'xu',xu,xu.shape
    print 'seeds',seeds,seeds.shape
    print 'xs_ref',xs_ref,xs_ref.shape
    print 'xs_grid',xs_grid,xs_grid.shape



    # Use plotting function:
#    if xs_fnll is not None: #should be calculated already, in case only one fscl is investigated
#        Statistics_And_Plotting(fscl, xs_fnll, xs_nnlo, info_values)
#        #plot here and save the two plots
#        print "Plotting done. (fscl=%s)" %fscl
#    elif xs_fnll_list is not None: #that list exists if all fscl from 1 to 7 are considered, otherwise it is None
    r_f2nl_list=[]
    a_f2nl_list=[]
    for fscl in range(1,8,1):
        xs_fnll = xs_grid[:,:,fscl-1]
        xs_nnlo = xs_ref[:,:,fscl-1]

        print 'grid shape', xs_grid.shape
        print 'ref shape', xs_ref.shape
        print 'fnll shape',xs_fnll.shape
        print 'nnlo shape',xs_nnlo.shape

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
##### Function to read bin bounds and cross section from datfiles
################################################################################
# Takes list of all matching datfiles (different seeds) and the no. of datfiles to use
# Always reads all seven scale-varied results to avoid accessing files multiple times
def Read_XS_from_datfiles(datfiles, nmax):
    # Prepare result arrays
    xl = []   # Lower bin border
    xu = []   # Upper bin border
    xs = []   # NNLOJET cross sections
    si = []   # 5-digit seed number index (with leading 's') for matching

    ndat = 0
    print "[fastnnlo_closure]: Reading datfiles ..."
    for datfile in datfiles:
        if _debug: print '[fastnnlo_closure]: Reading datfile no. ', ndat, ': ', datfile
        data = np.loadtxt(datfile,usecols=(0,2,3,5,7,9,11,13,15))
        if ndat == 0: # Use only once, since bin borders stay the same
            xl.append(data[:,0])
            xu.append(data[:,1])
        xs.append(data[:,2:9]) # Read seven scale-varied cross sections
        parts = datfile.split(".")
        seed  = parts[len(parts)-2]
        si.append(seed) # Store seed number from datfile name
        ndat += 1
        if ndat == nmax: # Stop with statistics if nmax reached
            break
    xl = np.array(xl)
    xu = np.array(xu)
    si = np.array(si)
    # TODO: Check from interpolation grid|log file whether this conversion is necessary!!!
    xs = np.array(xs)/1000. # Conversion of fb to pb
    return xl, xu, si, xs
################### End of function ###########################################
###############################################################################

###############################################################################
############### Function to read cross section from logfiles
###############################################################################
# Takes list of all matching logfiles (different seeds) and the no. of logfiles to use
# Always reads all seven nobs scale-varied results to avoid accessing files multiple times
# Crosschecks against seed number indices from datfiles
def Read_XS_from_logfiles(logfiles, nmax, nobs, si):
    # Prepare result array
    xs = []   # Interpolated cross sections

    nlog = 0
    print "[fastnnlo_closure]: Reading logfiles ..."
    for logfile in logfiles:
        xs_log = []
        if _debug: print '[fastnnlo_closure]: Reading logfile no. ', nlog, ': ', logfile
        if not si[nlog] in logfile:
            print >> sys.stderr, '[fastnnlo_closure]: ERROR! Mismatch in result sort order between NNLOJET and fastNLO. Aborted!'
            print >> sys.stderr, '[fastnnlo_closure]: si[',nlog,'] = ', si[nlog],', fastNLO file = ', logfile
            sys.exit(1)
        else:
            if _debug: print '[fastnnlo_closure]: NNLOJET and fastNLO result correctly matched. seed is ', si[nlog]

        xs_tmp = np.loadtxt(logfile,usecols=(6,),comments=['#',' #','C','L'])
        print "xs_tmp \n", xs_tmp

        # Subdivide list of all seven scale-varied results into nobs pieces
        for iobs in range(nobs):
            xobs = [xs_tmp[iobs+0*nobs], xs_tmp[iobs+1*nobs], xs_tmp[iobs+2*nobs], xs_tmp[iobs+3*nobs], xs_tmp[iobs+4*nobs], xs_tmp[iobs+5*nobs], xs_tmp[iobs+6*nobs]]
            xs_log.append(xobs)

# Subdivide list of all scale-varied results into seven pieces
#        for fscl in range(1,8):
#            indi = (fscl-1)*nobs # Skip lines of lower fscl
#            indf = indi + nobs   # Last line in logfile which belongs to certain fscl
#            xs.append(xs_tmp[indi:indf])

        xs.append(xs_log)
        nlog += 1
        if nlog == nmax: # Stop with statistics if nmax reached
            break
    xs = np.array(xs)
    return nlog, xs
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

    if args['filename'] is None:
        fignam = proc+'.'+jobn+'.'+obsv+'.'+'fscl'+str(fscl)+'.'+'approx_ratio_'+str(ndat)+'.png'
    else:
        fignam = args['filename']+'.fscl'+str(fscl)+'.approx_ratio_'+str(ndat)+'.png'
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

    if args['filename'] is None:
        fignam = proc+'.'+jobn+'.'+obsv+'.'+'fscl'+str(fscl)+'.'+'approx_asymm_'+str(ndat)+'.png'
    else:
        fignam = args['filename']+'.fscl'+str(fscl)+'.approx_asymm_'+str(ndat)+'.png' #should fscl necessarily be added, too?
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
    raw = mplpatches.Patch(color=col1, label='raw')
    wgt = mplpatches.Patch(color=col2, label='weighted')
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
