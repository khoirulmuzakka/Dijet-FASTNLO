#!/usr/bin/env python2
#-*- coding:utf-8 -*-
#
# Make statistical evaluation plots of ensemble of one-to-one comparisons
# between fastNLO interpolation tables and NNLOJET original results
#
# Version:
#
# created by K. Rabbertz: 13.07.2017
#
#-----------------------------------------------------------------------
#
import argparse
import glob, os, pylab, sys
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
# from copy import deepcopy
from matplotlib import cm
from fastnlo import fastNLOLHAPDF
from fastnlo import SetGlobalVerbosity
import re
from StringIO import StringIO
#import warnings
#warnings.filterwarnings("error")

parser = argparse.ArgumentParser()

#add arguments (three .dat and one .log files)
parser.add_argument('-d','--datfile', default='file.dat', required=True, nargs='?',
                    help='.dat file for evaluation. All other matching .dat files that '
                            'only differ in seed no. will be evaluated as well.')
parser.add_argument('-w','--weightfile', default='weight.txt', required=True, nargs='?',
                    help='.txt file containing weights.')
#parser.add_argument('-l', '--logfile', default='file.log', required=True, nargs='?',
#                    help='.log files, need one *0.log and one *6.log')
parser.add_argument('-f','--fscl', default=1, required=True, nargs='?', type=int,
                    help='1: fixed scale.  \n'
                            '0: flexible scale.') ### Has to be changed!!! Only allow: 1 to 7

parser.add_argument('-o', '--outputfilename', required=False, nargs='?', type=str,
                    help='Customise the first part of the output filename.'
                            'Default: Same structure as datfile name.')


#parse arguments
args = vars(parser.parse_args())
namesp = parser.parse_args()

''' ###NOT needed in this script, can be removed in final version
#take care of the input files
#log0file = args['logfile']
#create dictionary for order and correspoinding datfiles
datfiles = {}
for arg in args['datfiles']:
    pre = arg.split('.')[0]
    if pre in ['LO', 'NLO', 'NNLO']:
        datfiles[pre]=arg
    else:
        print "Problem with datfiles."
        sys.exit("Exit: Input ERROR.")

print '\n', "datfiles: ", datfiles, '\n'
if len(datfiles)!=3:
    print "Check whether all required datfiles are given."
    sys.exit("ERROR: Wrong amount of .dat files per order.")
    #for instance if user gives two times the NLO .dat file and no NNLO

'''



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
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large'}
pylab.rcParams.update(params)

#
# Start
#
print "\n######################################################################################"
print "# fastnnlo_approxtest.py: Plot statistical evaluation of fastNLO interpolation quality"
print "#######################################################################################\n"

#
# Default arguments
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
#arguments
print "Given datfile: ", datfile #, '\n'
datbase = os.path.basename(datfile)
datargs = datbase.split(".") #array containing filename-parts

print "Given weightfile: ", wgtfile #, '\n'
if fscl==1: print "Fixed scale. \n"
elif fscl==0: print "Flexible scale. \n"


#arguments from datfile
print 'arguments in datfile name: ', len(datargs)
print "datargs: ", datargs, '\n'
if len(datargs)==4:
    proc, jobn, obsv, ext = datargs
    #obsv = obsv[:-2]
elif len(datargs)==5:
    proc, jobn, kinn, obsv, ext = datargs
    #obsv = obsv[:-2] #cut the _0 part
elif len(datargs)==6:
    proc, jobn, kinn, obsv, seed, ext = datargs
    #seed = seed[:-2] #cut the _0 part --> is not there in datfiles



print 'proc: ', proc
print 'jobn: ', jobn
print 'kinn: ', kinn
print 'obsv: ', obsv
print 'seed: ', seed
print '\n'
#

'''
# TODO: Do some more elaborate argument treatment and help texts
print 'Usage:    fastnnlo_approxtest process jobname kinvegas obsv nmax fscl'
print 'Defaults:',proc, jobn, kinn, obsv, nmax, fscl
print 'Number of arguments given:', len(sys.argv), 'arguments.'
print 'Argument list:', str(sys.argv)
if len(sys.argv) > 1:
    proc = sys.argv[1]
if len(sys.argv) > 2:
    jobn = sys.argv[2]
if len(sys.argv) > 3:
    kinn = sys.argv[3]
if len(sys.argv) > 4:
    obsv = sys.argv[4]
if len(sys.argv) > 5:
    if sys.argv[5] != '_':
        nmax = int(sys.argv[5])
if len(sys.argv) > 6:
    fscl = int(sys.argv[6])
print 'Actual:  ',proc, jobn, kinn, obsv, nmax, fscl
'''


# Extract order/contribution from job type (substring before first '-')
order  = jobn.split('-')[0]

# Prepare result arrays
xl = []      # left bin border
xm = []      # bin "center"
xu = []      # right bin border
ndat = 0
xs_nnlo = [] # NNLOJET results
ntab = 0
xs_fnlt = [] # fastNLO results
nlog = 0
xs_fnll = [] # Pre-evaluated fastNLO results
nlin = 0
weights = [] # Weight factors
seeds = []   # Seed numbers for matching

# Read binning and cross sections from NNLOJET dat file
datglob  = order+'/'+proc+'.'+jobn+'.'+kinn+'.'+obsv+'.s*.dat' #PROBLEM: this is tied to a certain structure (quite inconvenient), should be more flexible according to given data file
print 'datglob old', datglob

#datglob = datfile[:-9]+'*.dat' #or better: cut extension and seed and add '.s*.dat'
datglob = datfile[:-(len(seed)+len(ext))]+'*.dat'
#datglob = order+'/'+datbase[:-9]+'*.dat'
print 'datglob new', datglob, '\n'


datfiles = glob.glob(datglob)
if not datfiles:
    print >> sys.stderr, 'No NNLOJET dat files matching', datglob ,'found, aborted!'
    sys.exit(1) 
datfiles.sort()
ixscol = 3 + 2 * (fscl-1)
for datfile in datfiles:
    print 'Datfile no. ', ndat, ' is ', datfile
    if ndat == 0:
        xl.append(np.loadtxt(datfile,usecols=(0,)))
        xu.append(np.loadtxt(datfile,usecols=(2,)))
    xs_nnlo.append(np.loadtxt(datfile,usecols=(ixscol,)))
    print "xs_nnlo \n", xs_nnlo ## just to see what's happening
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

# Read weights per file per bin from Alex
#wgtfile = 'Combined/Final/'+order+'.'+obsv+'.'+'APPLfast.txt'  #should be accessible more flexibly
#wgtfile = args['weightfile'] #handle this already in the beginning with the other args
wgtnams = np.genfromtxt(wgtfile, dtype=None, usecols=0)
wgttmps = np.loadtxt(wgtfile, usecols=(list(range(1,nobs+1))))
ntmp = len(wgtnams)
# Combine to key-value tuple ( name, weights[] ) and sort according to key=name
wgttup = zip(wgtnams,wgttmps)
wgttup.sort(key = lambda row: (row[0]))
# Unzip again
allnames,allweights = zip(*wgttup)
print 'datfile names in weightfile: \n', np.array(allnames), '\n'

#################################################
for dfile in datfiles:
    dirna = os.path.dirname(dfile)
    basena = os.path.basename(dfile)
    newna = './'+order+'/'+basena
    print 'dirname', dirna, ' --- basename ', basena
    print './+order+basename', newna, '\n' #that's the way they're listed in wgtfile

print 'datfiles array: \n', np.array(datfiles), '\n'

####################HERE IS THE CURRENT DIFFICULTY --> IMPLEMENT NEW HANDLING OF ORDER STRUCTURE
#make it possible to use less .dat files than entries in the weightfile
#select the corresponding weightfile-entries

###new attempt for different handling:
#ncount = 0 #used to be handled via nlin, but now nlin has new meaning
            #nlin refers to line number only temporarily
            #whereas ncount counts how many lines of the weightfiles are used
for dfile in datfiles:
    newna = './'+order+'/'+os.path.basename(dfile)
    indexlin = allnames.index(newna)
    print 'Weight file line no. ', indexlin, ' is for ' , allnames[indexlin]
    weights.append(allweights[indexlin])
    #ncount+=1
    nlin+=1
weights = np.array(weights)

print '\n', 'Using %s weight lines.' %nlin  #%ncount
print '\n','weights-array: \n', weights,'\n'



''' ##Do we always need all the seeds?
for name in allnames:
    print 'Weight file line no. ', nlin, ' is for ', name 
    if seeds[nlin] not in name:
        print 'seeds[',nlin,'] = ', seeds[nlin],', weight file name = ', name
        sys.exit('ERROR: Mismatch in result sort order between NNLOJET and NNLOJET weights. Aborted!')
    else:
        weights.append(allweights[nlin])
    nlin += 1
    if nlin == nmax:
        break
weights = np.array(weights)
print 'Using ', nlin, 'weight lines.'
'''

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
if fscl == 1:
    #fnlologs = glob.glob(order+'/'+proc+'.'+jobn+'.'+kinn+'.'+obsv+'.s*_0.log')
    
    #fnlologsname = order+'/'+proc+'.'+jobn+'.'+kinn+'.'+obsv+'.s*_0.log'
    #print fnlologsname
    #fnlologs_new = glob.glob(datfiles[0][:-3]+'_0.log')
    #print datfiles[0][:-(len(ext)+1)]
    #print fnlologs_new
    
    for dfile in datfiles:
        log_list.append(dfile[:-(len(ext)+1)]+'_0.log')
    fnlologs = np.array(log_list)
    print 'fnlologs: \n', fnlologs, '\n'
    
else:
    #fnlologs = glob.glob(order+'/'+proc+'.'+jobn+'.'+kinn+'.'+obsv+'.s*_6.log')
    for dfile in datfiles:
        log_list.append(dfile[:-(len(ext)+1)]+'_6.log')
    fnlologs = np.array(log_list)
    print 'fnlologs: \n', fnlologs, '\n'

fnlologs.sort()
for fnlolog in fnlologs:
    print 'fastNLO file no. ', nlog, ' is ', fnlolog
    if not seeds[nlog] in fnlolog:
        print 'Mismatch in result sort order between NNLOJET and fastNLO. Aborted!'
        print 'seeds[',nlog,'] = ', seeds[nlog],', fastNLO file = ', fnlolog
    else:
        print 'NNLOJET and fastNLO result correctly matched. seed is ', seeds[nlog]
    xs_tmp = np.loadtxt(fnlolog,usecols=(6,),comments=['#',' #','C','L'])
    if fscl == 1:
        indi = 0
    else:
        indi = (fscl-2)*nobs
    indf = indi + nobs
    xs_sub = xs_tmp[indi:indf]
    xs_fnll.append(xs_sub)
    nlog += 1
    if nlog==nmax: break
print 'Using ', nlog, 'pre-evaluated table files.'
xs_fnll = np.array(xs_fnll)

print 'xs_fnll', xs_fnll ####just to have a look at it

# Check on identical file numbers, either for table or log files and weight lines
if ndat == nlog == nlin and ndat*nlog*nlin != 0:
    print 'OK: Have equal number of dat files, fastNLO results, and weight lines: ndat = ', ndat, ', nlog = ', nlog, ', nlin = ', nlin
#if ndat == nlog == ncount and ndat*nlog*ncount != 0:
    #print 'OK: Have equal number of dat files, fastNLO results, and weight lines: ndat = ', ndat, ', nlog = ', nlog, ', ncount = ', ncount
else:
#    sys.exit('ERROR: No matching file found or file number mismatch! ndat = '+str(ndat)+', ntab = '+str(ntab))
    sys.exit('ERROR: No matching file found or file number mismatch! ndat = '+str(ndat)+', nlog = '+str(nlog)+', nlin = '+str(nlin))
    #sys.exit('ERROR: No matching file found or file number mismatch! ndat = '+str(ndat)+', nlog = '+str(nlog)+', ncount = '+str(ncount))
#r_f2nt = xs_fnlt/xs_nnlo
r_f2nl = xs_fnll/xs_nnlo
a_f2nl = (xs_fnll - xs_nnlo)/(xs_fnll + xs_nnlo)

print 'r_f2nl', r_f2nl ####just to have a look at it
print 'a_f2nl', a_f2nl ####just to have a look at it

# Ratio statistics
_ratio_stats  = hist_stats(r_f2nl, weights=weights, axis=0)
rave_f2nl    = _ratio_stats['mean']
rstd_f2nl    = _ratio_stats['stddev']
rwav_f2nl    = _ratio_stats['weighted_mean']
rwst_f2nl    = _ratio_stats['weighted_stddev']
rmed_f2nl    = _ratio_stats['median']
riqd_f2nl    = _ratio_stats['iqd2']

print 'rave_f2nl', rave_f2nl ####just to have a look at it

# Asymmetry statistics
_asymm_stats = hist_stats(a_f2nl, weights=weights, axis=0)
aave_f2nl = _asymm_stats['mean']
astd_f2nl = _asymm_stats['stddev']
awav_f2nl = _asymm_stats['weighted_mean']
awst_f2nl = _asymm_stats['weighted_stddev']
amed_f2nl = _asymm_stats['median']
aiqd_f2nl = _asymm_stats['iqd2']

print 'awav_f2nl', awav_f2nl ####just to have a look at it
print 'aave_f2nl', aave_f2nl ####just to have a look at it

# Ratio plot
limfs = 'x-large'
fig = plt.figure(figsize=(16,12))
ax  = fig.gca()
plt.title(r'Ratio: {} {} for {} using scale {}'.format(proc, order, obsv, fscl), fontsize='x-large', fontweight='bold', loc='left')
plt.xlabel('Observable bin index', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
plt.ylabel('fastNLO/NNLOJET', horizontalalignment='right', x=1.0, verticalalignment='bottom', y=1.0)
plt.axhline(y=1.001, linestyle='--', linewidth=1.0, color='black')
plt.axhline(y=0.999, linestyle='--', linewidth=1.0, color='black')
plt.fill_between([0.0,34.0],0.999,1.001, color='black', alpha=0.1)
plt.text(33.6,1.00085,u'+1‰')
plt.text(33.7,0.99885,u'–1‰')

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

plt.xlim(0.0,34.0)
#plt.ylim(0.99,1.01)
plt.ylim(-50, 12000) #to display points for fscl=0 (log6file)
if nlin>0:
    handles = [rmean,rmwgt,rmedian]
else:
    handles = [rmean,rmedian]
labels  = [h.get_label() for h in handles]
legend = ax.legend(handles, labels, title=r'No. of entries/bin = {}'.format(ndat), loc='upper left', numpoints=1, handlelength=0)
legend.get_title().set_fontsize(limfs)

if args['outputfilename'] is None:
    fignam = proc+'.'+jobn+'.'+obsv+'.'+'fscl'+str(fscl)+'.'+'approx_ratio_'+str(ndat)+'.png'
else:
    fignam = args['outputfilename']+'.approx_ratio_'+str(ndat)+'.png'
plt.savefig(fignam)
#plt.show()

# Asymmetry plot
limfs = 'x-large'
fig = plt.figure(figsize=(16,12))
ax  = fig.gca()
plt.title(r'Asymmetry: {} {} for {} using scale {}'.format(proc, order, obsv, fscl), fontsize='x-large', fontweight='bold', loc='left')
plt.xlabel('Observable bin index', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
plt.ylabel('(fastNLO-NNLOJET)/(fastNLO+NNLOJET)', horizontalalignment='right', x=1.0, verticalalignment='bottom', y=1.0)
plt.axhline(y=-0.001, linestyle='--', linewidth=1.0, color='black')
plt.axhline(y=+0.001, linestyle='--', linewidth=1.0, color='black')
plt.fill_between([0.0,34.0],-0.001,0.001, color='black', alpha=0.1)
plt.text(33.6,+0.00085,u'+1‰')
plt.text(33.7,-0.00115,u'–1‰')

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

plt.xlim(0.0,34.0)
#plt.ylim(-0.01,0.01)
plt.ylim(-1.02, 1.02) ## to display points for fscl=0 (log6file)
if nlin>0:
    handles = [amean,amwgt,amedian]
else:
    handles = [amean,amedian]
labels  = [h.get_label() for h in handles]
legend = ax.legend(handles, labels, title=r'No. of entries/bin = {}'.format(ndat), loc='upper left', numpoints=1, handlelength=0)
legend.get_title().set_fontsize(limfs)

if args['outputfilename'] is None:
    fignam = proc+'.'+jobn+'.'+obsv+'.'+'fscl'+str(fscl)+'.'+'approx_asymm_'+str(ndat)+'.png'
else:
    fignam = args['outputfilename']+'.approx_asymm_'+str(ndat)+'.png'
plt.savefig(fignam)
#plt.show()



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

################################################################################

exit(0)
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
