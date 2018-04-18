#!/usr/bin/env python2
#-*- coding:utf-8 -*-
import datetime
from datetime import timedelta
import glob, os, pylab, sys
import numpy as np
import re
from StringIO import StringIO

################################################################################
# Define method for histogram statistics
################################################################################
def hist_stats(x, weights, axis=None):

# Fill proper np.array in case only the number 1. is given as weight
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

# Minimum & maximum
    _min = np.min(x, axis=axis)
    _max = np.max(x, axis=axis)

    return dict(mean=_ave,
                stddev=_std,
                weighted_mean=_wave,
                weighted_stddev=_wstd,
                median=_med,
                iqd2=_med_err,
                mini=_min,
                maxi=_max)
################################################################################
# End of method for histogram statistics
################################################################################
# Define method to chop off microseconds from time differences
################################################################################
def chop_msec(delta):
    return delta - datetime.timedelta(microseconds=delta.microseconds)
################################################################################

# Default arguments
njobs = 100 # No. of jobs to evaluate

# Get from cmdline argument
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
if len(sys.argv) > 1:
    njobs = int(sys.argv[1])

# Prepare result arrays
mykeys = []
myvals = []

# Evaluate job info and vars files
ikeys = [ 'EVENTS', 'TIME', 'TIMESTAMP_EXECUTION_START', 'TIMESTAMP_EXECUTION_DONE', 'OUTPUT_FILE_0_SIZE' ]
okeys = [ 'Subprocess_', '# jobs_', '# events_[Mev.]', 'Tot. time_[h]', 'Tot. output_[GB]',
          'Ev./job_', 'Ev./hour_', '_min',
          '_-stddev', '<time/job>_hh:mm:ss', '_+stddev',
          '_-iqd/2', '{time/job}_hh:mm:ss', '_+iqd/2',
          '_max', 'Output/job_[MB]', 'Memory/job_[GB]' ]
jobdirs = glob.glob('output/job_*')
njob = 0
for jobdir in jobdirs:
    njob += 1
    jobnr = jobdir.split('_')[1]
    jobinfo = jobdir+'/'+'job.info'
    jobvars = jobdir+'/'+'job_'+jobnr+'.var'
    if njob<11:
        print 'Examining job info file ', jobinfo, ' and job var file ', jobvars
    else:
        print '.',
        sys.stdout.flush()

    myinfo = {}
    with open(jobinfo) as myfile:
        for line in myfile:
            name, var = line.partition("=")[::2]
            myinfo[name.strip()] = var.strip()

    myvars = {}
    with open(jobvars) as myfile:
        for line in myfile:
            myline = line.split(' ')[1]
            name, var = myline.partition("=")[::2]
            myvars[name.strip()] = var.strip()
        
    mykey = [ myinfo["JOBID"], myvars["CHN"].strip('"'), myvars["REG"].strip('"'), myvars["L"] ]
    myval = [ int(myvars[ikeys[0]].strip('"')), int(myinfo[ikeys[1]]), int(myinfo[ikeys[2]]), int(myinfo[ikeys[3]]), int(myinfo[ikeys[4]]) ]
    mykeys.append(mykey)
    myvals.append(myval)

    if njob > njobs:
        break

mykeys = np.array(mykeys)
myvals = np.array(myvals)

mysps = np.array([  "LO",   "R",   "V", "RR", "RR",  "RV", "VV" ])
myrgs = np.array([ "all", "all", "all",  "a",  "b", "all", "all" ])
myones = np.ones_like(myvals)
okeysa = [i.split('_')[0] for i in okeys]
okeysb = [i.split('_')[1] for i in okeys]
fmt = '{:>11} '*len(okeysa)

print
print '#==========================================================================================================================================================================================================='
head1 = fmt.format(*okeysa)
head2 = fmt.format(*okeysb)
print head1
print head2
print '#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
for i in range(mysps.size):
    spbl = (mykeys[:,1]==mysps[i]) * (mykeys[:,2]==myrgs[i])
    mysum = myvals[spbl].sum(axis=0)
    mycnt = myones[spbl].sum(axis=0)

    if mycnt[0] > 0:
        stats = {}
        if myrgs[i] == 'all':
            stats[okeys[0]] = mysps[i]
        else:
            stats[okeys[0]] = mysps[i]+myrgs[i]
        stats[okeys[1]] = mycnt[0]
        stats[okeys[2]] = mysum[0]/1000000
        stats[okeys[3]] = mysum[1]/3600
        stats[okeys[4]] = mysum[4]/1024**3
        stats[okeys[5]] = mysum[0]/mycnt[0]
        stats[okeys[6]] = mysum[0]/mysum[1]*3600
# Timing statistics
        col = myvals[:, [1]][spbl]
        _timing_stats = hist_stats(col, weights=1., axis=0)
        stats[okeys[7]]  = chop_msec(timedelta(seconds=_timing_stats['mini'][0]))
        stats[okeys[8]]  = chop_msec(timedelta(seconds=(_timing_stats['mean'][0] - _timing_stats['stddev'][0])))
        stats[okeys[9]]  = chop_msec(timedelta(seconds=_timing_stats['mean'][0]))
        stats[okeys[10]] = chop_msec(timedelta(seconds=(_timing_stats['mean'][0] + _timing_stats['stddev'][0])))
        stats[okeys[11]] = chop_msec(timedelta(seconds=(_timing_stats['median'][0] - _timing_stats['iqd2'][0])))
        stats[okeys[12]] = chop_msec(timedelta(seconds=_timing_stats['median'][0]))
        stats[okeys[13]] = chop_msec(timedelta(seconds=(_timing_stats['median'][0] + _timing_stats['iqd2'][0])))
        stats[okeys[14]] = chop_msec(timedelta(seconds=_timing_stats['maxi'][0]))
        stats[okeys[15]] = mysum[4]/mycnt[0]/1024**2
        stats[okeys[16]] = 'not av.:-('

# Prepare list output for printout
        out = []
        for key in okeys:
            stats[key] = str(stats[key])
            out.append(stats[key])
            
        line = fmt.format(*out)
        print line
print '#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'

exit(0)
