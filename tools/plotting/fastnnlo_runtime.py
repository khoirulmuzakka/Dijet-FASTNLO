#!/usr/bin/env python3
#-*- coding:utf-8 -*-

import glob
import argparse
import glob
import os
import re
import sys
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.ticker import (FormatStrFormatter, LogFormatter, NullFormatter, ScalarFormatter, AutoMinorLocator, MultipleLocator)
from matplotlib import cm
# We do not want any interactive plotting! Figures are saved to files instead.
# This also avoids the ANNOYANCE of frequently missing Tkinter/tkinter (python2/3) GUI backends!
# To produce scalable graphics for publication use eps, pdf, or svg as file format.
# For this to work we try the Cairo backend, which can do all of these plus the raster format png.
# If this is not usable, we fall back to the Agg backend capable only of png for nice web plots.
#ngbackends = mpl.rcsetup.non_interactive_bk
#print('[fastnnlo_runtime]: Non GUI backends are: ', ngbackends)
# 1st try cairo
backend = 'cairo'
usecairo = True
try:
    import cairocffi as cairo
except ImportError:
    try:
        import cairo
    except ImportError:
        usecairo = False
#        print('[fastnnlo_runtime]: Can not use cairo backend :-(')
#        print('                   cairocffi or pycairo are required to be installed')
    else:
        if cairo.version_info < (1, 11, 0):
            # Introduced create_for_data for Py3.
            usecairo = False
#            print('[fastnnlo_runtime]: Can not use cairo backend :-(')
#            print('                   cairo {} is installed; cairo>=1.11.0 is required'.format(cairo.version))
if usecairo:
    mpl.use('cairo')
else:
    backend = 'agg'
    useagg = True
    try:
        mpl.use(backend, force=True)
    except:
        useagg = False
        raise ImportError('[PlotRuntime]: Neither cairo nor agg backend found :-( Cannot produce any plots. Good bye!')
    mpl.use('agg')
import matplotlib.pyplot as plt
# numpy
import numpy as np


# Redefine ScalarFormatter
class ScalarFormatterForceFormat(ScalarFormatter):
    # Override function that finds format to use.
    def _set_format(self, vmin, vmax):
        self.format = "%1.2f"  # Give format here


# Action class to allow comma-separated (or empty) list in options
class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values:
            setattr(namespace, self.dest, values[0].split(','))
        else:
            setattr(namespace, self.dest, [''])

# Some global definitions
_debug = False
_formats = {'eps': 0, 'pdf': 1, 'png': 2, 'svg': 3}
_channels = ['LO', 'R', 'V', 'RRa', 'RRb', 'RV', 'VV', 'ALL']
_channel_number = {'LO': 0, 'R': 1, 'V': 2, 'RRa': 3, 'RRb': 4, 'RV': 5, 'VV': 6}
_channel_colors = ['tab:green', 'tab:cyan', 'tab:blue', 'tab:red', 'tab:orange', 'tab:pink', 'tab:purple']

#####################################################################################

def main():

    # Print header
    print("\n###########################################################################################")
    print("# fastnnlo_runtime:")
    print("# Plot accumulated runtime of grid production")
    print("###########################################################################################\n")

    # Parse arguments
    args = arguments()

    # extract correct paths for input and outputfiles
    files = get_files(args['logfiles'])
    logfiles = []
    runfiles = []
    for file in files:
        runfile = re.sub('.log$', '.run', file)
        if not os.path.isfile(runfile):
            print('[fastnnlo_runtime]: WARNING! No matching runcard found for log file {}, skipped!')
        else:
            logfiles.append(file)
            runfiles.append(runfile)
    outputpath = args['output']
    print('[fastnnlo_runtime]: Output path argument is: {}'.format(outputpath))
    outputname = args['filename']
    print('[fastnnlo_runtime]: Filename argument is: {}'.format(outputname))

    # Plot formats to use
    formats = args['format']
    if formats is None:
        formats = ['png']
    for fmt in formats:
        if fmt not in _formats:
            print('[fastnnlo_runtime]: Illegal format specified, aborted!')
            print('[fastnnlo_runtime]: Format list:', args['format'])
            exit(1)
        elif fmt != 'png' and not usecairo:
            print('[fastnnlo_runtime]: Vector format plots not possible without cairo backend, aborted!')
            print('[fastnnlo_runtime]: Format list:', args['format'])
            exit(1)

    # get all the information from logfiles as dict
    # dict contains: runtime, runtime_unit, channel, events
    loginformation = get_loginformation(logfiles)
    runinformation = get_runinformation(runfiles)
    info = {**loginformation, **runinformation}

    # plot all the information
    if args['CPUtime']:
        plot_elapsed_time(info, outputpath, outputname, formats)
    if args['Events']:
        plot_events_per_hour(info, outputpath, outputname, formats)
    if not args['CPUtime'] and not args['Events']:
        plot_elapsed_time(info, outputpath, outputname, formats)
        plot_events_per_hour(info, outputpath, outputname, formats)

    exit(0)

def arguments():

    # Define arguments and options
    parser = argparse.ArgumentParser(epilog='Skript to plot elapsed time of fastNLO channels', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # positional argument
    parser.add_argument('logfiles', nargs='+', type=str,
                        help='Either input is a single string from LAW Task or it is a list')

    # optional arguments
    parser.add_argument('-f', '--filename', type=str,
                        help='Redefine part of output filename')
    parser.add_argument('-o', '--output', type=str, default='./',
                        help='Set here the outputpath')
    parser.add_argument('--CPUtime', dest='CPUtime', action='store_true',
                        help='Plot only the elapsed time')
    parser.add_argument('--Events', dest='Events', action='store_true',
                        help='Plot only the events per hour')
    parser.add_argument('--format', required=False, nargs=1, type=str, action=SplitArgs,
                        help='Comma-separated list of plot formats to use: eps, pdf, png, svg. If nothing is chosen, png is used.')
    return vars(parser.parse_args())


def get_files(files):
    # check if logfiles argument is one filename glob or a list of files
    if len(files)==1:
        files = glob.glob(files[0])
        if len(files)==1:
            print('fastnnlo_runtime: ERROR! Aborted, only one log file found: {}'.format(files[0]))
            exit(3)

    # sort unsorted glob list
    files.sort()
    return files

def get_loginformation(files):

    runtimes = []

    for file in files:
        runtimes_temp = []

        with open(file) as origin:
            for line in origin:
                # extract elapsed time with time unit
                if 'Time elapsed' in line:
                    line = line.split(':')
                    hours = float(line[1])
                    minutes = float(line[2])
                    seconds = float(line[3])

                    if hours != 0.:
                        runtimes_temp.append(hours + minutes/60 + seconds/360)
                        unit = 'hours'
                    else:
                        runtimes_temp.append(minutes + seconds/60)
                        unit = 'minutes'

        runtimes.append(runtimes_temp[-1])

    runtimes = np.array(runtimes)

    information = {
        'runtime': runtimes,
        'runtime_unit': unit
    }

    return information

def get_runinformation(files):

    nevents  = []
    channels = []

    for file in files:
        with open(file) as origin:
            for line in origin:
                # extract total events
                if 'Number of events' in line:
                    line = line.split()
                    nevents.append(line[0])
                if 'Job name id' in line:
                    line  = line.split()
                    parts = line[0].split('-')
                    channels.append(parts[0])

    nevents  = np.array(nevents)
    channels = np.array(channels)

    information = {
        'events': nevents,
        'channels': channels
    }

    return information

def plot_elapsed_time(infodict, out_path, out_name, formats):

    times = infodict['runtime']
    unit = infodict['runtime_unit']
    channels = infodict['channels']
    events   = infodict['events']

    # prepare histogram input
    unique_channels = set(channels)
    if len(unique_channels) == 0:
        print('fastnnlo_runtime: ERROR! Aborted, no channel info found.')
        exit(11)

    bins = np.linspace(min(times), max(times), 100)
    # get relevant values
    mean = np.mean(times)
    std = np.std(times)
    median = np.median(times)
    iqd = np.subtract(*np.percentile(times, [75, 25], interpolation='linear'))/2.

    CPUtime = np.sum(times) / (1 if unit == 'hours' else 60)

    # set figure
    fig = plt.figure(figsize=(16, 12))
    ax = fig.gca()

    # plot histogram

    if len(unique_channels) > 1 or [unique_channels] == 'ALL':
        n, batches, _ = ax.hist(times, bins=20, color='deepskyblue', edgecolor='black', label='Total CPU time: {0:0.0f} hours'.format(CPUtime))
        ax.legend(loc='best', fontsize=20)
    else:
        plt.text
        # plot each unique number in different color
        for ev_num in set(events):
            n, batches, _ = ax.hist(times[events == ev_num], histtype='barstacked', log=True, stacked=True, bins=bins, edgecolor='black', label='# Events: {}'.format(ev_num))

        # plot mean and median
        ax.vlines(mean, 0, max(n), colors='red', linestyles='dashed', label=r'Mean: {0:0.1f}$\pm${1:0.1f}'.format(mean, std))
        ax.vlines(median, 0, max(n), colors='green', linestyles='dashdot', label=r'Median: {0:0.2f}$\pm${1:0.2f}'.format(median, iqd))
        ax.ticklabel_format(axis='x', style='plain', useOffset=False)
        ax.legend(title='Total CPU time: {0:0.0f}hours'.format(CPUtime), loc='best', fontsize=20, title_fontsize=20)

    # finish and save figure
    chnlabel = channels[0]
    if out_name:
        chnlabel = out_name
    ax.set_title('Elapsed time of ' + chnlabel + ' production', fontsize=20)
    ax.set_xlabel('CPU time [' + unit + ']', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, fontsize=20)
    ax.set_ylabel('# jobs', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, fontsize=20, labelpad=20)
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=20)

    ax.grid()
    ax.set_axisbelow(True)

    # set saving location
    basename = 'runtime'
    for fmt in formats:
        filename = out_path + ('' if out_path[-1] == '/' else '/')
        if out_name:
            filename += out_name + '.' + basename + '.' + fmt
        else:
            filename += channels[0]  + '.' + basename + '.' + fmt
        print('[fastnnlo_runtime]: Saving runtime plot {}'.format(filename))
        fig.savefig(filename)

def plot_events_per_hour(infodict, out_path, out_name, formats):

    # get input
    channels = infodict['channels']
    events   = infodict['events']
    times    = infodict['runtime']
    unit     = infodict['runtime_unit']

    # prepare input
    unique_channels = set(channels)
    if len(unique_channels) == 0:
        print('fastnnlo_runtime: ERROR! Aborted, no channel info found.')
        exit(11)
    eph = []
    for i, time in enumerate(times):
        if unit == 'hours':
            eph.append(float(events[i])/time)
        else:
            eph.append(float(events[i])/(time/60))
    ephchn = []
    for i, val in enumerate(eph):
        for j, chn in enumerate(_channels):
            if channels[i] == chn:
                ephchn.append([val, j])
    ephchn = np.array(ephchn)
    masks = []
    for i, chn in enumerate(_channels):
        masks.append(ephchn[:,1] == i)

    # get relevant values
    mean    = np.mean(eph)
    std     = np.std(eph)
    median  = np.median(eph)
    iqd     = np.subtract(*np.percentile(eph, [75, 25], interpolation='linear'))/2.
    ephmin  = np.min(eph)
    ephmax  = np.max(eph)
    logbins = np.geomspace(ephmin, ephmax, 100)
    CPUtime = np.sum(times) / (1 if unit == 'hours' else 60)

    # create figure
    fig = plt.figure(figsize=(16, 12))
    ax = fig.gca()

    # plot (multistack-)histogram
    evrs = []
    lastch = 'LO'
    for chn in _channels:
        if chn in unique_channels:
            evrs.append(ephchn[masks[_channel_number[chn]]][:,0])
            lastch = chn
    if len(unique_channels) == 1:

        # plot each unique number in different color
        for ev_num in set(events):
            n, batches, _ = ax.hist(evrs[0][events == ev_num], histtype='barstacked', log=True, stacked=True, bins=50, edgecolor='black', label='# Events: {}'.format(ev_num))

        # plot mean and median
        ax.vlines(mean, 0, max(n), colors='red', linestyles='dashed', label=r'Mean: {0:0.1e}$\pm${1:0.1e} events/hour'.format(mean, std))
        ax.vlines(median, 0, max(n), colors='green', linestyles='dashdot', label=r'Median: {0:0.2e}$\pm${1:0.2e} events/hour'.format(median, iqd))
        ax.ticklabel_format(axis='x', style='plain', useOffset=False)
    else:
        n, batches, _ = ax.hist(evrs, histtype='barstacked', log=True, stacked=True, bins=logbins, edgecolor='black', color=_channel_colors, label=_channels)
        ax.set_xlim(0.9*ephmin, 1.1*ephmax)
        ax.set_xscale('log')

    # finish and save figure
    chnlabel = channels[0]
    if out_name:
        chnlabel = out_name
    ax.set_title('Event rate of ' + chnlabel + ' production', fontsize=20)
    ax.set_xlabel('event rate [1/hour]', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, fontsize=20)
    ax.set_ylabel('# jobs', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, fontsize=20, labelpad=20)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.legend(loc='best', fontsize=20)
    ax.grid()
    ax.set_axisbelow(True)

    # set saving location
    basename = 'evtrate'
    for fmt in formats:
        filename = out_path + ('' if out_path[-1] == '/' else '/')
        if out_name:
            filename += out_name + '.' + basename + '.' + fmt
        else:
            filename += channels[0]  + '.' + basename + '.' + fmt
        print('[fastnnlo_runtime]: Saving event rate plot {}'.format(filename))
        fig.savefig(filename)

if __name__ == "__main__":
    main()
