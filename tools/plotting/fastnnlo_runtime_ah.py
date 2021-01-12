#!/usr/bin/env python3
#-*- coding:utf-8 -*-

import glob
import argparse
import glob
import sys
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.ticker import (FormatStrFormatter, LogFormatter, NullFormatter, ScalarFormatter, AutoMinorLocator, MultipleLocator)
from matplotlib import cm
# We do not want any interactive plotting! Figures are saved to files instead.
# This also avoids the ANNOYANCE of frequently missing Tkinter/tkinter (python2/3) GUI backends!
# To produce scalable graphics for publication use eps, pdf, or svg as file format.
# For this to work we try the Cairo backend, which can do all of these plus the raster format png.
# If this is not usable, we fall back to the Agg backend capable only of png for nice web plots.
#ngbackends = mpl.rcsetup.non_interactive_bk
#print('[fastnnlo_pdfunc]: Non GUI backends are: ', ngbackends)
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
#        print('[fastnnlo_pdfunc]: Can not use cairo backend :-(')
#        print('                   cairocffi or pycairo are required to be installed')
    else:
        if cairo.version_info < (1, 11, 0):
            # Introduced create_for_data for Py3.
            usecairo = False
#            print('[fastnnlo_pdfunc]: Can not use cairo backend :-(')
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



def main():

    # Parse arguments
    args = arguments()

    # extract correct paths for input and outputfiles
    logfiles = get_files(args['logfiles'])
    outputpath = args['output']
    all_plot = args['All']
    fmt = args['format']

    # get all the information from logfiles as dict
    # dict contains: runtime, runtime_unit, channel, events
    loginformation = get_loginformation(logfiles)

    # plot all the information
    if args['CPUtime']:
        plot_elapsed_time(loginformation, outputpath, all_plot, fmt)
    if args['Events']:
        plot_events_per_hour(loginformation, outputpath, all_plot, fmt)
    if not args['CPUtime'] and not args['Events']:
        plot_elapsed_time(loginformation, outputpath, all_plot, fmt)
        plot_events_per_hour(loginformation, outputpath, all_plot, fmt)

    exit(0)

def arguments():

    # Define arguments and options
    parser = argparse.ArgumentParser(epilog='Skript to plot elapsed time of fastNLO channels', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # positional argument
    parser.add_argument('logfiles', nargs='+', type=str, help='Either input is a single string from LAW Task or it is a list')

    # optional arguments
    parser.add_argument('-o', '--output', nargs=1, type=str, default='./',
                        help='Set here the outputpath')
    parser.add_argument('--CPUtime', dest='CPUtime', action='store_true',
                        help='Plot only the elapsed time')
    parser.add_argument('--Events', dest='Events', action='store_true',
                        help='Plot only the events per hour')
    parser.add_argument('--All', dest='All', action='store_true',
                        help='Set plotoptions according to plot for all channels')
    parser.add_argument('--format', required=False, nargs=1, type=str, default='png',
                        help='Comma-separated list of plot formats to use: eps, pdf, png, svg. If nothing is chosen, png is used.')

    return vars(parser.parse_args())


def get_files(files):
    # check if logfiles argument is one filename glob or a list of files
    if len(files)==1:
        files = glob.glob(files[0])
        if len(files)==1:
            print('fastnnlo_runtime: ERROR! Aborted, only one log file found: {}'.format(files[0]))
            exit(3)

    return files

def get_loginformation(files):

    run_time = []
    number_events = []
    channel = None

    for file in files:
        event = False
        run_time_temp = []

        with open(file) as origin:
            for line in origin:
                # extract elapsed time with time unit
                if 'Time elapsed' in line:
                    line = line.split(':')
                    hours = float(line[1])
                    minutes = float(line[2])
                    seconds = float(line[3])

                    if hours == True:
                        run_time_temp.append(hours + minutes/60 + seconds/360)
                        unit = 'hours'
                    else:
                        run_time_temp.append(minutes + seconds/60)
                        unit = 'minutes'

                # extract channel name
                if 'Tablename' in line and not channel:
                    line = line.split()
                    tablename = line[2].split('.')
                    channel = tablename[0] + '.' + tablename[1]
                # extract total events
                if 'ncalltot=' in line and not event:
                    line = line.split(',')
                    number_events.append(float(line[4][10:]))
                    event = True

        run_time.append(run_time_temp[-1])


    run_time = np.array(run_time)
    number_events = np.array(number_events)

    information = {
        'runtime': run_time,
        'runtime_unit': unit,
        'channel': channel,
        'events': number_events
    }

    return information


def plot_elapsed_time(informationdict, out_path, plot_all=False, format='png'):

    time = informationdict['runtime']
    unit = informationdict['runtime_unit']
    channel = informationdict['channel']

    if plot_all:
        channel = channel.split('.')
        channel = channel[0]

    # get relevant values
    mean = np.mean(time)
    std = np.std(time)
    median = np.median(time)
    iqd = np.subtract(*np.percentile(time, [75, 25], interpolation='linear'))/2.

    CPUtime = np.sum(time) / (1 if unit == 'hours' else 60)

    # set saving location
    filename = out_path[0] + ('' if out_path[0][-1] == '/' else '/')
    filename += channel + '.Hist_Elapsed_time.' + format

    # set figure
    fig = plt.figure(figsize=(16, 12))
    ax = fig.gca()

    # plot histogram
    n, batches, _ = ax.hist(time, bins=20, color='deepskyblue', edgecolor='black', label='Total CPU time: {0:0.0f} hours'.format(CPUtime))

    if not plot_all:
        # plot mean and median
        ax.vlines(mean, 0, max(n), colors='red', linestyles='dashed', label=r'Mean: {0:0.1f}$\pm${2:0.1f} {1}'.format(mean, unit, std))
        ax.vlines(median, 0, max(n), colors='green', linestyles='dashed', label=r'Median: {0:0.1f}$\pm${2:0.1f} {1}'.format(median, unit, iqd))

    # finish and save figure
    ax.set_title('Elapsed time of ' + channel + ' production', fontsize=20)
    ax.set_xlabel('CPU time [' + unit + ']', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, fontsize=20, labelpad=15)
    ax.set_ylabel('frequency', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, fontsize=20, labelpad=15)
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=20)

    ax.legend(loc='best', fontsize=20)
    ax.grid()
    ax.set_axisbelow(True)

    fig.savefig(filename)

def plot_events_per_hour(informationdict, out_path, plot_all=False, format='png'):

    time = informationdict['runtime']
    unit = informationdict['runtime_unit']
    channel = informationdict['channel']
    events = informationdict['events']

    if unit == 'hours':
        eph = events/time
    else:
        eph = events/(time/60)

    if plot_all:
        channel = channel.split('.')
        channel = channel[0]

    # get relevant values
    mean = np.mean(eph)
    std = np.std(eph)
    median = np.median(eph)
    iqd = np.subtract(*np.percentile(eph, [75, 25], interpolation='linear'))/2.

    CPUtime = np.sum(time) / (1 if unit == 'hours' else 60)


    # set saving location
    filename = out_path[0] + ('' if out_path[0][-1] == '/' else '/')
    filename += channel + '.Hist_Events_per_hour.' + format

    # set figure
    fig = plt.figure(figsize=(16, 12))
    ax = fig.gca()

    # plot histogram
    n, batches, _ = ax.hist(eph, bins=20, color='deepskyblue', edgecolor='black', label='Total CPU time: {0:0.0f} hours'.format(CPUtime))

    if not plot_all:
        # scientific format
        f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
        g = lambda x,pos : "${}$".format(f._formatSciNotation('%0.2e' % x))
        sci = mticker.FuncFormatter(g)

        # plot mean and median
        ax.vlines(mean, 0, max(n), colors='red', linestyles='dashed', label=r'Mean: {0}$\pm${1} events/hour'.format(sci(mean), sci(std)))
        ax.vlines(median, 0, max(n), colors='green', linestyles='dashed', label=r'Median: {0}$\pm${1} events/hour'.format(sci(median), sci(iqd)))

    # finish and save figure
    ax.set_title('Events per hour of ' + channel + ' production', fontsize=20)
    ax.set_xlabel('events/hour', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, fontsize=20, labelpad=15)
    ax.set_ylabel('frequency', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, fontsize=20, labelpad=15)
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=20)

    ax.legend(loc='best', fontsize=20)
    ax.grid()
    ax.set_axisbelow(True)

    fig.savefig(filename)


if __name__ == "__main__":
    main()
