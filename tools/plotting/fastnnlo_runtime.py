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

# Action class to allow comma-separated (or empty) list in options
class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values:
            setattr(namespace, self.dest, values[0].split(','))
        else:
            setattr(namespace, self.dest, [''])

# Some global definitions
_debug = True
_formats = {'eps': 0, 'pdf': 1, 'png': 2, 'svg': 3}

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
    logfiles = get_files(args['logfiles'])
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

    # plot all the information
    if args['CPUtime']:
        plot_elapsed_time(loginformation, outputpath, outputname, formats)
    if args['Events']:
        plot_events_per_hour(loginformation, outputpath, outputname, formats)
    if not args['CPUtime'] and not args['Events']:
        plot_elapsed_time(loginformation, outputpath, outputname, formats)
        plot_events_per_hour(loginformation, outputpath, outputname, formats)

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

                    if hours != 0.:
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

def plot_elapsed_time(informationdict, out_path, out_name, formats):

    time = informationdict['runtime']
    unit = informationdict['runtime_unit']
    channel = informationdict['channel']
    basename = 'runtime'

    # get relevant values
    mean = np.mean(time)
    std = np.std(time)
    median = np.median(time)
    iqd = np.subtract(*np.percentile(time, [75, 25], interpolation='linear'))/2.

    CPUtime = np.sum(time) / (1 if unit == 'hours' else 60)

    # set figure
    fig = plt.figure(figsize=(16, 12))
    ax = fig.gca()

    # plot histogram
    n, batches, _ = ax.hist(time, bins=20, color='deepskyblue', edgecolor='black', label='Total CPU time: {0:0.0f} hours'.format(CPUtime))

    # plot mean and median
    ax.vlines(mean, 0, max(n), colors='red', linestyles='dashed', label=r'Mean: {0:0.1f}$\pm${2:0.1f} {1}'.format(mean, unit, std))
    ax.vlines(median, 0, max(n), colors='green', linestyles='dashed', label=r'Median: {0:0.1f}$\pm${2:0.1f} {1}'.format(median, unit, iqd))

    # finish and save figure
    chnlabel = channel
    if out_name:
        chnlabel = out_name
    ax.set_title('Elapsed time of ' + chnlabel + ' production', fontsize=20)
    ax.set_xlabel('CPU time [' + unit + ']', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, fontsize=20)
    ax.set_ylabel('frequency', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, fontsize=20)
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=20)

    ax.legend(loc='best', fontsize=20)
    ax.grid()
    ax.set_axisbelow(True)

    # set saving location
    for fmt in formats:
        filename = out_path + ('' if out_path[-1] == '/' else '/')
        if out_name:
            filename += out_name + '.' + basename + '.' + fmt
        else:
            filename += channel  + '.' + basename + '.' + fmt
        print('[fastnnlo_runtime]: Saving runtime plot {}'.format(filename))
        fig.savefig(filename)

def plot_events_per_hour(informationdict, out_path, out_name, formats):

    time = informationdict['runtime']
    unit = informationdict['runtime_unit']
    channel = informationdict['channel']
    events = informationdict['events']
    basename = 'evtrate'

    if unit == 'hours':
        eph = events/time
    else:
        eph = events/(time/60)

    # get relevant values
    mean = np.mean(eph)
    std = np.std(eph)
    median = np.median(eph)
    iqd = np.subtract(*np.percentile(eph, [75, 25], interpolation='linear'))/2.

    CPUtime = np.sum(time) / (1 if unit == 'hours' else 60)

    # set figure
    fig = plt.figure(figsize=(16, 12))
    ax = fig.gca()

    # plot histogram
    n, batches, _ = ax.hist(eph, bins=20, color='deepskyblue', edgecolor='black', label='Total CPU time: {0:0.0f} hours'.format(CPUtime))

    # plot mean and median
    ax.vlines(mean, 0, max(n), colors='red', linestyles='dashed', label=r'Mean: {0:0.1e}$\pm${1:0.1e} events/hour'.format(mean, std))
    ax.vlines(median, 0, max(n), colors='green', linestyles='dashed', label=r'Median: {0:0.2e}$\pm${1:0.2e} events/hour'.format(median, iqd))

    # finish and save figure
    ax.set_title('Event rate of ' + channel + ' production', fontsize=20)
    ax.set_xlabel('event rate [1/hour]', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, fontsize=20)
    ax.set_ylabel('frequency', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, fontsize=20)
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=20)

    ax.legend(loc='best', fontsize=20)
    ax.grid()
    ax.set_axisbelow(True)

    # set saving location
    for fmt in formats:
        filename = out_path + ('' if out_path[-1] == '/' else '/')
        if out_name:
            filename += out_name + '.' + basename + '.' + fmt
        else:
            filename += channel  + '.' + basename + '.' + fmt
        print('[fastnnlo_runtime]: Saving event rate plot {}'.format(filename))
        fig.savefig(filename)

if __name__ == "__main__":
    main()
