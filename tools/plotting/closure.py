#!/usr/bin/env python2
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os,sys
import fastnlo



def main():
    """
    PLotting script to create closure plots


    Enter nnlojet file and fastNLO table for comparison.
    Contributions can be deactivated in the script.
    
    PDF = 'NNPDF31_nnlo_as_0118'

    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=main.__doc__ )
    parser.add_argument('-f', '--nnlojetfile', help = 'nnlojet file (.dat) to be evaluated')
    parser.add_argument('-t', '--table', help = 'fastNLO table (.tab.gz) to be evaluated')
    args = parser.parse_args()

    fastnlo_table = args.table
    nnlojet_file = args.nnlojetfile

    fastnlo = get_fastnlo(str(fastnlo_table))
    nnlojet, unc = get_nnlojet(str(nnlojet_file))
    bins = []

    for i in range(len(fastnlo)):
        bins.append(i)

    #Title
    save_as = str(nnlojet_file)
    while(save_as.count('/')>0):
        index = save_as.find('/')
        save_as = save_as[index + 1:]
        save_as = save_as.replace('.dat','')

    color_opt = 'b'
    marker_opt = '*'
    fig, axs = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0) # Remove horizontal space between axes
    plt.subplot(211)
    plt.title(save_as)
    plt.yscale('log')
    plt.ylabel('Cross section [pb]')
    plt.plot(bins, fastnlo, marker=marker_opt, color='r', linestyle='dashed', linewidth = 0, markersize=5, label = '.tab.gz')
    bins_shift = []     #Data points get slightly shifted
    for i in range(len(bins)):
        bins_shift.append(bins[i]+0.2)
    plt.plot(bins_shift, nnlojet, marker=marker_opt, color=color_opt, linestyle='dashed', linewidth = 0, markersize=5, label = '.dat')
    plt.errorbar(bins_shift, nnlojet, yerr=unc, xerr=None, ecolor = color_opt, elinewidth = 1, capsize = 5, linestyle='')
    plt.legend()

    ratio(fastnlo, nnlojet, unc, bins, color_opt, marker_opt)

    save_as = save_as+'.eps'
    #save_as = 'closure/'+save_as
    plt.savefig(save_as)
    print('File saved as:')
    print(save_as)


def get_fastnlo(fnlofile):

    fnlo = fastnlo.fastNLOLHAPDF(fnlofile, 'NNPDF31_nnlo_as_0118', 0)
    fnlo.CalcCrossSection()
    fnlo.PrintCrossSections()

    #Set contribuions
    fnlo.SetContributionON(fastnlo.kFixedOrder,0,True)     #LO, True = On, False = Off
    fnlo.SetContributionON(fastnlo.kFixedOrder,1,True)     #NLO
    fnlo.SetContributionON(fastnlo.kFixedOrder,2,True)     #NNLO
    fnlo.CalcCrossSection()
    #fnlo.SetConstantAsScaleMuR/MuF

    xs = fnlo.GetCrossSection()
    return xs


def get_nnlojet(nnlofile):
    x,y = np.loadtxt(nnlofile, comments='#', usecols=(3,4), unpack= True)
    for i in range(len(x)):   #Cross section in pbarn
        x[i]=x[i]/1000
        y[i]=y[i]/1000
    
    return x,y


def ratio(fastnlo, nnlojet, uncertainty, bins, color_opt, marker_opt): 
    #Ratio plot including stat. uncertainties
    plt.subplot(212)
    plt.yscale('linear')
    plt.xlabel('Bin number')
    plt.ylabel('Ratio to nnlojet')
    value = []
    inverse = []
    bins_shift = []
    for i in range(len(nnlojet)):
        uncertainty[i] = uncertainty[i]/nnlojet[i]
        inverse.append(1/(nnlojet[i]/fastnlo[i]))    #Ratio tab.gz to dat file
        value.append(1)                              #Value set to 1, for comparison
    for i in range(len(bins)):
        bins_shift.append(bins[i]+0.2)
    plt.plot(bins, inverse, color =  'r', marker = marker_opt, linewidth = 0, markersize=5, label = '.tab.gz')
    plt.plot(bins_shift, value, color =  color_opt, marker = marker_opt, linewidth = 0, markersize=5, label = '.dat')
    plt.errorbar(bins_shift, value, yerr=uncertainty, xerr=None, ecolor = color_opt, elinewidth = 1, capsize = 5, linestyle='')

    







if __name__== "__main__":
  main()

