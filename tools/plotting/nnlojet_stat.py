#!/usr/bin/env python3

#----------------------------------------------

#Python script to create plots

#Plot includes 4 panels:
# 1) Cross section vs. Bin
# 2) Ratio plot including a fit (polynomial, third degree)
# 3) Ratio compared to the fit
# 4) Statistical uncertainites

#Figure is saved as eps under the name of the NNLO file

#----------------------------------------------


import argparse
import numpy as np
from io import StringIO
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec
import textwrap


def crosssection(number, name, color_opt, marker_opt):

    x,y = np.loadtxt(name, comments='#', usecols=(3,4), unpack= True)
    bins = []
    for i in range(len(x)): 
        bins.append(i+int(number)*0.15)

    #Takes care of correct labeling
    if(name.count('/NLO.')==1):
        label = 'NLO'
    elif(name.count('NNLO')==1):
        label = 'NNLO'
        for i in range(7):
            string = 'y'+str(i)
            if(name.count(string)!=0):
                title = str(name)
                while(title.count('/')>0):
                    index = title.find('/')
                    title = title[index + 1:]
                title = title.replace('.dat','')
        plt.title(title)
        plt.ylabel('Cross section')

    else:
        label = 'LO'

    #Plot the data with stat uncertainies
    plt.plot(bins, x, marker=marker_opt, color=color_opt, linestyle='dashed', linewidth = 0, markersize=5, label = label)
    plt.errorbar(bins, x, yerr=y, xerr=None, ecolor = color_opt, elinewidth = 1, capsize = 5, linestyle='')
    

def ratio(number, default, name, color_opt, marker_opt, color_default, marker_default, order):
    x_d,y_d = np.loadtxt(default, comments='#', usecols=(3,4), unpack= True)
    x,y = np.loadtxt(name, comments='#', usecols=(3,4), unpack= True)
    bins=[]

    if(name.count('NLO.')==1):
        label = 'NLO'
    elif(name.count('NNLO')==1):
        label = 'NNLO'
    else:
        label = 'LO'
    
    for i in range(len(x)):
        bins.append(i)

    if(int(number)==int(order)):          #Ratio to the specified order
        plt.ylabel('Ratio to '+ label)
        for i in range(len(x)):
            y[i]=y[i]/x[i]
            x[i]=1
            bins[i] += number*0.15
    else:
        for i in range(len(x)):
            y[i]=y[i]/x[i]
            x[i]=x[i]/x_d[i]
            bins[i] += number*0.15

    plt.plot(bins, x, color =  color_opt, marker = marker_opt, linewidth = 0, markersize=5, label = label)
    plt.errorbar(bins, x, yerr=y, xerr=None, ecolor = color_opt, elinewidth = 1, capsize = 5, linestyle='')

    plt.subplot(414)
    plt.bar(bins, y, color=color_opt, width=0.8-number*0.2, alpha = 0.2*(number+1), edgecolor = color_opt, align='center', label = None)
    plt.ylabel('Rel. stat. unc.')
    plt.subplot(412)
    
    if(number!=order):
        function = fit_polynomial(bins, x, y, color_opt, marker_opt, label)
        plt.plot(bins, function, color =  color_opt, marker = None, linewidth = 1)


def fit_polynomial(x, y, y_error, color_opt, marker_opt, label):            #Fits a polynomial of third degree to the ratio plot
    plt.subplot(413)
    series = np.polynomial.polynomial.Polynomial.fit(x,y,3)
    #print(series.convert().coef)    #print fit coefficients
    coeffs = series.convert().coef
    function = 0
    bins = []

    for i in range(len(x)):
        bins.append(i)
    z = np.linspace(int(np.amin(bins)), int(np.amax(bins)), int(np.amax(bins))-int(np.amin(bins))+1)
    for i in range(len(coeffs)):
        function += coeffs[i]*z**i
    function_const=[]
    for i in  range(len(function)):
        y[i] = y[i]/function[i]
        function_const.append(1.) 
            
    plt.plot(x, y, color =  color_opt, marker = marker_opt, linewidth = 0, markersize=5, label = label)
    plt.plot(z, function_const, color =  color_opt, linestyle = '-',linewidth = 1, marker = None, label = None)
    plt.errorbar(x, y, yerr=y_error, xerr=None, ecolor = color_opt, elinewidth = 1, capsize = 5, linestyle='')
    plt.ylabel('Ratio/fit')
    plt.subplot(412)
    return function
    

def  main():
    """
        ----------------------------------------------
        Python script to create plots

        Plot includes 4 panels:
        1) Cross section vs. Bin
        2) Ratio plot including a fit (polynomial, third degree)
        3) Ratio compared to the fit
        4) Statistical uncertainites

        Figure is saved as eps under the name of the NNLO file                                                                      
        
        ----------------------------------------------                                                                                                                                  

        Example input:
        ./plot_5.py -f <path_to_file>/ALL.xm12_y1.dat -f  <path_to_file>/NLO.xm12_y1.dat -f  <path_to_file>/NNLO.xm12_y1.dat -r 0
                                                      
           \'-f\': required before every file
           \'-r\': File for the ratio plot has to be specified

        ----------------------------------------------   
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=main.__doc__ )
    parser.add_argument('-f', '--files', action = 'append', help = 'Enter files you want to analyze, \'-f\' required before every file')
    parser.add_argument('-r', '--ratio_to', help = 'Enter index of the file for the ratio plot (mandatory), starting from 0')
    args = parser.parse_args()
    
    filelist = args.files
    number = int(args.ratio_to)
    print('----------------------------------------------')
    print('File for the ratio plot: '+str(number))

    #Sets the plotting style
    #plotting_option = ('bo','rv','k*','g^','c+','ms')
    color = ('b','r','k','g','c','m')
    marker = ('o','v','*','^','+','s')

    #Cross section plot
    fig, axs = plt.subplots(4, 1, sharex=True)
    fig.subplots_adjust(hspace=0) # Remove horizontal space between axes

    plt.subplot(411)
    plt.yscale('log')
    for i in range(len(filelist)):
        print('Processing file ',str(i+1),'/',str(len(filelist)),': ',str(filelist[i]))
        crosssection(i, str(filelist[i]), color[i], marker[i])    
    plt.legend()   

    #Ratio plot
    plt.subplot(412)
    plt.yscale('linear')
    for i in range(len(filelist)):
        ratio(i, str(filelist[number]), str(filelist[i]), color[i], marker[i], color[number], marker[number], number)

    plt.subplot(414)
    plt.xlabel('Bin')
    
    #Save figure
    save_as = ''
    for i in range(len(filelist)):
        if(filelist[i].count('NNLO')==1):
            save_as = str(filelist[i])
    if(save_as==''):
        for i in range(len(filelist)):
            if(filelist[i].count('NLO')==1):
                save_as = str(filelist[i])
    while(save_as.count('/')>0):
        index = save_as.find('/')
        save_as = save_as[index + 1:]
    #save_as = 'CMS_plots07/'+save_as
    save_as = save_as.replace('.dat','.eps')
    plt.savefig(save_as)
    print('File saved as:')
    print(save_as)

    print('----------------------------------------------')


if __name__== "__main__":
  main()
