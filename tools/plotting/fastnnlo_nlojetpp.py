#!/usr/bin/env python2
#-*- coding:utf-8 -*-
import fastnlo
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

# Default arguments
proc = '1jet'
jobn = 'LO-CMS7'
kinn = 'vBa'
obsv = 'fnl2332d_xptji_y1'
fnlo = 1
# Get from cmdline argument
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
if len(sys.argv) > 1:
    proc = sys.argv[1]
if len(sys.argv) > 2:
    jobn = sys.argv[2]
if len(sys.argv) > 3:
    kinn = sys.argv[3]
if len(sys.argv) > 4:
    obsv = sys.argv[4]
if len(sys.argv) > 5:
    fnlo = sys.argv[5]

# Extract order/contribution from job type (substring before first '-')
order  = jobn.split('-')[0]
exper  = jobn.split('-')[1]
ordcol = 6
if order == 'NLO':
    ordcol = 7
elif order == 'NNLO':
    ordcol = 8

# Prepare binning arrays
xl = []      # left bin border
xm = []      # bin "center"
xu = []      # right bin border

# Read binning and cross sections from NNLOJET dat file
if kinn == '_':
    datfile = 'Combined/Final/'+order+'.'+obsv+'.dat'
else:
    datfile = 'Combined/Final/'+order+'.'+kinn+'.'+obsv+'.dat'
xs_scls = []
dxs_scls = []
print 'Reading from NNLOJET dat file ', datfile
xl.append(np.loadtxt(datfile,usecols=(0,)))
xu.append(np.loadtxt(datfile,usecols=(2,)))
xl = np.array(xl)
xu = np.array(xu)
# Determine no. of observable bins
nobs = xl.size
print 'Number of observable bins: ', nobs
xs_scls.append(np.loadtxt(datfile,usecols=(3,)))
dxs_scls.append(np.loadtxt(datfile,usecols=(4,)))
xs_nnlo   = np.array(xs_scls)/1000.  # Conversion of fb to pb
dxs_nnlo  = np.array(dxs_scls)/1000.
xs_tmp1 = np.array(xs_nnlo)
xs_tmp2 = np.array(dxs_nnlo)
xs1 = xs_tmp1.reshape(xs_tmp1.size)
xs2 = xs_tmp2.reshape(xs_tmp2.size)
xs_nnlo = xs1[0:nobs]
dxs_nnlo = xs2[0:nobs]
dst_nnlo = np.divide(dxs_nnlo, xs_nnlo, out=np.ones_like(dxs_nnlo), where=xs_nnlo!=0)
#print 'xs_nnlo', xs_nnlo
print 'dst_nnlo', dst_nnlo

# Evaluate cross sections from NNLOJET fastNLO tables
# INFO=0, WARNING=1
SetGlobalVerbosity(1)
xs_fnnl1 = []
xs_fnnl2 = []
nnlotab = 'Combined/Final/'+proc+'.'+order+'-'+exper+'.'+obsv+'.tab.gz'
print 'NNLOJET fastNLO table is ', nnlotab
fnlo = fastNLOLHAPDF(nnlotab)
fnlo.SetLHAPDFFilename('CT14nnlo')
fnlo.SetLHAPDFMember(0)
fnlo.SetMuFFunctionalForm(0); # kScale1=0
fnlo.SetMuRFunctionalForm(0); # kScale2=1
fnlo.SetContributionON(fastnlo.kFixedOrder,0,True);  #Switch LO on     
if 'NLO' in order:
    fnlo.SetContributionON(fastnlo.kFixedOrder,1,True); #Switch NLO on
else:
    fnlo.SetContributionON(fastnlo.kFixedOrder,1,False); #Switch NLO off
if 'NNLO' in order:
    fnlo.SetContributionON(fastnlo.kFixedOrder,2,True); #Switch NNLO on
else:
    fnlo.SetContributionON(fastnlo.kFixedOrder,2,False); #Switch NNLO off
fnlo.CalcCrossSection()
xs_fnnl1.append(fnlo.GetCrossSection())
fnlo.SetMuFFunctionalForm(1); # kScale1=0
fnlo.SetMuRFunctionalForm(1); # kScale2=1
fnlo.CalcCrossSection()
xs_fnnl2.append(fnlo.GetCrossSection())
xs_tmp1 = np.array(xs_fnnl1)
xs_tmp2 = np.array(xs_fnnl2)
xs1 = xs_tmp1.reshape(xs_tmp1.size)
xs2 = xs_tmp2.reshape(xs_tmp2.size)
xs_fnnl1 = xs1[0:nobs]
xs_fnnl2 = xs2[0:nobs]

# Evaluate cross sections from NLOJet++ fastNLO tables
xs_fnl1 = []
xs_fnl2 = []
nlotab1 = 'NLOJet++/fnl2332d_I1208923.tab.gz'
nlotab2 = 'NLOJet++/fnl2332dptmax.tab.gz'
print 'NLOJet++ fastNLO tables are ', nlotab1, ' and ', nlotab2
fnlo = fastNLOLHAPDF(nlotab1)
fnlo.SetLHAPDFFilename('CT14nnlo')
fnlo.SetLHAPDFMember(0)
fnlo.SetContributionON(fastnlo.kFixedOrder,0,True);  #Switch LO on     
if 'NLO' in order:
    fnlo.SetContributionON(fastnlo.kFixedOrder,1,True); #Switch NLO on
else:
    fnlo.SetContributionON(fastnlo.kFixedOrder,1,False); #Switch NLO off
if 'NNLO' in order:
    print 'WARNING! NNLO not available from NLOJet++, ignored.'
fnlo.CalcCrossSection()
xs_fnl1.append(fnlo.GetCrossSection())
fnlo = fastNLOLHAPDF(nlotab2)
fnlo.SetLHAPDFFilename('CT14nnlo')
fnlo.SetLHAPDFMember(0)
fnlo.SetContributionON(fastnlo.kFixedOrder,0,True);  #Switch LO on     
if 'NLO' in order:
    fnlo.SetContributionON(fastnlo.kFixedOrder,1,True); #Switch NLO on
else:
    fnlo.SetContributionON(fastnlo.kFixedOrder,1,False); #Switch NLO off
if 'NNLO' in order:
    print 'WARNING! NNLO not available from NLOJet++, ignored.'
fnlo.CalcCrossSection()
xs_fnl2.append(fnlo.GetCrossSection())
xs_tmp1 = np.array(xs_fnl1)
xs_tmp2 = np.array(xs_fnl2)
xs1 = xs_tmp1.reshape(xs_tmp1.size)
xs2 = xs_tmp2.reshape(xs_tmp2.size)
xs_fnl1 = xs1[0:nobs]
xs_fnl2 = xs2[0:nobs]

# Read statistical uncertainties from NLOJet++ fastNLO log files
if 'NLO' in order:
    nlolog1 = 'NLOJet++/fnl2332d-hhc-nlo-2jet_statunc_ct14nnlo.log'
    nlolog2 = 'NLOJet++/fnl2332dptmax-hhc-nlo-2jet_statunc_ct14nnlo.log'
else:
    nlolog1 = 'NLOJet++/fnl2332d-hhc-born-2jet_statunc_ct14nnlo.log'
    nlolog2 = 'NLOJet++/fnl2332dptmax-hhc-born-2jet_statunc_ct14nnlo.log'
print 'NLOJet++ fastNLO statlogs are ', nlolog1, ' and ', nlolog2
xs_tmp1 = np.loadtxt(nlolog1,usecols=(3,),comments=['#',' #','C','L'])
xs_tmp2 = np.loadtxt(nlolog2,usecols=(3,),comments=['#',' #','C','L'])
xs1 = xs_tmp1.reshape(xs_tmp1.size)
xs2 = xs_tmp2.reshape(xs_tmp2.size)
dst_fnl1 = xs1[0:nobs]
dst_fnl2 = xs2[0:nobs]
print "dst_fnl1", dst_fnl1
print "dst_fnl2", dst_fnl2

# Calculate ratios with statistical uncertainties
r_fnnl1  = np.divide(xs_fnnl1, xs_nnlo, out=np.ones_like(xs_fnnl1), where=xs_nnlo!=0)
r_fnnl2  = np.divide(xs_fnnl2, xs_nnlo, out=np.ones_like(xs_fnnl2), where=xs_nnlo!=0)
dr_fnnl1 = np.multiply(r_fnnl1,dst_nnlo)
dr_fnnl2 = np.multiply(r_fnnl2,dst_nnlo)
r_fnl1   = np.divide(xs_fnl1, xs_nnlo, out=np.ones_like(xs_fnl1), where=xs_nnlo!=0)
r_fnl2   = np.divide(xs_fnl2, xs_nnlo, out=np.ones_like(xs_fnl2), where=xs_nnlo!=0)
dst_fnl1 = np.square(dst_fnl1) + np.square(dst_nnlo)
dst_fnl2 = np.square(dst_fnl2) + np.square(dst_nnlo)
dst_fnl1 = np.sqrt(dst_fnl1)
dst_fnl2 = np.sqrt(dst_fnl2)
dr_fnl1  = np.multiply(r_fnl1,dst_fnl1)
dr_fnl2  = np.multiply(r_fnl2,dst_fnl2)

# Plot
titwgt = 'bold'
limfs = 'x-large'
fig = plt.figure()
ax  = fig.gca()
plt.title(r'Ratios to NNLOJET: {} {} {}'.format(proc, jobn, obsv), fontweight=titwgt)
plt.xlabel('Observable bin index', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0)
plt.ylabel('Ratio to NNLOJET', horizontalalignment='right', x=1.0, verticalalignment='top', y=1.0, labelpad=20)
plt.axhline(y=1.001, linestyle='--', linewidth=1.0, color='black')
plt.axhline(y=0.999, linestyle='--', linewidth=1.0, color='black')
plt.fill_between([0.0,34.0],0.999,1.001, color='black', alpha=0.1)
plt.text(34.6,1.0,u'$\pm$1‰',fontsize=limfs)

dx = 1./4.
x  = np.arange(1   , nobs+1.e-6)
xa = np.arange(1-dx, nobs-dx+1.e-6)
xb = np.arange(1+dx, nobs+dx+1.e-6)

# Don't show statistical uncertainties for NNLOJET grid/NNLOJET with identical events
ratio1 = plt.errorbar(x, r_fnnl1, yerr=0.*dr_fnnl1, marker='s', linestyle='none', label=r'NNLOJET grid, $\mu_r=\mu_f=p_{\rm T,jet}$', color='blue')
ratio2 = plt.errorbar(x, r_fnnl2, yerr=0.*dr_fnnl2, marker='^', linestyle='none', label=r'NNLOJET grid, $\mu_r=\mu_f=p_{\rm T,max}$', color='orange')
# Show quadratically added statistical uncertainties for NLOJet++ / NNLOJET
ratio3 = plt.errorbar(x, r_fnl1, yerr=dr_fnl1, marker='o', linestyle='none', label=r'NLOJet++ grid, $\mu_r=\mu_f=p_{\rm T,jet}$', color='magenta')
ratio4 = plt.errorbar(x, r_fnl2, yerr=dr_fnl2, marker='v', linestyle='none', label=r'NLOJet++ grid, $\mu_r=\mu_f=p_{\rm T,max}$', color='green')

plt.xlim(0.0,34.0)
if 'NLO' in order:
    plt.ylim(0.95,1.1)
else:
    plt.ylim(0.99,1.01)

handles = [ratio1,ratio2,ratio3,ratio4]
labels  = [h.get_label() for h in handles]

legend = ax.legend(handles, labels, loc='upper left', numpoints=1, handlelength=0)
legend = ax.legend(handles, labels, title=r'APPLfast / NNLOJET', loc='upper left', numpoints=1, handlelength=0)
legend.get_title().set_fontsize(limfs)

fignam = proc+'.'+jobn+'.'+obsv+'.'+'nlocomp'+'.png'
plt.savefig(fignam)

exit(0)
