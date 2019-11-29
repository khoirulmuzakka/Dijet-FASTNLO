#!/usr/bin/env python2
#-*- coding:utf-8 -*-

###########################################
#
# Plot alpha_s(M_Z) values for comparison
#
#
# Created by K. Rabbertz, 23.05.2019
#
###########################################
#
import argparse
import glob
import os
import re
import sys
# Use matplotlib with Cairo offline backend for eps, pdf, png, or svg output
import matplotlib as mpl
# mpl.use('Agg')
mpl.use('Cairo')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.ticker import (
    FormatStrFormatter, ScalarFormatter, AutoMinorLocator, MultipleLocator)
from matplotlib import cm
# numpy
import numpy as np

# Redefine ScalarFormatter


class ScalarFormatterForceFormat(ScalarFormatter):
    # Override function that finds format to use.
    def _set_format(self, vmin, vmax):
        self.format = "%1.2f"  # Give format here


# Some global definitions 0.1181, 11
_myred = '#aa062e'
_myblue = '#2e429a'
_asmz_label = r'$\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$'
_asmz_dlim = 0.0800
_text_dlim = 0.0810
_asmz_ulim = 0.1300
_ypos = [1, 3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16, 18, 20, 21, 22, 23]
_pdg = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
_cms = [0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
_atlas = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
_asmz = [0.1181,
         0.1135, 0.1144, 0.1139, 0.1151,
         0.1196, 0.1195, 0.1162, 0.1173,
         0.1127,
         0.1150, 0.1148, 0.1171, 0.1199,
         0.1185, 0.1164, 0.1192, 0.1185]
_dasmz_oth_up = [0,
                 18, 25, 23, 27,
                 22, 24, 21, 20,
                 36,
                 22, 23, 28, 15,
                 19, 29, 23, 35]
_dasmz_oth_dn = [0,
                 17, 25, 23, 26,
                 22, 24, 21, 20,
                 19,
                 22, 23, 28, 16,
                 26, 33, 19, 35]
_dasmz_scl_up = [11,
                 11, 16, 14, 9,
                 61, 60, 76, 63,
                 52,
                 50, 50, 69, 31,
                 22, 53, 24, 53]
_dasmz_scl_dn = [11,
                 5, 20, 1, 8,
                 13, 15, 61, 20,
                 19,
                 0, 50, 40, 19,
                 18, 28, 39, 24]

gs = gridspec.GridSpec(1, 1)
fig = plt.figure(figsize=(7, 9))
ax = plt.subplot(gs[0])
# ax.grid(False)
ax.grid(b=True, which='major', linestyle='-')
ax.grid(b=True, which='minor', linestyle=':')
ax.set_yticks([])
ax.set_xticks([0.11, 0.12, 0.13])
ax.set_xticks([0.108, 0.112, 0.114, 0.116, 0.118,
               0.122, 0.124, 0.126, 0.128], minor=True)
ax.tick_params(direction='in', which='both', labelsize='x-large')
ax.tick_params(axis='x', which='minor', length=4)
ax.tick_params(axis='x', which='major', length=8)
ax.set_title('%s results from            and         ' %
             _asmz_label, fontsize='xx-large')
ax.text(0.1070, 25.3, 'ATLAS', color=_myblue, fontsize='xx-large')
ax.text(0.1190, 25.3, 'CMS', color=_myred, fontsize='xx-large')
ax.set_xlim([_asmz_dlim, _asmz_ulim])
ax.set_ylim([0, _ypos[-1]+2])
ax.set_xlabel('%s' % _asmz_label, horizontalalignment='right',
              x=1.0, verticalalignment='top', y=1.0, fontsize='xx-large')
#xfmt = ScalarFormatterForceFormat()
# xfmt.set_powerlimits((0,0))
# ax.xaxis.set_major_formatter(xfmt)

#_ypos = np.array(_ypos)-0.5
ax.text(_text_dlim, 24, 'Inclusive jets, NLO',
        fontsize='large', fontweight='bold', va='center')
ax.text(_text_dlim,
        _ypos[-1], r'7 TeV, $p_\mathrm{T}$, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$', fontsize='large', va='center')
ax.text(_text_dlim,
        _ypos[-2], r'7 TeV, $p_\mathrm{T}$, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$+PDF', fontsize='large', va='center')
ax.text(_text_dlim,
        _ypos[-3], r'8 TeV, $p_\mathrm{T}$, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$', fontsize='large', va='center')
ax.text(_text_dlim,
        _ypos[-4], r'8 TeV, $p_\mathrm{T}$, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$+PDF', fontsize='large', va='center')
ax.text(_text_dlim, 19, 'Dijets, NLO', fontsize='large',
        fontweight='bold', va='center')
ax.text(_text_dlim,
        _ypos[-5], r'8 TeV, $\left<p_\mathrm{1,2}\right>$, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$+PDF', fontsize='large', va='center')
ax.text(_text_dlim, 17, '3-jets, NLO', fontsize='large',
        fontweight='bold', va='center')
ax.text(_text_dlim,
        _ypos[-6], r'7 TeV, $M_3$, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$', fontsize='large', va='center')
ax.text(_text_dlim,
        _ypos[-7], r'7 TeV, $R_{3/2}$, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$', fontsize='large', va='center')
ax.text(_text_dlim,
        _ypos[-8], r'8 TeV, $R_{3/2}$, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$, prel.', fontsize='large', va='center')
ax.text(_text_dlim,
        _ypos[-9], r'8 TeV, $R_{\Delta\phi}$, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$', fontsize='large', va='center')
ax.text(_text_dlim, 12, r'EEC, NLO',
        fontsize='large', fontweight='bold', va='center')
ax.text(_text_dlim,
        _ypos[-10], r'7 TeV, TEEC, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$', fontsize='large', va='center')
ax.text(_text_dlim,
        _ypos[-11], r'8 TeV, TEEC, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$', fontsize='large', va='center')
ax.text(_text_dlim,
        _ypos[-12], r'7 TeV, ATEEC, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$', fontsize='large', va='center')
ax.text(_text_dlim,
        _ypos[-13], r'8 TeV, ATEEC, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$', fontsize='large', va='center')
ax.text(_text_dlim,  7, r'top-antitop',
        fontsize='large', fontweight='bold', va='center')
ax.text(_text_dlim,
        _ypos[-14], r'7 TeV, $\sigma_\mathrm{tot}$, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$, NNLO+NNLL', fontsize='large', va='center')
ax.text(_text_dlim,
        _ypos[-15], r'13 TeV, $\sigma_\mathrm{tot}$, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$, NNLO', fontsize='large', va='center')
ax.text(_text_dlim,
        _ypos[-16], r'13 TeV, $\sigma_\mathrm{diff}$, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$, NLO, prel.', fontsize='large', va='center')
ax.text(_text_dlim,
        _ypos[-17], r'13 TeV, $\sigma_\mathrm{diff}$, $\alpha_\mathrm{s}(\mathrm{M}_\mathrm{Z})$+PDF, NLO, prel.', fontsize='large', va='center')
ax.text(_text_dlim, 2, '(Inner uncertainty: All except scale)',
        fontsize='small', va='center')
ax.text(_text_dlim, 1,
        'World average [PDG 2016]', fontsize='large', va='center')

xerr_oth = [_dasmz_oth_dn, _dasmz_oth_up]
xerr_scl = [_dasmz_scl_dn, _dasmz_scl_up]
xerr_oth = np.array(xerr_oth)/10000.
xerr_scl = np.array(xerr_scl)/10000.
xerr_tot = np.sqrt(xerr_oth*xerr_oth+xerr_scl*xerr_scl)

_pdg = np.array(_pdg)
_cms = np.array(_cms)
_atlas = np.array(_atlas)

ax.errorbar(_asmz*_pdg, _ypos*_pdg, xerr=xerr_oth*_pdg, marker='o', markersize=8,
            capsize=4, color='black', linestyle='none')
ax.errorbar(_asmz*_pdg, _ypos*_pdg, xerr=xerr_tot*_pdg, marker='', markersize=8,
            capsize=0, color='black', linestyle='none')
ax.errorbar(_asmz*_cms, _ypos*_cms, xerr=xerr_oth*_cms, marker='o', markersize=8,
            capsize=4, color=_myred, linestyle='none')
ax.errorbar(_asmz*_cms, _ypos*_cms, xerr=xerr_tot*_cms, marker='', markersize=8,
            capsize=0, color=_myred, linestyle='none')
ax.errorbar(_asmz*_atlas, _ypos*_atlas, xerr=xerr_oth*_atlas, marker='o', markersize=8,
            capsize=4, color=_myblue, linestyle='none')
ax.errorbar(_asmz*_atlas, _ypos*_atlas, xerr=xerr_tot*_atlas, marker='', markersize=8,
            capsize=0, color=_myblue, linestyle='none')
ax.axvspan(_asmz[0]-xerr_tot[0, 0], _asmz[0]+xerr_tot[1, 0], color='gold')

figname = 'asmz_lhc'
fig.savefig(figname+'.png')
fig.savefig(figname+'.pdf')

exit(0)
