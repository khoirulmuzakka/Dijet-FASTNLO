#!/usr/bin/env python2
#-*- coding:utf-8 -*-
from fastnlo import fastNLOLHAPDF
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def main():

    fnlo = fastNLOLHAPDF('fnl5332gy0_v23_fix_I1459051.tab.gz')
    fnlo.SetLHAPDFFilename('CT10nlo')
    fnlo.SetLHAPDFMember(0)

    # Enable Hoppet muf scale variations
    fnlo.UseHoppetScaleVariations(True)
    # Generate an array with values from 0.8 to 1.2 with a stepsize of 0.02
    hop_muf_vars = []
    hop_xs_vars = []
    # Loop over that array
    for muf in np.arange(0.25, 2.5, 0.1):
        fnlo.SetScaleFactorsMuRMuF(1.0, muf)
        fnlo.CalcCrossSection()
        # Just use the absolute summed cross section for simplicity
        hop_muf_vars.append(muf)
        hop_xs_vars.append(np.sum(np.array(fnlo.GetCrossSection())))

    # Do the same without hoppet
    fnlo.UseHoppetScaleVariations(False)
    def_muf_vars = []
    def_xs_vars  = []

    for muf in [0.5, 1.0, 2.0]:
        fnlo.SetScaleFactorsMuRMuF(1.0, muf)
        fnlo.CalcCrossSection()
        def_muf_vars.append(muf)
        def_xs_vars.append(np.sum(np.array(fnlo.GetCrossSection())))


    # Produce a very simple plot
    # Create a figure (canvas)
    fig = plt.figure()
    # Add a single subplot axis
    ax = fig.add_subplot(111)

    # plot connects the dots with a line
    ax.plot(hop_muf_vars, hop_xs_vars, label='Hoppet scale var.', color='green')
    # scatter draws the individual points
    ax.scatter(def_muf_vars, def_xs_vars, label='Default scale var.', color='blue')


    ax.set_ylabel('Cross Section $\sigma$')
    ax.set_xlabel('Scale factor $\mu_F$')

    # Add a legend at upper right location.
    ax.legend(loc='upper right')

    # Save to file, extension defines backend.
    fig.savefig('scale.png')

# If called as a script directly __name__ is main and will execute the main function
# if imported into another script then main function is not automatically called.
if __name__ == '__main__':
    main()
