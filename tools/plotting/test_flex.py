#! /usr/bin/env python2

from fastnlo import fastNLOLHAPDF

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import numpy as np

fnlo = fastNLOLHAPDF('InclusiveNJets_fnl2342b_I902309_v23_flex.tab')
fnlo.SetLHAPDFFilename('CT10nlo.LHgrid')
fnlo.SetLHAPDFMember(0)

mufs = np.arange(0.1, 1.5, 0.10)
murs = np.arange(0.1, 1.5, 0.10)
xs = np.zeros((mufs.size, murs.size))

for i, muf in enumerate(mufs):
    for j, mur in enumerate(murs):
        fnlo.SetScaleFactorsMuRMuF(mur, muf)
        fnlo.CalcCrossSection()
        xs[i][j] = np.sum(np.array(fnlo.GetCrossSection()))

fig = plt.figure(figsize=(13,13))
ax = fig.gca(projection='3d')
murs, mufs = np.meshgrid(murs, mufs)
ax.plot_surface(murs, mufs, xs, cmap=cm.jet, rstride=1, cstride=1, alpha=0.8)
ax.set_ylabel('Scale factor $\mu_F$')
ax.set_xlabel('Scale factor $\mu_R$')
ax.set_zlabel('Cross Section [pb/GeV]')
scatter1_proxy = matplotlib.lines.Line2D([0],[0], linestyle="none", c='black', marker = 'o')
ax.legend([scatter1_proxy], ['fastNLO flex table'], numpoints = 1)
ax.view_init(elev=15., azim=100.)
plt.show()
