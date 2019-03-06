import numpy as N
import matplotlib.pyplot as P
import cmocean
from .rspec import rstokes
# Module with functions to visualize zdi data

def DynamicStokesI(filename):
    # Read measured data and difine sorted grid
    mdata = rstokes(filename)
    mphase = mdata[0]; mvr = mdata[1]; mstokesI = mdata[2]
    mdy=1./mphase.shape[0]
    mdx=mvr[1]-mvr[0]
    my, mx = N.mgrid[slice(0, 1 + mdy, mdy),
                    slice(mvr[0], mvr[-1] + mdx, mdx)]
    mindex_sort_phase=(mphase%1).argsort()
    msort_stokesI = mstokesI[mindex_sort_phase,:]
    # Read syntetic data and define sorted grid
    sdata = rstokes(filename[:-2]+'s1')
    sphase = sdata[0]; svr = sdata[1]; sstokesI = sdata[2]
    sdy=1./sphase.shape[0]
    sdx=svr[1]-svr[0]
    sy, sx = N.mgrid[slice(0, 1 + sdy, sdy),
                    slice(svr[0], svr[-1] + sdx, sdx)]
    sindex_sort_phase=(sphase%1).argsort()
    ssort_stokesI = sstokesI[sindex_sort_phase,:]
    # Plot dynamic spectra
    P.figure()
    ax1 = P.subplot(131)
    P.pcolormesh(mx, my, msort_stokesI, cmap = cmocean.cm.thermal)
    P.ylabel(r'Phase')
    P.xlabel(r'$v_r (km/s) $')
    ax2 = P.subplot(132, sharey = ax1)
    P.pcolormesh(sx, sy, ssort_stokesI, cmap = cmocean.cm.thermal)
    P.xlabel(r'$v_r (km/s) $')
    ax3 = P.subplot(133, sharey = ax1)
    P.pcolormesh(sx, sy, (msort_stokesI-ssort_stokesI)*(10.**4)/msort_stokesI, cmap = 'seismic')
    P.ylabel(r'Residual $(10^{-4})$')
    P.xlabel(r'$v_r (km/s) $')
    P.colorbar(format='%2d')
    P.setp(ax2.get_yticklabels(), visible=False)
    P.setp(ax3.get_yticklabels(), visible=False)
    P.tight_layout()
    P.show(False)
