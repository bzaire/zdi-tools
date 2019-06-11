import numpy as N
import matplotlib.pyplot as P
import matplotlib.colors
import cmocean
from .rspec import rstokes
# Module with functions to visualize zdi data

def DynamicStokesI_pcolor(filename):
    P.style.use('seaborn-notebook')
    # Read measured data and difine sorted grid
    mdata = rstokes(filename)
    mphase = mdata[0]; mvr = mdata[1]; msnI = mdata[2]; mstokesI = mdata[3]
    print(mstokesI.shape)
    mindex_sort_phase=(mphase%1).argsort()
    mx, my = N.meshgrid(mvr, (mphase%1)[mindex_sort_phase])
    msort_stokesI = mstokesI[mindex_sort_phase,:]
    # Read syntetic data and define sorted grid
    sdata = rstokes(filename[:-2]+'s1')
    sphase = sdata[0]; svr = sdata[1]; ssnI = sdata[2]; sstokesI = sdata[3]
    sdy=1./sphase.shape[0]
    sdx=svr[1]-svr[0]
    sindex_sort_phase=(sphase%1).argsort()
    sx, sy = N.meshgrid(svr, (sphase%1)[sindex_sort_phase])
    # Plot dynamic spectra
    P.figure()
    ax1 = P.subplot(131)
    cmap = cmocean.cm.balance
    P.pcolormesh(mx, my, mstokesI/mstokesI.mean(axis=0), cmap = cmap)
    P.ylabel(r'Phase')
    P.xlabel(r'$v_r (km/s) $')
    ax2 = P.subplot(132, sharey = ax1)
    P.pcolormesh(sx, sy, sstokesI/mstokesI.mean(axis=0), cmap = cmap)
    P.xlabel(r'$v_r (km/s) $')
    ax3 = P.subplot(133, sharey = ax1)
    P.pcolormesh(sx, sy, (mstokesI-sstokesI)*(10.**4), vmin=-30, vmax=30, cmap = cmap)
    P.ylabel(r'Residual $(10^{-4})$')
    P.xlabel(r'$v_r (km/s) $')
    P.colorbar(format='%2d')
    P.setp(ax2.get_yticklabels(), visible=False)
    P.setp(ax3.get_yticklabels(), visible=False)
    P.tight_layout()
    P.savefig('dynspec.png')
    P.show(False)

def DynamicStokesI_old(filename):
    P.style.use('seaborn-notebook')
    # Read measured data and difine sorted grid
    mdata = rstokes(filename)
    mphase_old = mdata[0]; mvr = mdata[1]; msnI = mdata[2]; mI_old = mdata[3]
    # index for phase folding
    mphase_fold = N.float32(mphase_old%1)
    # create a new grid
    dp = N.ones_like(mphase_fold)
    temp = mphase_fold.copy()
    temp.sort()
    for i in range(dp.shape[0]-1):
        dp[i] = temp[i+1] - temp[i]
    npoints = int(1./dp.min())
    mphase = N.linspace(0, 1, npoints, dtype=N.float32, endpoint=True)
#    for ip in mphase_fold:
#        if ip not in mphase:
#           mphase = N.concatenate((mphase,[ip]))
    mI = N.full((npoints, mI_old.shape[1]), N.nan) 
    for i in range(npoints):
        ip = mphase[i]
        jpf = N.bitwise_and((mphase_fold > (ip-dp.min()/2.)), (mphase_fold < (ip+dp.min()/2.)))
        if True in jpf:
            mI[i,:] = mI_old[jpf,:]  
    # sort in phase
    smp = mphase.argsort()
    mphase = mphase[smp]
    mI = mI[smp,:]
    ## Read zdi data and difine sorted grid
    cdata = rstokes(filename[:-2]+'s1')
    cphase_old = cdata[0]; cvr = cdata[1]; csnI = cdata[2]; cI_old = cdata[3]
    # index for phase folding
    cphase_fold = N.float32(cphase_old%1)
    # create a new grid
    cphase = N.linspace(0, 1, npoints, dtype=N.float32, endpoint=True)
    for ip in cphase_fold:
        if ip not in cphase:
           cphase = N.concatenate((cphase,[ip]))
    cI = N.full((npoints, cI_old.shape[1]), N.nan) 
    for i in range(npoints):
        ip = cphase[i]
        jpf = N.bitwise_and((cphase_fold > (ip-dp.min()/2.)), (cphase_fold < (ip+dp.min()/2.)))
        if True in jpf:
            cI[i,:] = cI_old[jpf,:]  
    # sort in phase
    cphase = cphase[smp]
    cI = cI[smp,:]
    # Plot dynamic spectra normalized by the mean line profile
    cmap = cmocean.cm.balance
    cmap.set_bad(color='gray', alpha=3./3)
    dpth = 0.004

    fig, ax = P.subplots(1,3, figsize=(7, 7), dpi=110)
    ax[0].ishow(mI[::-1,:]/mI_old.mean(axis=0), vmin=1.-dpth, vmax=1.+dpth, cmap = cmap, shape=mI.shape, aspect='auto', interpolation = 'none')
    ax[0].set_ylabel(r'Phase')
    ax[0].set_xlabel(r'$v_r (km/s) $')
    ax[0].set_xticks([0,mvr.shape[0]-1])
    ax[0].set_xticklabels([int(mvr[0]),int(mvr[-1])])
    ax[0].set_yticks([0,mphase.shape[0]-1])
    ax[0].set_yticklabels([mphase[-1], mphase[0]])
    ax[1].ishow(cI[::-1,:]/cI_old.mean(axis=0), vmin=1.-dpth, vmax=1.+dpth, shape=cI.shape, cmap = cmap,  aspect='auto', interpolation = 'none')
    imgplot = ax[2].ishow((mI-cI)[::-1,:]*(1.e4), vmin = -20, vmax = 20., cmap = cmap, shape=cI.shape, aspect='auto', interpolation = 'none')
    fig.colorbar(imgplot)
    ax[1].set_yticks([])
    ax[2].set_yticks([])
    ax[1].set_xlabel(r'$v_r (km/s) $')
    ax[1].set_xticks([0,mvr.shape[0]-1])
    ax[1].set_xticklabels([int(mvr[0]),int(mvr[-1])])
    ax[2].set_xlabel(r'$v_r (km/s) $')
    ax[2].set_xticks([0,mvr.shape[0]-1])
    ax[2].set_xticklabels([int(mvr[0]),int(mvr[-1])])
#    for iax in ax:
#        ones = N.ones_like(mphase)
#        iax.plot(ones*90., mphase, '--w')
#        iax.plot(ones*(-90.), mphase, '--w')
    P.tight_layout()
    P.savefig('dynspec.png')
    P.show()

def DynamicStokesI_sort(filename):
    P.style.use('seaborn-notebook')
    # Read measured data and difine sorted grid
    mdata = rstokes(filename)
    mphase_old = mdata[0]; mvr = mdata[1]; msnI = mdata[2]; mI_old = mdata[3]
    mphase = (mphase_old%1)
#    dp = N.ones_like(mphase_fold)
#    for i in range(dp.shape[0]-1):
#        dp[i] = mphase_fold[i+1] - mphase_fold[i]
#    npoints = int((mphase_old.max() - mphase_old.min())/dp.min())
#    # create a new grid
#    mphase = N.linspace(mphase_old.min(),mphase_old.max(),npoints,dtype=N.float64, endpoint=True)
#    for ip in mphase_fold:
#        if ip not in mphase:
#           mphase = N.concatenate((mphase,[ip]))
#    nmphase = mphase.shape[0]
#    mI = N.full((nmphase, mI_old.shape[1]), N.nan) 
#    for i in range(nmphase):
#        ip = mphase[i]
#        if ip in mphase_fold:
#            jpf = (mphase_fold==ip)
#            mI[i,:] = mI_old[jpf,:]  
#    # sort in phase
    smp = mphase.argsort()
    mphase = mphase[smp]
    mI = mI_old[smp,:]
    ## Read zdi data and difine sorted grid
    cdata = rstokes(filename[:-2]+'s1')
    cphase_old = cdata[0]; cvr = cdata[1]; csnI = cdata[2]; cI_old = cdata[3]
    cphase = (cphase_old%1)
#    # index for phase folding
#    cphase_fold = N.float64(cphase_old)
#    # create a new grid
#    cphase = N.linspace(cphase_old.min(),cphase_old.max(),npoints,dtype=N.float64, endpoint=True)
#    for ip in cphase_fold:
#        if ip not in cphase:
#           cphase = N.concatenate((cphase,[ip]))
#    ncphase = cphase.shape[0]
#    cI = N.full((ncphase, cI_old.shape[1]), N.nan) 
#    for i in range(ncphase):
#        ip = cphase[i]
#        if ip in cphase_fold:
#            jpf = (cphase_fold==ip)
#            cI[i,:] = cI_old[jpf,:]  
#    # sort in phase
    scp = cphase.argsort()
    cphase = cphase[scp]
    cI = cI_old[scp,:]
    # Plot dynamic spectra normalized by the mean line profile
    cmap = cmocean.cm.balance
    cmap.set_bad(color='gray', alpha=3./3)
    dpth = 0.004

    fig, ax = P.subplots(1,3, figsize=(7, 7), dpi=110)
    ax[0].matshow(mI[::-1,:]/mI_old.mean(axis=0), vmin=1.-dpth, vmax=1.+dpth, cmap = cmap, shape=mI_old.shape, aspect='auto', interpolation = 'none')
    #ax[0].imshow(mI[::-1,:], vmin=1.-dpth, vmax=1.+dpth, cmap = cmap, shape=mI.shape, aspect='auto', interpolation = 'none')
    ax[0].set_ylabel(r'Phase')
    ax[0].set_xlabel(r'$v_r (km/s) $')
    ax[0].set_xticks([0,mvr.shape[0]-1])
    ax[0].set_xticklabels([int(mvr[0]),int(mvr[-1])])
    ax[0].set_yticks([0,mphase.shape[0]-1])
    ax[0].set_yticklabels([mphase[-1], mphase[0]])
    ax[1].matshow(cI[::-1,:]/cI_old.mean(axis=0), vmin=1.-dpth, vmax=1.+dpth, shape=cI_old.shape, cmap = cmap,  aspect='auto', interpolation = 'none')
    #ax[1].imshow(cI[::-1,:], vmin=1.-dpth, vmax=1.+dpth, shape=cI.shape, cmap = cmap,  aspect='auto', interpolation = 'none')
    imgplot = ax[2].matshow((mI-cI)[::-1,:]*(1.e4), vmin = -20, vmax = 20., cmap = cmap, shape=cI_old.shape, aspect='auto', interpolation = 'none')
    fig.colorbar(imgplot)
    ax[1].set_yticks([])
    ax[2].set_yticks([])
    ax[1].set_xlabel(r'$v_r (km/s) $')
    ax[1].set_xticks([0,mvr.shape[0]-1])
    ax[1].set_xticklabels([int(mvr[0]),int(mvr[-1])])
    ax[2].set_xlabel(r'$v_r (km/s) $')
    ax[2].set_xticks([0,mvr.shape[0]-1])
    ax[2].set_xticklabels([int(mvr[0]),int(mvr[-1])])
#    for iax in ax:
#        ones = N.ones_like(mphase)
#        iax.plot(ones*90., mphase, '--w')
#        iax.plot(ones*(-90.), mphase, '--w')
    P.tight_layout()
    P.savefig('dynspec.png')
    P.show()

def DynamicStokesI(filename):
    P.style.use('seaborn-notebook')
    # Read measured data and difine sorted grid
    mdata = rstokes(filename)
    mphase_old = mdata[0]; mvr = mdata[1]; msnI = mdata[2]; mI_old = mdata[3]
#    mphase_fold = N.float64(mphase_old)
#    dp = N.ones_like(mphase_fold)
#    for i in range(dp.shape[0]-1):
#        dp[i] = mphase_fold[i+1] - mphase_fold[i]
#    npoints = int((mphase_old.max() - mphase_old.min())/dp.min())
#    # create a new grid
#    mphase = N.linspace(mphase_old.min(),mphase_old.max(),npoints,dtype=N.float64, endpoint=True)
#    for ip in mphase_fold:
#        if ip not in mphase:
#           mphase = N.concatenate((mphase,[ip]))
#    nmphase = mphase.shape[0]
#    mI = N.full((nmphase, mI_old.shape[1]), N.nan) 
#    for i in range(nmphase):
#        ip = mphase[i]
#        if ip in mphase_fold:
#            jpf = (mphase_fold==ip)
#            mI[i,:] = mI_old[jpf,:]  
#    # sort in phase
#    smp = mphase.argsort()
#    mphase = mphase[smp]
#    mI = mI[smp,:]
    ## Read zdi data and difine sorted grid
    cdata = rstokes(filename[:-2]+'s1')
    cphase_old = cdata[0]; cvr = cdata[1]; csnI = cdata[2]; cI_old = cdata[3]
#    # index for phase folding
#    cphase_fold = N.float64(cphase_old)
#    # create a new grid
#    cphase = N.linspace(cphase_old.min(),cphase_old.max(),npoints,dtype=N.float64, endpoint=True)
#    for ip in cphase_fold:
#        if ip not in cphase:
#           cphase = N.concatenate((cphase,[ip]))
#    ncphase = cphase.shape[0]
#    cI = N.full((ncphase, cI_old.shape[1]), N.nan) 
#    for i in range(ncphase):
#        ip = cphase[i]
#        if ip in cphase_fold:
#            jpf = (cphase_fold==ip)
#            cI[i,:] = cI_old[jpf,:]  
#    # sort in phase
#    scp = cphase.argsort()
#    cphase = cphase[scp]
#    cI = cI[scp,:]
    # Plot dynamic spectra normalized by the mean line profile
    cmap = cmocean.cm.balance
    cmap.set_bad(color='gray', alpha=3./3)
    dpth = 0.004

    fig, ax = P.subplots(1,3, figsize=(7, 7), dpi=110)
    ax[0].matshow(mI_old[::-1,:]/mI_old.mean(axis=0), vmin=1.-dpth, vmax=1.+dpth, cmap = cmap, shape=mI_old.shape, aspect='auto', interpolation = 'none')
    #ax[0].imshow(mI[::-1,:], vmin=1.-dpth, vmax=1.+dpth, cmap = cmap, shape=mI.shape, aspect='auto', interpolation = 'none')
    ax[0].set_ylabel(r'Phase')
    ax[0].set_xlabel(r'$v_r (km/s) $')
    ax[0].set_xticks([0,mvr.shape[0]-1])
    ax[0].set_xticklabels([int(mvr[0]),int(mvr[-1])])
    ax[0].set_yticks([0,mphase_old.shape[0]-1])
    ax[0].set_yticklabels([mphase_old[-1], mphase_old[0]])
    ax[1].matshow(cI_old[::-1,:]/cI_old.mean(axis=0), vmin=1.-dpth, vmax=1.+dpth, shape=cI_old.shape, cmap = cmap,  aspect='auto', interpolation = 'none')
    #ax[1].imshow(cI[::-1,:], vmin=1.-dpth, vmax=1.+dpth, shape=cI.shape, cmap = cmap,  aspect='auto', interpolation = 'none')
    imgplot = ax[2].matshow((mI_old-cI_old)[::-1,:]*(1.e4), vmin = -20, vmax = 20., cmap = cmap, shape=cI_old.shape, aspect='auto', interpolation = 'none')
    fig.colorbar(imgplot)
    ax[1].set_yticks([])
    ax[2].set_yticks([])
    ax[1].set_xlabel(r'$v_r (km/s) $')
    ax[1].set_xticks([0,mvr.shape[0]-1])
    ax[1].set_xticklabels([int(mvr[0]),int(mvr[-1])])
    ax[2].set_xlabel(r'$v_r (km/s) $')
    ax[2].set_xticks([0,mvr.shape[0]-1])
    ax[2].set_xticklabels([int(mvr[0]),int(mvr[-1])])
#    for iax in ax:
#        ones = N.ones_like(mphase)
#        iax.plot(ones*90., mphase, '--w')
#        iax.plot(ones*(-90.), mphase, '--w')
    P.tight_layout()
    P.savefig('dynspec.png')
    P.show()

def DynamicSpec(phase_old, vr, I_old):
    P.style.use('seaborn-notebook')
    # index for phase folding
    phase_fold = N.float64(phase_old%1)
    # create a new grid
#    phase = N.linspace(0,1,int(2*phase_old.shape[0]),dtype=N.float64, endpoint=True)
    phase = N.linspace(phase_old.min(),phase_old.max(),10000,dtype=N.float32, endpoint=True)
    for ip in phase_fold:
        if ip not in phase:
           phase = N.concatenate((phase,[ip]))
    nphase = phase.shape[0]
    I = N.full((nphase,I_old.shape[1]), N.nan) 
    for i in range(nphase):
        ip = phase[i]
        if ip in phase_fold:
            jpf = (phase_fold==ip)
            I[i,:] = I_old[jpf,:]  
    # sort in phase
    sp = phase.argsort()
    phase = phase[sp]
    # Plot dynamic spectra normalized by the mean line profile
    fig, ax = P.subplots(1,1, figsize=(4, 6), dpi=110)
    z = I[sp,:]/I_old.mean(axis=0)
    cmap = cmocean.cm.balance
    cmap.set_bad(color='gray', alpha=3./3)
    dpth = 0.004
    imgplot = ax.imshow(z[::-1,:], vmin=1-dpth, vmax=1.+dpth, shape=z.shape, cmap = cmap, aspect='auto', interpolation = 'none')
    #imgplot.cmap.set_under(color="gray")
    fig.colorbar(imgplot)
    ax.set_ylabel(r'Phase')
    ax.set_xlabel(r'$v_r (km/s) $')
    ax.set_xticks([0,vr.shape[0]-1])
    ax.set_xticklabels([int(vr[0]),int(vr[-1])])
    ax.set_yticks([0,phase.shape[0]-1])
    ax.set_yticklabels([phase[-1], phase[0]])
    P.tight_layout()
    P.savefig('dynspec.png')
    P.show()
