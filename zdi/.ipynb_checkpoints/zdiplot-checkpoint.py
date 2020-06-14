import numpy as N
import matplotlib.pyplot as P
import matplotlib.colors
import matplotlib
import cmocean
from .rspec import rstokes
# Module with functions to visualize zdi data
P.style.use(['seaborn-white', 'seaborn-paper', 'seaborn-ticks'])
matplotlib.rc("font", family="Times New Roman", size=20)
matplotlib.rc('xtick', labelsize='medium')
matplotlib.rc('ytick', labelsize='medium')


def surf_brightness(vr,I):
    from scipy import interpolate
    flat_vr = list(i for a in vr for i in a)
    flat_I = list(i for a in I for i in a)
    f = interpolate.interp1d(flat_vr, flat_I)
    return f

def interp_data(x, y, z):
    yvals = (y%1)
    iy = (y%1).argsort()
    vmin = min([min(i) for i in x])
    vmax = max([max(i) for i in x])
    xvals = N.arange(vmin, vmax, 1.8)
    X, Y = N.meshgrid(xvals, yvals[iy])
    Z = N.zeros_like(X)
    for i in range(y.shape[0]):
        Z[i, :] = N.interp(xvals, x[iy[i]], z[iy[i]])
    return Z

def interp_data1(x, y, z):
    vmin = min([min(i) for i in x])
    vmax = max([max(i) for i in x])
    xvals = N.arange(vmin, vmax, 1.8)
    X, Y = N.meshgrid(xvals, y)
    Z = N.zeros_like(X)
    for i in range(y.shape[0]):
        Z[i, :] = N.interp(xvals, x[i], z[i])
    return xvals, Z

def interp_data2(x, y, z):
    from scipy.interpolate import RegularGridInterpolator
    X, Y = N.meshgrid(x, y)
    Z = interp_data(x, y, z)
    my_interpolating_function = RegularGridInterpolator((X.flatten(), Y.flatten()), Z.flatten())
    vmin = min([min(i) for i in x])
    vmax = max([max(i) for i in x])
    xvals = N.arange(vmin, vmax, 1.8)
    tyvals = N.arange(0,1,0.03)
    yvals = N.unique(N.concatenate((tyvals, y)))
    nX, nY = N.meshgrid(xvals, yvals)
    nZ = my_interpolating_function(xvals, yvals)
    return nX, nY, N.ma.where(nY not in y, nZ, N.nan)


def extend_y(y):
    dy = N.array([y[i+1]-y[i] for i in range(y.shape[0]-1)])
    igaps = N.where(dy > 0.01)[0]
    new_y = y.copy()
    count = 0
    for j in igaps:
        #extra_y = N.arange(y[j], y[j+1], N.median(dy)*2) 
        extra_y = N.arange(y[j], y[j+1], N.min(dy)) 
        new_y = N.concatenate((new_y, extra_y)) 
        count += 1
    return new_y 

def refine(y, z):
    new_y = extend_y(y) 
    iy = new_y.argsort()
    new_z = N.zeros((new_y.shape[0], z.shape[1]))
    count = 0
    for i in range(new_y.shape[0]):
        if new_y[iy[i]] in y:
            ipos = N.where(new_y[iy[i]] == y)[0]
            new_z[i,:] = z[ipos,:] 
            count += 1
        else:
            new_z[i,:] = N.zeros(z.shape[1])
    return new_y[iy], new_z

def SpecV(filename, iSort=False, title='',save=''):
    # Read data
    cycleI, vrI, snI, I, cycleV, vrV, snV, V =  rstokes(filename)
    vmin = min([min(i) for i in vrV])
    vmax = max([max(i) for i in vrV])
    new_vr = N.arange(vmin, vmax, 1.8)

    # Plot dynamic spectra
    cmap = cmocean.cm.gray_r  # Colormap used in plot
    P.figure(figsize=(4,5))
    ax1 = P.subplot(111)
    pc = 0.0005
    if iSort:
        tZ = interp_data(vrV, cycleV, V)
        phaseV = (cycleV%1)
        iPhase = (cycleV%1).argsort()
        new_phase, Z = refine(phaseV[iPhase], tZ)
        X, Y = N.meshgrid(new_vr, new_phase)
        masked_Z = N.ma.array(Z, mask = Z == 0) 
        P.pcolormesh(X, Y, masked_Z, vmin=-pc, vmax=pc, cmap=cmap)
        P.ylabel(r'Rotation phase', fontsize='large')
        P.title(r'%s' %title, fontsize='large')
        P.ylim((0,1))
    else:
        tvr, tZ = interp_data1(vrV, cycleV, V)
        new_cycle, Z = refine(cycleV, tZ)
        X, Y = N.meshgrid(new_vr, new_cycle)
        masked_Z = N.ma.array(Z, mask = Z == 0)
        P.pcolormesh(X, Y, masked_Z, vmin=-pc, vmax=pc, cmap=cmap)
        P.ylabel(r'Rotation cycle')
    P.xlabel(r'$v_\mathrm{r} (\mathrm{km/s}) $', fontsize='large')
    P.xticks(fontsize='medium')
    P.yticks(fontsize='medium')
    P.tight_layout()
    cmap = P.colorbar(extend='both', fraction=0.15, shrink=0.5)
    cmap.ax.tick_params(labelsize='medium')
    if save != '':
        P.savefig(save+'.pdf')
    P.show()

def SpecI(filename, iSort=False, title='',save=''):
    # Read data
    cycleI, vrI, snI, I, cycleV, vrV, snV, V =  rstokes(filename)
    vmin = min([min(i) for i in vrI])
    vmax = max([max(i) for i in vrI])
    new_vr = N.arange(vmin, vmax, 1.8)

    # Plot dynamic spectra
    cmap = cmocean.cm.gray  # Colormap used in plot
    P.figure(figsize=(4,5))
    ax1 = P.subplot(111)
    pc = 0.003
    if iSort:
        tZ = interp_data(vrI, cycleI, I)
        phaseI = (cycleI%1)
        iPhase = (cycleI%1).argsort()
        new_phase, Z = refine(phaseI[iPhase], tZ)
        X, Y = N.meshgrid(new_vr, new_phase)
        rZ = N.zeros_like(X)
        meanZ = N.ma.array(Z, mask = Z == 0).mean(axis=0)
        for j in range(new_phase.shape[0]):
            rZ[j, :] = Z[j,:]/meanZ
        masked_Z = N.ma.array(rZ, mask = Z == 0) 
        P.pcolormesh(X, Y, masked_Z, vmin=1.-pc, vmax=1.+pc, cmap=cmap)
        P.ylabel(r'Rotation phase', fontsize='large')
        P.title(r'%s' %title, fontsize='large')
        P.ylim((0,1))
    else:
        tvr, tZ = interp_data1(vrI, cycleI, I)
        new_cycle, Z = refine(cycleI, tZ)
        X, Y = N.meshgrid(new_vr, new_cycle)
        rZ = N.zeros_like(X)
        rtZ = N.zeros_like(tZ)
        meanZ = N.ma.array(Z, mask = Z == 0).mean(axis=0)
        meantZ = tZ.mean(axis=0)
        for j in range(new_cycle.shape[0]):
            rZ[j, :] = Z[j,:]/meanZ
        for j in range(cycleI.shape[0]):
            rtZ[j, :] = tZ[j,:]/meantZ
        masked_Z = N.ma.array(rZ, mask = Z == 0)
        P.pcolormesh(X, Y, masked_Z, vmin=1.-pc, vmax=1.+pc, cmap=cmap)
#        P.imshow(rtZ, vmin=1.-pc, vmax=1.+pc, cmap=cmap, origin='lower', extent=[tvr.min(), tvr.max(), cycleI.min(), cycleI.max()], aspect=30.)
#        P.yticks(cycleI)
        P.ylabel(r'Rotation cycle')
    P.xlabel(r'$v_\mathrm{r} (\mathrm{km/s}) $', fontsize='large')
    P.xticks(fontsize='medium')
    P.yticks(fontsize='medium')
    P.tight_layout()
    cmap = P.colorbar(extend='both', fraction=0.15, shrink=0.5)
    cmap.ax.tick_params(labelsize='medium')
    if save != '':
        P.savefig(save+'.pdf')
    P.show()

def SpecHa(filename, iSort=False,  title='',save=''):
    # Read data
    cycleI, vrI, snI, I, cycleV, vrV, snV, V =  rstokes(filename)
    vmin = min([min(i) for i in vrI])
    vmax = max([max(i) for i in vrI])
    new_vr = N.arange(vmin, vmax, 1.8)

    # Plot dynamic spectra
    cmap = cmocean.cm.matter_r  # Colormap used in plot
    P.figure(figsize=(4,5))
    ax1 = P.subplot(111)
    pc = 0.3
    if iSort:
        tZ = interp_data(vrI, cycleI, I)
        N.save('hamatriz',tZ)
        phaseI = (cycleI%1)
        iPhase = (cycleI%1).argsort()
        new_phase, Z = refine(phaseI[iPhase], tZ)
        X, Y = N.meshgrid(new_vr, new_phase)
        masked_Z = N.ma.array(Z, mask = Z == 0) 
        P.pcolormesh(X, Y, masked_Z, vmin=1.-pc, vmax=1.+pc, cmap=cmap)
        P.ylabel(r'Rotation phase', fontsize='large')
    else:
        tvr, tZ = interp_data1(vrI, cycleI, I)
        new_cycle, Z = refine(cycleI, tZ)
        X, Y = N.meshgrid(new_vr, new_cycle)
        masked_Z = N.ma.array(Z, mask = Z == 0)
        P.pcolormesh(X, Y, masked_Z, vmin=1.-pc, vmax=1.+pc, cmap=cmap)
        P.ylabel(r'Rotation cycle', fontsize='large')
    P.xlim((-400.,400.))
    P.title(r'%s' %title, fontsize='large')
    P.ylim((0.,1.))
    P.xlabel(r'$v_\mathrm{r} (\mathrm{km/s}) $', fontsize='large')
    P.xticks(fontsize='medium')
    P.yticks(fontsize='medium')
    P.tight_layout()
    cmap = P.colorbar(extend='both', fraction=0.15, shrink=0.7)
    cmap.ax.tick_params(labelsize='medium')
    if save != '':
        P.savefig(save+'.pdf')
    P.show()

def DynSpecI(filename, filename2=None, save=''):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.tri as tri
    matplotlib.rc('xtick', labelsize='medium') 
    matplotlib.rc('ytick', labelsize='medium')
    # Read data
    cycleI, vrI, snI, I, cycleV, vrV, snV, V =  rstokes(filename)
    if filename2 is not None:
        cycleI1, vrI1, snI1, I1, cycleV1, vrV1, snV1, V1 =  rstokes(filename2)
    else:
        cycleI1, vrI1, snI1, I1, cycleV1, vrV1, snV1, V1 =  rstokes(filename[0:-1]+'1')
    # Interpolate
    vmin = min([min(i) for i in vrI])
    vmax = max([max(i) for i in vrI])
    new_vr = N.arange(vmin, vmax, 1.8)
    phaseI = (cycleI%1)
    iPhase = (cycleI%1).argsort()
    tZ = interp_data(vrI, cycleI, I) 
    new_phase, Z = refine(phaseI[iPhase], tZ)
    X, Y = N.meshgrid(new_vr, new_phase)
    phaseI1 = (cycleI1%1)
    iPhase1 = (cycleI1%1).argsort()
    tZ1 = interp_data(vrI1, cycleI1, I1) 
    new_phase1, Z1 = refine(phaseI1[iPhase1], tZ1)
    X1, Y1 = N.meshgrid(new_vr, new_phase1)

    cycleI_uspot, vrI_uspot, snI_uspot, I_uspot, cycleV_uspot, vrV_uspot, snV_uspot, V_uspot =  rstokes(filename[0:-3]+'_unspotted.s1')
    phaseI_uspot = (cycleI_uspot%1)
    iPhase_uspot = (cycleI_uspot%1).argsort()
    tZ_uspot = interp_data(vrI_uspot, cycleI_uspot, I_uspot) 
    new_phase_uspot, Z_uspot = refine(phaseI_uspot[iPhase_uspot], tZ_uspot)
    meanZ = N.ma.array(Z_uspot, mask = Z_uspot == 0).mean(axis=0)

    rZ = N.zeros_like(X)
    rZ1 = N.zeros_like(X1)
    for j in range(new_phase.shape[0]):
        rZ[j, :] = Z[j,:]/meanZ
        rZ1[j, :] = Z1[j,:]/meanZ
    masked_Z = N.ma.array(rZ, mask = Z == 0)
    masked_Z1 = N.ma.array(rZ1, mask = Z1 == 0)

    cmap1 = cmocean.cm.gray_r
    cmap2 = cmocean.cm.tarn
    # Plot dynamic spectra
    fig, ax = P.subplots(1,3,figsize=(8, 4), sharex=True, sharey=True)
    pc=0.003
    im0 = ax[0].pcolormesh(X, Y, masked_Z, vmin=1.-pc, vmax=1.+pc, cmap=cmap1, alpha=1.)
    P.ylim((0,1))
    P.yticks((0,0.2,0.4,0.6,0.8,1.))
    im1 = ax[1].pcolormesh(X1, Y1, masked_Z1, vmin=1.-pc, vmax=1.+pc, cmap=cmap1, alpha=1.)
    im2 = ax[2].pcolormesh(X, Y, 1.0e4*(masked_Z - masked_Z1), vmin=-15, vmax=15, cmap=cmap2, alpha=1.)

    divider1 = make_axes_locatable(ax[1])
    cax1 = divider1.append_axes('right', size='5%', pad=0.05) 
    fig.colorbar(im1, cax=cax1, orientation='vertical', extend='both')

    divider2 = make_axes_locatable(ax[2])
    cax2 = divider2.append_axes('right', size='5%', pad=0.05) 
    fig.colorbar(im2, cax=cax2, orientation='vertical', extend='both')

    divider0 = make_axes_locatable(ax[0])
    cax0 = divider0.append_axes('right', size='0%', pad=0.005) 
    cax0.axis('off')

    ax[0].set_ylabel(r'Rotation cycle', fontsize='large')
    ax[0].set_ylabel(r'Rotation phase', fontsize='large')
    ax[0].set_title(r'Data', fontsize='large')
    ax[1].set_title(r'Model', fontsize='large')
    ax[2].set_title(r'Residuals $\times 10^{4}$', fontsize='large')
    ax[0].set_xlabel(r'$v_\mathrm{r} (\mathrm{km/s}) $', fontsize='large')
    ax[1].set_xlabel(r'$v_\mathrm{r} (\mathrm{km/s}) $', fontsize='large')
    ax[2].set_xlabel(r'$v_\mathrm{r} (\mathrm{km/s}) $', fontsize='large')
    P.tight_layout()
    if save != '':
        P.savefig(save+'.pdf')
    P.show()
