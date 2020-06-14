import numpy as N 
import matplotlib.pyplot as P
import glob 
from math import modf 
from .rspec import rstokes
from .norm import norm, normHa
import os
import cmocean 

# CONSTANTS
lamb_0 = 625.          # central lambda in nm
c = 299792.            # speed of light in km/s
period =  0.5211833875 # from Vaccaro 2015 
To = 2445821.898291


def wmatrix(file, matrix):
    ms = int(matrix.shape[0]/8)
    rms = matrix.shape[0]%8
    for i in range(ms):
        line = tuple(matrix[8*i:8*(i+1)])
        file.write("{:>12.5e} {: >12.5e} {: >12.5e} {: >12.5e} {: >12.5e} {: >12.5e} {: >12.5e} {: >12.5e}\n".format(*line))
    if rms!= 0:
        for j in range(rms):
            file.write("{:>12.5e} ".format(matrix[-rms+j]))
        file.write('\n')
    file.write('\n')

def wfile(wi, phase, vr, snI, I):
    lmin = lamb_0*(vr.min()/c + 1)
    lmax = lamb_0*(vr.max()/c + 1)
    ncont = 1
    dlamb = lamb_0/65000.
    count = 0.
    wi.write('Mean    %1.9e  0   1 0 0 0\n' % lamb_0)
    wi.write('%d  %.6f  %.6f    %1.9e    %1.9e  %1.6e  %1.6e  %1.6e  %.6f\n'% (len(I), phase, phase, lmin, lmax, snI, ncont, dlamb, count))
    wi.write('\n')
    wmatrix(wi, I)

def wfile_V(wi, phase, vr, snV, V):
    lmin = lamb_0*(vr.min()/c + 1)
    lmax = lamb_0*(vr.max()/c + 1)
    ncont = 1
    dlamb = lamb_0/65000.
    count = 0.
    wi.write('Mean    %1.9e  0   0 0 0 1\n' % lamb_0)
    wi.write('%d  %.6f  %.6f    %1.9e    %1.9e  %1.6e  %1.6e  %1.6e  %.6f\n'% (len(V), phase, phase, lmin, lmax, snV, ncont, dlamb, count))
    wi.write('\n')
    wmatrix(wi, V)

def create_night_lsd(year, phase_shift, vr_shift, vr_max, iNorm=False, iSearch=False, caimI=1., caimV=1.):
    if iSearch:
        filename = "v471tau_%s_a%3.3fp%1.6fv%2.3f.ss" %(year, vr_max, abs(phase_shift), vr_shift)
    elif iNorm:
        filename = 'v471tau_'+year+'.ss'
    else:
        filename = 'v471tau_'+year+'.s'
    with open(filename, 'w') as wi:
        files_int = N.genfromtxt('observed_data.txt', dtype=str, comments='*')
        nt_int = int(len(files_int))
        print('Nint observations = ', nt_int)
        files_int.sort()
        t_int = N.zeros(nt_int)
        count = 0
        nCycles = 0
        nexcluded = 0
        for ifile in files_int:
            read_t = open(ifile, 'r')
            t_int[count] = float(read_t.readline())
            data = N.genfromtxt(ifile, skip_header=3)
            cycle = (t_int[count] - To)/period 
            if count == 0:
                nCycles = int(cycle)
            vr = data[:,0] - (vr_max*N.sin(2.*N.pi*(cycle + phase_shift)) + vr_shift)
            I = data[:,1]; sI = data[:,2]
            ics = N.abs(vr) < 110.
            cvr = vr[ics]; cI = I[ics]; csI = sI[ics]
            snI = cI/csI
            cycle -= nCycles #subtract some cycles to make things more readable
            count += 1
            if iNorm:
                normI = norm(cvr, cI, cont=0.9993, fac=5.13) # norm data
            else:
                normI = cI
            # store result
            wfile(wi, cycle, cvr, snI.mean()*N.sqrt(1./caimI), normI)
    with open(filename, 'r') as file:
        # read a list of lines into data
        data = file.readlines()
    # and write everything back
    with open(filename, 'w') as file:
        file.write('LSD profiles of v471tau in 20'+year+'\n')
        file.write('%d\n' %(nt_int))
        file.write('\n')
        file.writelines(data)


def Ha_night_lsd(year, phase_shift, vr_shift, vr_max):
    filename = 'v471tau_'+year+'.ss'
    with open(filename, 'w') as wi:
        files_int = N.genfromtxt('observed_data.txt', dtype=str, comments='*')
        nt_int = int(len(files_int))
        print('Nint observations = ', nt_int)
        files_int.sort()
        t_int = N.zeros(nt_int)
        count = 0
        nCycles = 0
        nexcluded = 0
        for ifile in files_int:
            read_t = open(ifile, 'r')
            t_int[count] = float(read_t.readline())
            data = N.genfromtxt(ifile, skip_header=3)
            cycle = (t_int[count] - To)/period 
            if count == 0:
                nCycles = int(cycle)
            vr = data[:,0] - (vr_max*N.sin(2.*N.pi*(cycle + phase_shift)) + vr_shift)
            I = data[:,1]; sI = data[:,2]
            ics = N.abs(vr) < 250.
            cvr = vr[ics]; cI = I[ics]; csI = sI[ics]
            snI = cI/csI
            cycle -= nCycles #subtract some cycles to make things more readable
            count += 1
            normI = normHa(cvr, cI) 
            P.plot(cvr, cI, 'k',linewidth=1) 
            P.text(-100., .95,s=r'$\phi$ = %1.4f, $v_\mathrm{rad}$ = %2.1f, k = %3.1f' %(phase_shift, vr_shift, vr_max), fontsize=10)
            P.plot(N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
            P.plot(-N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
            # store result
            wfile(wi, cycle, cvr, snI.mean()*N.sqrt(1.), normI)
    with open(filename, 'r') as file:
        # read a list of lines into data
        data = file.readlines()
    # and write everything back
    with open(filename, 'w') as file:
        file.write('Ha LSD profiles of v471tau in 20'+year+'\n')
        file.write('%d\n' %(nt_int))
        file.write('\n')
        file.writelines(data)

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

def SpecI_night(filename, iSort=False, title='',save=''):
    cycle_tmp, vr_tmp, snI_tmp, I_tmp, cycleV_tmp, vrV_tmp, snV_tmp, V_tmp =  rstokes('../%s' %filename)
    Z_tmp = interp_data(vr_tmp, cycle_tmp, I_tmp)
    phase_tmp = (cycle_tmp%1)
    iPhase_tmp = (cycle_tmp%1).argsort()
    new_phase_tmp, Z_tmp2 = refine(phase_tmp[iPhase_tmp], Z_tmp)
    masked_Z_tmp = N.ma.array(Z_tmp2, mask = Z_tmp2 == 0) 
    meanZ = masked_Z_tmp.mean(axis=0)

    # Read data
    cycleI, vrI, snI, I, cycleV, vrV, snV, V =  rstokes(filename)
    vmin = min([min(i) for i in vrI])
    vmax = max([max(i) for i in vrI])
    new_vr = N.arange(vmin, vmax, 1.8)

    # Plot dynamic spectra
    cmap = cmocean.cm.gray_r  # Colormap used in plot
    P.figure(figsize=(4,5))
    ax1 = P.subplot(111)
    pc = 0.003
    tZ = interp_data(vrI, cycleI, I)
    phaseI = (cycleI%1)
    iPhase = (cycleI%1).argsort()
    new_phase, Z = refine(phaseI[iPhase], tZ)
    X, Y = N.meshgrid(new_vr, new_phase)
    rZ = N.zeros_like(X)
    for j in range(new_phase.shape[0]):
        rZ[j, :] = Z[j,:]/meanZ
    masked_Z = N.ma.array(rZ, mask = Z == 0) 
    P.pcolormesh(X, Y, masked_Z, vmin=1.-pc, vmax=1.+pc, cmap=cmap)
    P.plot( 78*N.ones(15), N.linspace(0.,1.,15), 'r--')
    P.plot(-78*N.ones(15), N.linspace(0.,1.,15), 'r--')
    P.ylabel(r'Rotation phase', fontsize='large')
    P.title(r'%s' %title, fontsize='large')
    P.ylim((0,1))
    P.xlabel(r'$v_\mathrm{r} (\mathrm{km/s}) $', fontsize='large')
    P.xticks(fontsize='medium')
    P.yticks(fontsize='medium')
    P.tight_layout()
    cmap = P.colorbar(extend='both', fraction=0.15, shrink=0.7)
    cmap.ax.tick_params(labelsize='medium')
    if save != '':
        P.savefig(save+'.pdf')
    P.show()
