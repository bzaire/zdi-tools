import numpy as N 
import matplotlib.pyplot as P
import glob 
from math import modf 
from .rspec import rstokes
from .norm import norm, norm_search
import os

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

def create_lsd(year, phase_shift, vr_shift, vr_max, iNorm=False):
    if iNorm:
        filename = 'v471tau_'+year+'.ss'
    else:
        filename = 'v471tau_'+year+'.s'
    with open(filename, 'w') as wi:
        PATH = '/Users/bzaire/lsd/spectra/v471tau/v471tau_int' + year
        os.chdir(PATH)
        files_int = N.genfromtxt('observed_data.txt', dtype=str)
        nt_int = int(len(files_int))
        print('Nint observations = ', nt_int)
        files_int.sort()
        t_int = N.zeros(nt_int)
        P.figure(figsize=(8,5))
        count = 0
        nexcluded = 0
        with open('observed_data_excluded.txt', 'w') as wud:
            for ifile in files_int:
                read_t = open(ifile, 'r')
                t_int[count] = float(read_t.readline())
                data = N.genfromtxt(ifile, skip_header=3)
                cycle = (t_int[count] - To)/period + phase_shift 
                vr = data[:,0] - (vr_max*N.sin(2.*N.pi*cycle) + vr_shift)
                I = data[:,1]; sI = data[:,2]
                snI = I/sI
                ics = N.abs(vr) < 110.
                cvr = vr[ics]; cI = I[ics]
                cycle -= 14420. #subtract some cycles to make things more readable
                count += 1
                if(snI[0:25].mean()>=400.):
                    P.plot(cvr, cI, 'k',linewidth=1) 
                    P.text(-100., .95,s=r'$\phi$ = %1.4f, $v_\mathrm{rad}$ = %2.1f, k = %3.1f' %(phase_shift, vr_shift, vr_max), fontsize=10)
                    P.plot(N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
                    P.plot(-N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
                    if iNorm:
                        normI = norm(cvr, cI, cont=0.9993, fac=5.13) # norm data
                    else:
                        normI = cI
                    # store result
                    wfile(wi, cycle, cvr, snI.mean()*N.sqrt(1.), normI)
                else:
                    wud.write(ifile+'\n')
                    nexcluded += 1
        PATH = '/Users/bzaire/lsd/spectra/v471tau/v471tau_pol' + year
        os.chdir(PATH)
        files_pol = N.genfromtxt('observed_data.txt', dtype=str)
        nt_pol = int(len(files_pol))
        print('Npol observations = ', nt_pol)
        files_pol.sort()
        t_pol = N.zeros(nt_pol)
        count = 0
        with open('observed_data_excluded.txt', 'w') as wud:
            for ifile in files_pol:
                read_t = open(ifile, 'r')
                t_pol[count] = float(read_t.readline())
                data = N.genfromtxt(ifile, skip_header=3)
                cycle = (t_pol[count] - To)/period + phase_shift 
                vr = data[:,0] - (vr_max*N.sin(2.*N.pi*cycle) + vr_shift)
                V = data[:,3]; sV = data[:,4]
                snV = 1./sV
                Null = data[:,5]; sN = data[:,6]
                snN = 1./sN
                cV = V[ics]; cNull = Null[ics]
                cycle -= 14420. #subtract some cycles to make things more readable
                if(snV.mean()>=1500.):
                    # store result
                    wfile_V(wi, cycle, cvr, snV.mean(), cV)
                else:
                    nexcluded += 1
    PATH = '/Users/bzaire/lsd/spectra/v471tau/'
    os.chdir(PATH)
    print('Excluding %d observations in 20%2s' %(nexcluded, year))
    with open(filename, 'r') as file:
        # read a list of lines into data
        data = file.readlines()
    # and write everything back
    with open(filename, 'w') as file:
        file.write('LSD profiles of v471tau in 20'+year+'\n')
        file.write('%d\n' %(nt_int+nt_pol - nexcluded))
        file.write('\n')
        file.writelines(data)
    print(nt_int+nt_pol-nexcluded)
