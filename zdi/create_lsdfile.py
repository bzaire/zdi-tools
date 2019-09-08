import numpy as N 
import matplotlib.pyplot as P
import glob 
from math import modf 
from .rspec import rstokes
from .norm import norm, norm_search

# constants
lamb_0 = 625. # central lambda in nm
c = 299792. # speed of light in km/s
period =  0.521183398 #days
To = 2440610.06406

def wmatrix(file, matrix):
    ms = int(matrix.shape[0]/8)
    rms = matrix.shape[0]%8
    for i in range(ms):
        line = tuple(matrix[8*i:8*(i+1)])
        file.write("{:>12.5e} {: >12.5e} {: >12.5e} {: >12.5e} {: >12.5e} {: >12.5e} {: >12.5e} {: >12.5e}\n".format(*line))
#        file.write('   %1.5e     %1.5e     %1.5e     %1.5e     %1.5e     %1.5e     %1.5e     %1.5e\n' % line)
    if rms!= 0:
        for j in range(rms):
            file.write("{:>12.5e} ".format(matrix[-rms+j]))
        file.write('\n')
    file.write('\n')

#def wmatrix(file, temp_matrix):
#    matrix = temp_matrix
#    ms = int(matrix.shape[0]/8)
#    rms = matrix.shape[0]%8
#    for i in range(ms):
#        line = tuple(matrix[8*i:8*(i+1)])
#        file.write('   %1.5e     %1.5e     %1.5e     %1.5e     %1.5e     %1.5e     %1.5e     %1.5e\n' % line)
#    if rms!= 0:
#        for j in range(rms):
#            file.write('   %.5e  ' % matrix[-rms+j])
#        file.write('\n')
#    file.write('\n')

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

def wfile_V(wi, phase, vr, snI, I, snV, V):
    lmin = lamb_0*(vr.min()/c + 1)
    lmax = lamb_0*(vr.max()/c + 1)
    ncont = 1
    dlamb = lamb_0/65000.
    count = 0.
    wi.write('Mean    %1.9e  0   1 0 0 0\n' % lamb_0)
    wi.write('%d  %.6f  %.6f    %1.9e    %1.9e  %1.6e  %1.6e  %1.6e  %.6f\n'% (len(I), phase, phase, lmin, lmax, snI, ncont, dlamb, count))
    wi.write('\n')
    wmatrix(wi, I)
    wi.write('Mean    %1.9e  0   0 0 0 1\n' % lamb_0)
    wi.write('%d  %.6f  %.6f    %1.9e    %1.9e  %1.6e  %1.6e  %1.6e  %.6f\n'% (len(V), phase, phase, lmin, lmax, snV, ncont, dlamb, count))
    wi.write('\n')
    wmatrix(wi, V)

def create_lsd_05(filename, phase_shift, vr_shift, vr_max, stokesV=False):
    if stokesV:
        files = glob.glob('lsd05_*')
        nt = 2*len(files)
    else:
        files = glob.glob('I_lsd05_*')
        nt = len(files)
    files.sort()
    print(nt)
    t = N.zeros(nt)
    count = 0
    for ift in files:
        if stokesV:
            with open('info'+ift[3:], 'r') as rt:
                t[count] = float(rt.readline())
                count += 1
        else:
            with open('I_info'+ift[5:], 'r') as rt:
                t[count] = float(rt.readline())
                count += 1
    P.figure(figsize=(8,5))
    count = 0
    nexcluded = 0
    with open(filename, 'w') as wi:
        wi.write('LSD profiles of v471tau in 2005\n')
        wi.write('%d\n' % (nt))
        wi.write('\n')
        for ifile in files:
            data = N.genfromtxt(ifile, skip_header=2)
            vr = data[:,0] - (vr_max*N.sin(2.*N.pi*((t[count] - To)/period + phase_shift)) + vr_shift)
            I = data[:,1]; sI = data[:,2]
            snI = I/sI
            ics = N.abs(vr) < 110.
            cvr = vr[ics]; cI = I[ics]
            cycle = (t[count] - To)/period + phase_shift - 25151
            count += 1
            if stokesV:
                V = data[:,3]; sV = data[:,4]
                snV = 1./sV
                Null = data[:,5]; sN = data[:,6]
                snN = 1./sN
                cV = V[ics]; cNull = Null[ics]
                if((snI[0:25].mean()>=500.) and (snV.mean()>=1000.)):
                    P.plot(cvr, cI, 'k',linewidth=1) 
                    P.text(-100., .95,s=r'$\phi$ = %1.4f, $v_\mathrm{rad}$ = %2.1f, k = %3.1f' %(phase_shift, vr_shift, vr_max), fontsize=10)
                    P.plot(N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
                    P.plot(-N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
                    # store result
                    #wfile_V(wi, cycle, cvr, snI.mean()*N.sqrt(1./0.35), cI, snV.mean(), cV)
                    wfile_V(wi, cycle, cvr, snI.mean()*N.sqrt(1.), cI, snV.mean(), cV)
                else:
                    nexcluded += 2
            else: 
                if(snI[0:25].mean()>=300.):
                    P.plot(cvr, cI, 'k',linewidth=1) 
                    P.text(-100., .95,s=r'$\phi$ = %1.4f, $v_\mathrm{rad}$ = %2.1f, k = %3.1f' %(phase_shift, vr_shift, vr_max), fontsize=10)
                    P.plot(N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
                    P.plot(-N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
                    # store result
                    wfile(wi, cycle, cvr, snI.mean()*N.sqrt(1.), cI)
                    #wfile(wi, cycle, cvr, snI.mean()*N.sqrt(1./0.35), cI)
                else:
                    nexcluded += 1
    if nexcluded != 0:
        print('Excluding %d nights in 2005' %nexcluded)
        with open(filename, 'r') as file:
            # read a list of lines into data
            data = file.readlines()
        # Changing second line
        data[1] = '%d\n' %(nt - nexcluded)
        # and write everything back
        with open(filename, 'w') as file:
            file.writelines(data)
    print(nt-nexcluded)

def create_lsd_04(filename, phase_shift, vr_shift, vr_max, stokesV=False):
    if stokesV:
        files = glob.glob('lsd04_*')
        print(len(files))
        nt = int(2.*len(files))
        print(nt)
    else:
        files = glob.glob('I_lsd04_*')
        nt = len(files)
    files.sort()
    print(nt)
    t = N.zeros(nt)
    count = 0
    for ift in files:
        if stokesV:
            with open('info'+ift[3:], 'r') as rt:
                t[count] = float(rt.readline())
                count += 1
        else:
            with open('I_info'+ift[5:], 'r') as rt:
                t[count] = float(rt.readline())
                count += 1
    P.figure(figsize=(8,5))
    count = 0
    nexcluded = 0
    with open(filename, 'w') as wi:
        wi.write('LSD profiles of v471tau in 2004\n')
        wi.write('%d\n' % (nt))
        wi.write('\n')
        for ifile in files:
            data = N.genfromtxt(ifile, skip_header=2)
            vr = data[:,0] - (vr_max*N.sin(2.*N.pi*((t[count] - To)/period + phase_shift)) + vr_shift)
            I = data[:,1]; sI = data[:,2]
            snI = I/sI
            ics = N.abs(vr) < 110.
            cvr = vr[ics]; cI = I[ics]
            cycle = (t[count] - To)/period + phase_shift - 24420
            count += 1
            if stokesV:
                V = data[:,3]; sV = data[:,4]
                snV = 1./sV
                Null = data[:,5]; sN = data[:,6]
                snN = 1./sN
                cV = V[ics]; cNull = Null[ics]
                if((snI[0:25].mean()>=500.) and (snV.mean()>=1000.)):
                    P.plot(cvr, cI, 'k',linewidth=1) 
                    P.text(-100., .95,s=r'$\phi$ = %1.4f, $v_\mathrm{rad}$ = %2.1f, k = %3.1f' %(phase_shift, vr_shift, vr_max), fontsize=10)
                    P.plot(N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
                    P.plot(-N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
                    # store result
#                    wfile_V(wi, cycle, cvr, snI.mean()*N.sqrt(1./0.35), cI, snV.mean(), cV)
                    wfile_V(wi, cycle, cvr, snI.mean()*N.sqrt(1.), cI, snV.mean(), cV)
                else:
                    nexcluded += 2
            else: 
                if(snI[0:25].mean()>=0.009):
                    P.plot(cvr, cI, 'k',linewidth=1) 
                    P.text(-100., .95,s=r'$\phi$ = %1.4f, $v_\mathrm{rad}$ = %2.1f, k = %3.1f' %(phase_shift, vr_shift, vr_max), fontsize=10)
                    P.plot(N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
                    P.plot(-N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
                    # store result
                    #wfile(wi, cycle, cvr, snI.mean()*N.sqrt(1./0.35), cI)
                    wfile(wi, cycle, cvr, snI.mean()*N.sqrt(1.), cI)
                else:
                    nexcluded += 1
    if nexcluded != 0:
        print('Excluding %d nights in 2004' %nexcluded)
        with open(filename, 'r') as file:
            # read a list of lines into data
            data = file.readlines()
        # Changing second line
        data[1] = '%d\n' %(nt - nexcluded)
        # and write everything back
        with open(filename, 'w') as file:
            file.writelines(data)
    print(nt-nexcluded)

    
def norm_lsd_05(filename, fac=5.12, stokesV=False):
    if stokesV:
        phase, vr, snI, I, snV, V = rstokes(filename)
        nObs = len(phase)
        rnI = norm(vr, I, cont=0.9993, fac=fac)
        with open(filename+"s", 'w') as wi:
            wi.write('Spectrum of v471tau in 2005\n')
            wi.write('%d\n' % (2*nObs))
            wi.write('\n')
            for i in range(nObs):
                wfile_V(wi, phase[i], N.array(vr[i]), snI[i], N.array(rnI[i]), snV[i], N.array(V[i]))
    else:
        phase, vr, snI, I = rstokes(filename)
        nObs = len(I)
        rnI = norm(vr, I, cont=0.9993, fac=fac)
        with open(filename+"s", 'w') as wi:
            wi.write('Spectrum of v471tau in 2005\n')
            wi.write('%d\n' % (nObs))
            wi.write('\n')
            for i in range(nObs):
                wfile(wi, phase[i], N.array(vr[i]), snI[i], N.array(rnI[i]))

def norm_lsd_04(filename, fac=5.12, stokesV=False):
    if stokesV:
        phase, vr, snI, I, snV, V = rstokes(filename)
        nObs = len(phase)
        rnI = norm(vr, I, cont=0.9993, fac=fac)
        with open(filename+"s", 'w') as wi:
            wi.write('Spectrum of v471tau in 2005\n')
            wi.write('%d\n' % (2*nObs))
            wi.write('\n')
            for i in range(nObs):
                wfile_V(wi, phase[i], N.array(vr[i]), snI[i], N.array(rnI[i]), snV[i], N.array(V[i]))
    else:
        phase, vr, snI, I = rstokes(filename)
        nObs = len(phase)
        rnI = norm(vr, I, cont=0.9993, fac=fac)
        with open(filename+"s", 'w') as wi:
            wi.write('Spectrum of v471tau in 2004\n')
            wi.write('%d\n' % (nObs))
            wi.write('\n')
            for i in range(nObs):
                wfile(wi, phase[i], N.array(vr[i]), snI[i], N.array(rnI[i]))


def norm_lsd_search(filename, fit_cont, area):
    phase, vr, snI, I = rstokes(filename)
    nObs = len(I)
    rnI = norm_search(vr, I, fit_cont, area, cont=0.9993, fac=5.12)
    with open(filename+"s", 'w') as wi:
        wi.write('Spectrum of v471tau in 2004\n')
        wi.write('%d\n' % (nObs))
        wi.write('\n')
        for i in range(nObs):
            wfile(wi, phase[i], N.array(vr[i]), snI[i], N.array(rnI[i]))
