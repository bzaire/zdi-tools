import numpy as N 
import matplotlib.pyplot as P
import glob 
from math import modf 
from .rspec import rstokes
from .norm import norm

# constants
lamb_0 = 625. # central lambda in nm
c = 299792. # speed of light in km/s
period =  0.521183398 #days
To = 2440610.06406

def wmatrix(file, temp_matrix):
    matrix = temp_matrix
    ms = int(matrix.shape[0]/8)
    rms = matrix.shape[0]%8
    for i in range(ms):
        line = tuple(matrix[8*i:8*(i+1)])
        file.write('   %1.5e     %1.5e     %1.5e     %1.5e     %1.5e     %1.5e     %1.5e     %1.5e\n' % line)
    if rms!= 0:
        for j in range(rms):
            file.write('   %.5e  ' % matrix[-rms+j])
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

def create_lsd(filename, phase_shift, vr_shift, vr_max):
    files = glob.glob('I_lsd05_*')
    files.sort()
    nt = len(files)
    print(nt)
    t = N.zeros(nt)
    count = 0
    for ift in files:
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
            #phase = N.mod((t[count] - To), period)/period# + phase_shift 
            cycle = (t[count] - To)/period + phase_shift - 25151
            count += 1
            if (snI[0:25].mean()>=500.):
                P.plot(cvr, cI, 'k',linewidth=1) 
                P.text(-100., .95,s=r'$\phi$ = %1.4f, $v_\mathrm{rad}$ = %2.1f, k = %3.1f' %(phase_shift, vr_shift, vr_max), fontsize=10)
                P.plot(N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
                P.plot(-N.ones(10)*110.,N.linspace(0.98*cI.min(),1.02*cI.max(),10), '--k')
                # store result
                wfile(wi, cycle, cvr, snI.mean()*N.sqrt(1./0.4), cI)
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

    
def norm_lsd(filename):
    phase, vr, snI, I = rstokes(filename)
    nObs = len(I)
    rnI = norm(vr, I, cont=0.9993, fac=5.12)
    with open(filename+"s", 'w') as wi:
        wi.write('Spectrum of v471tau in 2005\n')
        wi.write('%d\n' % (nObs))
        for i in range(nObs):
            wfile(wi, phase[i], N.array(vr[i]), snI[i], N.array(rnI[i]))
