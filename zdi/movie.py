import numpy as N 
import matplotlib.pyplot as P
import glob 
import os
import matplotlib.animation as animation
from pylab import *

P.style.use('ggplot')
P.rcParams['font.family'] = 'arial'

# CONSTANTS
lamb_0 = 625.          # central lambda in nm
c = 299792.            # speed of light in km/s
period =  0.5211833875 # from Vaccaro 2015 
To = 2445821.898291

def build_matrix(year, phase_shift):
    PATH = '/Users/bzaire/lsd/spectra/v471tau/v471tau_int' + year
    os.chdir(PATH)
    list_obs = N.genfromtxt('observed_data.txt', dtype=str)
    remove_obs = N.genfromtxt('observed_data_excluded.txt', dtype=str)
    files_int = list(filter(lambda i: i not in remove_obs, list_obs))  # Only take data in a give level of signal to noise
    nt_int = int(len(files_int))
    files_int.sort()
    cycle = []
    I = []
    vr = []
    count = 0
    # Read data at each file
    for ifile in files_int:
        read_t = open(ifile, 'r')
        t_int = float(read_t.readline())
        data = N.genfromtxt(ifile, skip_header=3)
        cycle.append((t_int - To)/period + phase_shift) 
        if count == 0:
            nCycles = int(cycle[count]) # Subtracts a given number of cycles
        I.append(data[:,1])
        vr.append(data[:,0])
        count += 1
    cycle = N.array(cycle) # make it a numpy vector
    cycle -= nCycles #subtract some cycles to make things more readable
    return cycle, vr, I

def movie(year, phase_shift):
    cycle, vr, I = build_matrix(year, phase_shift)
    # Args sort in phase
    isp = (cycle%1).argsort()
    # Start to produce the movie
    dpi = 100
    fig = P.figure()
    ax = P.subplot(1,1,1)
    line, = ax.plot([],[], 'k')
    P.xlim(-350, 350)
    P.ylim(0.95,1.01) 
    P.xlabel(r'$v_\mathrm{r} (km/s)$', fontsize = 25)
    P.ylabel(r'$\mathrm{stokes}-I$', fontsize = 25, rotation = 90)
    label_text = ax.text(0, 0.952, '', color='red', fontsize=20, horizontalalignment='center', verticalalignment='center')
    y = N.linspace(0.98, 1.005, 15)
    vr_min = -149.11 + 34.57
    vr_max = 149.11 + 34.57 
    P.plot(N.ones(15)*vr_max, y, '--k', alpha=1./3)
    P.plot(N.ones(15)*vr_min, y, '--k', alpha=1./3)
    ax.text(vr_min, 1.007, r'$\phi=0.75$', color='red', fontsize=18, horizontalalignment='center', verticalalignment='center')
    ax.text(vr_max, 1.007, r'$\phi=0.25$', color='red', fontsize=18, horizontalalignment='center', verticalalignment='center')


    fig.set_size_inches([5,5])
    P.tight_layout()

    def update_img(n):
        line.set_data(vr[isp[n]], I[isp[n]])
        #P.legend([r'Cycle = %1.4f' %(cycle[isp[n]]%1)], loc='lower right', fontsize=14, fancybox=None, shadow=None, edgecolor=None) 
        label_text.set_text(r'Phase = %1.4f' %(cycle[isp[n]]%1))
        return [line,] + [label_text,]

    ani = animation.FuncAnimation(fig, update_img, frames=cycle.shape[0], repeat=True)
    writer = animation.writers['ffmpeg'](fps=20)

    ani.save('stokes_I_%s.mp4' %year, writer=writer,dpi=dpi)
    return ani

