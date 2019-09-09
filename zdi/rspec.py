import numpy as N
# Routines to read the data
# __init__ file should contain the following line:
# from rspec import rstokes

def get_ml(npoints):
    # Local function.
    # Defines the number of lines in each matrix.
    # # WARNING: This function assumes that the standard
    #        output with 8 numbers per line is preserved.
    if (npoints % 8) != 0:
        return int(npoints/8 + 1)
    else:
        return int(npoints/8)

def read_block(ofile):
    # Local function.
    # Read block with data. 
    info = [next(ofile) for j in range(4)]
    which_stokes = int(info[1].split()[-1])  # Cycle of observation
    cycle = float(info[2].split()[1])  # Cycle of observation
    nvr = int(info[2].split()[0])  # Number of points in radial velocity
    nlines = get_ml(nvr)  # Find the numbers of lines which contain the data
    lc = float(info[1].split()[1])   # Central wavelength
    lmin, lmax = float(info[2].split()[3]), float(info[2].split()[4])  # Min, Max wavelengths
    c = 299792. #  Speed of light in km/s
    vrmin, vrmax = (lmin-lc)*c/lc, (lmax-lc)*c/lc  # Min, Max radial velocities
    vr = N.linspace(vrmin, vrmax, nvr, endpoint=True)
    sn = float(info[2].split()[5])    # Signal to Noise
    tmp = [next(ofile).split() for i in range(nlines)]
    matriz = [float(tmp[i][j]) for i in range(nlines) for j in range(len(tmp[i]))]
    return vr, matriz, cycle, sn, which_stokes

def rstokes(filedata):
    # Global function
    """ Read the cycle and the stokes parameters available (I and/or V).
        Standard format produced by .ss and .s1 files.
        Parameters:
            filedata: str
                Name of the file with the stokes data.
        Return
            out: ndarray
                Phase of measurement, radial velocity and stokes parameters.
    """
    with open(filedata, 'r') as ofile:
        # Reading stokes parameters
        head = [next(ofile) for j in range(2)]
        nObs = int(head[1].split()[0])
        del(head)
        # Set list with variables
        cycleI = []; vrI = []
        snI = []; I = []
        cycleV = []; vrV = []
        snV = []; V = []
        for iObs in range(nObs):
            tempVr, tempStokes, tempCycle, tempSn, which_stokes = read_block(ofile)
            if which_stokes == 0:
                vrI.append(tempVr)
                cycleI.append(tempCycle)
                snI.append(tempSn)
                I.append(tempStokes)
            else:
                vrV.append(tempVr)
                cycleV.append(tempCycle)
                snV.append(tempSn)
                V.append(tempStokes)
        return cycleI, vrI, snI, I, cycleV, vrV, snV, V
