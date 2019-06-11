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

def read_matriz(ofile, nlines):
    # Local function.
    # Read matrix
    info = [next(ofile) for j in range(4)]
    phase = float(info[2].split()[1])
    sn = float(info[2].split()[5])
    tmp = [next(ofile).split() for i in range(nlines)]
    matriz = [float(tmp[i][j]) for i in range(nlines) for j in range(len(tmp[i]))]
    return matriz, phase, sn

def rstokes(filedata):
    # Global function
    """ Read the phase and the stokes parameters available (I and/or V).
        Standard format produced by .ss and .s1 files.
        Parameters:
            filedata: str
                Name of the file with the stokes data.
        Return
            out: ndarray
                Phase of measurement, radial velocity and stokes parameters.
    """
    with open(filedata, 'r') as ofile:
        # Reading basic informations to format the data
        head = [next(ofile) for j in range(200)]
        nvr = int(head[4].split()[0])
        lmin, lmax = float(head[4].split()[3]), float(head[4].split()[4])
        lc = float(head[3].split()[1])
        c = 299792. #  Speed of light in km/s
        vrmin, vrmax = (lmin-lc)*c/lc, (lmax-lc)*c/lc
        vr = N.linspace(vrmin, vrmax, nvr, endpoint=True)
        nlines = get_ml(nvr)
        nstokes = int(1 + int(head[7+nlines].split()[-1]))
        nphase = int(int(head[1].split()[0])/nstokes)
        del(head, lmin, lmax, lc, c)

    with open(filedata, 'r') as ofile:
        # Reading stokes parameters
        head = [next(ofile) for j in range(2)]
        del(head)
        phase = N.zeros(nphase)
        snI = N.zeros(nphase)
        stokesI = N.zeros((nphase, nvr))
        if nstokes == 2:
            stokesV = N.zeros((nphase, nvr))
            snV = N.zeros(nphase)
            for ip in range(nphase):
                stokesI[ip,:], phase[ip], snI[ip] = read_matriz(ofile, nlines)
                stokesV[ip,:], phase[ip], snV[ip] = read_matriz(ofile, nlines)
            return phase, vr, snI, stokesI, snV, stokesV
        elif nstokes ==1:
            for ip in range(nphase):
                stokesI[ip,:], phase[ip], snI[ip] = read_matriz(ofile, nlines)
            return phase, vr, snI, stokesI
        del(nphase, nvr, nstokes, vrmin, vrmax)
