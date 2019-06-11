import numpy as N
import matplotlib.pyplot as P

def norm2(vr, matrix):
    norm_matrix = N.zeros_like(matrix)
    nObs = matrix.shape[0]
    area = N.zeros(nObs)
    rescale = N.zeros_like(matrix)
    jpoints = [0,1,2,3,-4,-3,-2,-1]
    lev=0.99
    for i in range(nObs):
        stokes = matrix[i,:]
        fit_cont = N.poly1d(N.polyfit(vr[jpoints], stokes[jpoints], deg=1))
        norm_matrix[i,:] = stokes/fit_cont(vr)
        area[i] = N.trapz(1/norm_matrix[i,:]-1, dx=1.8)
    fac = area.mean()/area
    for i in range(nObs):
        rescale[i,:] = norm_matrix[i,:]*fac[i]
        shift_cont = N.polyfit(vr[jpoints], rescale[jpoints], deg=0)[0]
    return rescale - shift_cont + 1.

def norm1(vr, matrix):
    norm_matrix = N.zeros_like(matrix)
    nObs = matrix.shape[0]
    area = N.zeros(nObs)
    rescale = N.zeros_like(matrix)
    for i in range(nObs):
        stokes = matrix[i,:]
        rs = stokes[vr>0]
        ls = stokes[vr<0]
        rmax = rs.max()
        lmax = ls.max()
        lev=0.99
        ir, = N.where( stokes == rs[rs>lev*rmax][0] )
        il, = N.where( stokes == ls[ls>lev*lmax][-1] )
        print(ir,il)
        ipoints = [il[-1],ir[-1]]
        fit_cont = N.poly1d(N.polyfit(vr[ipoints], stokes[ipoints], deg=1))
        norm_matrix[i,:] = lev*stokes/fit_cont(vr)
        area[i] = N.trapz(1/norm_matrix[i,:]-1, dx=1.8)
    fac = area.mean()/area
    for i in range(nObs):
        rescale[i,:] = norm_matrix[i,:]*fac[i]
        jpoints = [0,1,2,3,-4,-3,-2,-1]
        shift_cont = N.polyfit(vr[jpoints], rescale[jpoints], deg=0)[0]
    return rescale - shift_cont + 1.

def norm(vr, matrix, cont=0.9992, fac = 5.1):
    norm_matrix = N.zeros_like(matrix)
    nObs = matrix.shape[0]
    area = N.zeros(nObs)
    areaf = N.zeros(nObs)
    rescale = N.zeros_like(matrix)
    vr_max = N.max(N.abs(vr))
    for i in range(nObs):
        stokes = matrix[i,:]
        ir = [0,1,2]; il = [-1,-2,-3]
        vrr = vr[ir].mean(); vrl = vr[il].mean()
        stokesr = stokes[ir].mean(); stokesl = stokes[il].mean()
        fit_cont = N.poly1d(N.polyfit((vrr, vrl), (stokesr, stokesl), deg=1))
        norm_matrix[i,:] = stokes/fit_cont(vr)
        area[i] = N.trapz(1-norm_matrix[i,:], dx=1.8)
        rescale[i,:] = cont - (1 - norm_matrix[i,:])*(fac/area[i])
#        areaf[i] = N.trapz(cont - rescale[i,:], dx=1.8)
    return rescale
