import numpy as N
import matplotlib.pyplot as P

def norm(vr, matrix, cont=0.9992, fac = 5.1):
    norm_matrix = N.zeros_like(matrix)
    nObs = matrix.shape[0]
    area = N.zeros(nObs)
    areaf = N.zeros(nObs)
    rescale = N.zeros_like(matrix)
    vr_max = N.max(N.abs(vr))
    for i in range(nObs):
        stokes = matrix[i,:]
        ir = [0,1,2,3,4,5]; il = [-1,-2,-3,-4,-5,-6]
        vrr = vr[ir].mean(); vrl = vr[il].mean()
        stokesr = stokes[ir].mean(); stokesl = stokes[il].mean()
        fit_cont = N.poly1d(N.polyfit((vrr, vrl), (stokesr, stokesl), deg=1))
        norm_matrix[i,:] = stokes/fit_cont(vr)
        area[i] = N.trapz(1-norm_matrix[i,:], dx=1.8)
        rescale[i,:] = cont - (1 - norm_matrix[i,:])*(fac/area[i])
#        areaf[i] = N.trapz(cont - rescale[i,:], dx=1.8)
    return rescale

def norm_search(vr, matrix, fit_cont, area, cont=0.9992, fac = 5.1):
    # input:
    #        fit_cont : list with function containing the 1th degree polynomial for each observation.
    #        area : list with the areas computed for each observation.
    norm_matrix = N.zeros_like(matrix)
    nObs = matrix.shape[0]
    rescale = N.zeros_like(matrix)
    for i in range(nObs):
        stokes = matrix[i,:]
        norm_matrix[i,:] = stokes/fit_cont[i](vr)
        rescale[i,:] = cont - (1 - norm_matrix[i,:])*(fac/area[i])
    return rescale

def find_cont(vr, matrix, cont=0.9992, fac = 5.1):
    norm_matrix = N.zeros_like(matrix)
    nObs = matrix.shape[0]
    area = N.zeros(nObs)
    fit_cont = N.zeros(nObs)
    rescale = N.zeros_like(matrix)
    vr_max = N.max(N.abs(vr))
    for i in range(nObs):
        stokes = matrix[i,:]
        ir = [0,1,2,3,4,5]; il = [-1,-2,-3,-4,-5,-6]
        vrr = vr[ir].mean(); vrl = vr[il].mean()
        stokesr = stokes[ir].mean(); stokesl = stokes[il].mean()
        fit_cont[i] = N.poly1d(N.polyfit((vrr, vrl), (stokesr, stokesl), deg=1))
        norm_matrix[i,:] = stokes/fit_cont(vr)
        area[i] = N.trapz(1-norm_matrix[i,:], dx=1.8)
        rescale[i,:] = cont - (1 - norm_matrix[i,:])*(fac/area[i])
    return fit_cont, area, rescale
