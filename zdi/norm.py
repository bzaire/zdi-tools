import numpy as N
import matplotlib.pyplot as P

def norm(vr_matrix, matrix, cont=0.9992, fac = 5.1):
    nObs = len(matrix)
    rescale = []
    for i in range(nObs):
        stokes = N.array(matrix[i])
        vr = N.array(vr_matrix[i])
        ir = [0,1,2]; il = [-1,-2,-3]
        vrr = vr[ir].mean(); vrl = vr[il].mean()
        stokesr = stokes[ir].mean(); stokesl = stokes[il].mean()
        fit_cont = N.poly1d(N.polyfit((vrr, vrl), (stokesr, stokesl), deg=1))
        area = N.trapz(1- stokes/fit_cont(vr), dx=1.8)
        rescale.append(cont - (1 - stokes/fit_cont(vr))*(fac/area))
    return rescale

def norm_search(vr_matrix, matrix, fit_cont, area, cont=0.9992, fac = 5.1):
    # input:
    #        fit_cont : list with function containing the 1th degree polynomial for each observation.
    #        area : list with the areas computed for each observation.
    nObs = len(matrix)
    rescale = []
    for i in range(nObs):
        stokes = N.array(matrix[i])
        vr = N.array(vr_matrix[i])
        ir = [0,1,2]; il = [-1,-2,-3]
        rescale.append(cont - (1 - stokes/fit_cont[i](vr))*(fac/area[i]))
    return rescale

def find_cont(vr_matrix, matrix, cont=0.9992, fac = 5.1):
    nObs = len(matrix)
    area = N.zeros(nObs)
    rescale = []
    for i in range(nObs):
        stokes = N.array(matrix[i])
        vr = N.array(vr_matrix[i])
        ir = [0,1,2]; il = [-1,-2,-3]
        vrr = vr[ir].mean(); vrl = vr[il].mean()
        stokesr = stokes[ir].mean(); stokesl = stokes[il].mean()
        area[i] = N.trapz(1 - stokes/fit_cont[i](vr), dx=1.8)
        rescale.append(cont - (1 - stokes/fit_cont[i](vr))*(fac/area[i]))
    return area
