import numpy as N
import matplotlib.pyplot as P


def norm(vr, I, cont=0.9992, fac = 5.1):
    ir = [0,1,2]; il = [-1,-2,-3]
    vrr = vr[ir].mean(); vrl = vr[il].mean()
    stokesr = I[ir].mean(); stokesl = I[il].mean()
    fit_cont = N.poly1d(N.polyfit((vrr, vrl), (stokesr, stokesl), deg=1))
    area = N.trapz(1 - I/fit_cont(vr), dx=1.8)
    return cont - (1 - I/fit_cont(vr))*(fac/area)

def normHa(vr, I):
    ir = [0,1,2,3,4,5,6,7,8,9]; il = [-1,-2,-3,-4,-5,-6,-7,-8,-9]
    vrr = vr[ir].mean(); vrl = vr[il].mean()
    stokesr = I[ir].mean(); stokesl = I[il].mean()
    fit_cont = N.poly1d(N.polyfit((vrr, vrl), (stokesr, stokesl), deg=1))
    return I/fit_cont(vr)
