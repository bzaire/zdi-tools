import lmfit
import numpy as N
from .norm import norm 

# create gaussian model
gmodel = lmfit.models.GaussianModel()  

def GetVrad(vr, I):
    J = 1./I.copy() - 1.
    #indx = N.bitwise_or(vr<-70., vr>70.)
    #result = gmodel.fit(data=J[indx], x=vr[indx], amplitude=0.05, center=0, sigma=30)
    result = gmodel.fit(data=J, x=vr, amplitude=0.05, center=0, sigma=30)
    vrad = result.best_values['center']
    return vrad

def GetVrad1(vr, I):
    indx = N.bitwise_or(vr>-90., vr<90.)
    normI = norm(vr[indx], I[indx], cont=0.9993, fac=5.13) # norm data
    J = 1. - normI 
    area = N.trapz(J, dx=1.8)
    warea = N.trapz(J*vr[indx], dx=1.8)
    vrad = warea/area
    return vrad

