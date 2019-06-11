import lmfit
import numpy as N

# create gaussian model
gmodel = lmfit.models.GaussianModel()  

def GetVrad(vr, I):
    J = 1./I.copy() - 1.
    #indx = N.bitwise_or(vr<-70., vr>70.)
    #result = gmodel.fit(data=J[indx], x=vr[indx], amplitude=0.05, center=0, sigma=30)
    result = gmodel.fit(data=J, x=vr, amplitude=0.05, center=0, sigma=30)
    vrad = result.best_values['center']
    return vrad
