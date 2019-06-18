import numpy as N
import matplotlib.pyplot as P
import cmocean
# Read txt file generated when searching the best parameters

def read_spi(*args, **kwargs):
    obj = SearchParamsi(*args, **kwargs)
    obj.__get_data__(*args)
    return obj
def read_spv(*args, **kwargs):
    obj = SearchParamsv(*args, **kwargs)
    obj.__get_data__(*args)
    return obj
def read_spdr(*args, **kwargs):
    obj = SearchParamsDR(*args, **kwargs)
    obj.__get_data__(*args)
    return obj
def read_spt(*args, **kwargs):
    obj = SearchParamsThree(*args, **kwargs)
    obj.__get_data__(*args)
    return obj

class SearchParamsDR(object):
    def __init__(self, filename):
        with open(filename, 'r') as file:
            line = file.readline().split()
            self.incl = float(line[3]); self.vsin = float(line[7])
        self.beta = None
        self.gamma = None
        self.chi = None
        self.s = None
        self.sp_ph = None
        self.test = None
        self.cool = None
        self.hot = None
        self.bbeta = None
        self.bgamma = None


    def __get_data__(self, filename):
        data = N.genfromtxt(filename, skip_header = 2)
        self.beta = data[:,0]
        self.gamma = data[:,1]
        self.chi = data[:,2]
        self.s = data[:,3]
        self.sp_ph = data[:,4]
        self.test = data[:,5]
        self.cool = data[:,6]
        self.hot = data[:,7]

    def plotCHI2(self, cmax = None):
        if cmax != None:
            cmax = self.chi < cmax
        P.scatter(self.beta[cmax], self.gamma[cmax], c=self.chi[cmax], cmap=cmocean.cm.thermal)
        P.xlim((self.beta.min(), self.beta.max()))
        P.ylim((self.gamma.min(), self.gamma.max()))
        P.xlabel(r'$\beta$')
        P.ylabel(r'$\gamma$')
        P.colorbar()
        P.tight_layout()
        P.show()

    def plotTest(self, tmax = None):
        if tmax != None:
            tmax = self.test < tmax
        P.scatter(self.beta[tmax], self.gamma[tmax], c=self.test[tmax], cmap=cmocean.cm.thermal)
        P.xlim((self.beta.min(), self.beta.max()))
        P.ylim((self.gamma.min(), self.gamma.max()))
        P.xlabel(r'$\beta$')
        P.ylabel(r'$\gamma$')
        P.colorbar()
        P.tight_layout()
        P.show()

    def plotSp(self, cmax = 1.0003, tmax = 0.003):
        if cmax != None and tmax != None:
            icmax = self.chi < cmax
            itmax = self.test < tmax
            imax =[icmax[i] and itmax[i] for i in range(self.chi.shape[0])]     
        elif cmax != None:
            imax = self.chi < cmax
        elif tmax != None:
            imax = self.test < tmax
        P.scatter(self.beta[imax], self.gamma[imax], c=self.sp_ph[imax], cmap=cmocean.cm.thermal_r)
        P.xlim((self.beta.min(), self.beta.max()))
        P.ylim((self.gamma.min(), self.gamma.max()))
        P.xlabel(r'$\beta$')
        P.ylabel(r'$\gamma$')
        P.colorbar()
        P.tight_layout()
        P.show()

    def plotS(self, cmax = 1.0003, tmax = 0.003):
        if cmax != None and tmax != None:
            icmax = self.chi < cmax
            itmax = self.test < tmax
            imax =[icmax[i] and itmax[i] for i in range(self.chi.shape[0])]     
        elif cmax != None:
            imax = self.chi < cmax
        elif tmax != None:
            imax = self.test < tmax
        P.scatter(self.beta[imax], self.gamma[imax], c=self.s[imax], cmap=cmocean.cm.thermal)
        P.xlim((self.beta.min(), self.beta.max()))
        P.ylim((self.gamma.min(), self.gamma.max()))
        P.xlabel(r'$\beta$')
        P.ylabel(r'$\gamma$')
        P.colorbar()
        P.tight_layout()
        P.show()

    def plot3D(self, field, smin = None, tmax = 0.003, cmax = 1.0003):
        from mpl_toolkits.mplot3d import Axes3D
        fig = P.figure()
        ax = fig.add_subplot(111, projection='3d')
        icmax = self.chi < cmax
        itmax = self.test < tmax
        icb = N.bitwise_and(icmax, itmax)
        beta = self.beta[icb].copy(); gamma = self.gamma[icb].copy(); s = self.s[icb].copy(); test = self.test[icb].copy(); field = field[icb]
        if smin != None:
            isc = self.s > smin
            beta = beta[isc]; gamma = gamma[isc]; s = s[isc]; test = test[isc]; field = field[isc]
        ax.scatter(beta, gamma, field, c=field, cmap=cmocean.cm.thermal_r)
        P.xlabel(r'$\beta$')
        P.ylabel(r'$\gamma$')
        P.show()

    def PolFitS(self, cmax, tmax, smin):
        icmax = self.chi < cmax
        beta = self.beta[icmax].copy(); gamma = self.gamma[icmax].copy(); sp_ph = self.sp_ph[icmax].copy(); s = self.s[icmax].copy(); test = self.test[icmax].copy()
        itmax = test < tmax
        beta = beta[itmax].copy(); gamma = gamma[itmax].copy(); sp_ph = sp_ph[itmax].copy(); s = s[itmax].copy(); test = test[itmax]
        ismax = s > smin
        beta = beta[ismax].copy(); gamma = gamma[ismax].copy(); sp_ph = sp_ph[ismax].copy(); s = s[ismax].copy(); test = test[ismax].copy()

        #fit = N.polyfit(beta, s, deg = 2)
        fit = N.polyfit(beta, sp_ph, deg = 2)
        self.bbeta = - fit[1]/(2.*fit[0])
        #gfit = N.polyfit(gamma, s, deg = 2)
        gfit = N.polyfit(gamma, sp_ph, deg = 2)
        self.bgamma = - gfit[1]/(2.*gfit[0])

class SearchParamsThree(object):
    def __init__(self, filename):
        with open(filename, 'r') as file:
            line = file.readline().split()
            self.beta = float(line[3]); self.gamma = float(line[7])
        self.incl = None
        self.vsin = None
        self.vrad = None
        self.chi = None
        self.s = None
        self.sp_ph = None
        self.test = None
        self.cool = None
        self.hot = None


    def __get_data__(self, filename):
        data = N.genfromtxt(filename, skip_header = 2)
        self.incl = data[:,0]
        self.vsin = data[:,1]
        self.vrad = data[:,2]
        self.chi = data[:,3]
        self.s = data[:,4]
        self.sp_ph = data[:,5]
        self.test = data[:,6]
        self.cool = data[:,7]
        self.hot = data[:,8]

    def plotCHI2(self, cmax = None):
        if cmax != None:
            cmax = self.chi < cmax
        P.scatter(self.incl[cmax], self.vsin[cmax], c=self.chi[cmax], cmap=cmocean.cm.thermal)
        P.xlim((self.incl.min(), self.incl.max()))
        P.ylim((self.vsin.min(), self.vsin.max()))
        P.xlabel(r'$i$')
        P.ylabel(r'$v\sin(i)$')
        P.colorbar()
        P.tight_layout()
        P.show()

    def plotTest(self, tmax = None):
        if tmax != None:
            tmax = self.test < tmax
        P.scatter(self.beta[tmax], self.gamma[tmax], c=self.test[tmax], cmap=cmocean.cm.thermal)
        P.xlim((self.incl.min(), self.incl.max()))
        P.ylim((self.vsin.min(), self.vsin.max()))
        P.xlabel(r'$i$')
        P.ylabel(r'$v\sin(i)$')
        P.colorbar()
        P.tight_layout()
        P.show()

    def plot3D(self, cmax = None):
        from mpl_toolkits.mplot3d import Axes3D
        fig = P.figure()
        ax = fig.add_subplot(111, projection='3d')
        if cmax != None:
            cmax = self.chi < cmax
            ax.scatter(self.incl[cmax], self.vsin[cmax], self.chi[cmax], c=self.chi[cmax], cmap=cmocean.cm.thermal)
        else:
            ax.scatter(self.incl, self.vsin, self.chi, c=self.chi, cmap=cmocean.cm.thermal)
        P.show()

class SearchParamsi(object):
    def __init__(self, filename):
        with open(filename, 'r') as file:
            line1 = file.readline().split()
            line = file.readline().split()
            self.vsin = float(line[3]); self.vrad = float(line[7]); self.beta = float(line[11]); self.gamma = float(line[15])
        self.beta = None
        self.gamma = None
        self.chi = None
        self.s = None
        self.sp_ph = None
        self.test = None
        self.cool = None
        self.hot = None
        self.bincl = None

    def __get_data__(self, filename):
        data = N.genfromtxt(filename, skip_header = 3)
        self.incl = data[:,0]
        self.chi = data[:,1]
        self.s = data[:,2]
        self.sp_ph = data[:,3]
        self.test = data[:,4]
        self.cool = data[:,5]
        self.hot = data[:,6]

    def PolFit(self):
        fit = N.polyfit(self.incl, self.chi, deg = 2)
        self.bincl = - fit[1]/(2.*fit[0])
        inew = N.linspace(self.incl.min(),self.incl.max(),3*self.incl.shape[0]) 
        return inew, N.polyval(fit, inew)
        
    def plotCHI2(self, cmax = None):
        if cmax != None:
            cmax = self.chi < cmax
        P.scatter(self.incl[cmax], self.chi[cmax])
        P.xlim((self.incl.min(), self.incl.max()))
        P.ylim((self.chi.min(), self.chi.max()))
        P.xlabel(r'$i$')
        P.ylabel(r'$\chi^2$')
        P.tight_layout()
        P.show(False)

class SearchParamsv(object):
    def __init__(self, filename):
        with open(filename, 'r') as file:
            line1 = file.readline().split()
            line = file.readline().split()
            self.incl = float(line[3]); self.vrad = float(line[7]); self.beta = float(line[11]); self.gamma = float(line[15])
        self.vsin = None
        self.chi = None
        self.s = None
        self.sp_ph = None
        self.test = None
        self.cool = None
        self.hot = None
        self.bvsin = None

    def __get_data__(self, filename):
        data = N.genfromtxt(filename, skip_header = 3)
        self.vsin = data[:,0]
        self.chi = data[:,1]
        self.s = data[:,2]
        self.sp_ph = data[:,3]
        self.test = data[:,4]
        self.cool = data[:,5]
        self.hot = data[:,6]

    def PolFit(self):
        isp = self.s > -50.
        fit = N.polyfit(self.vsin[isp], self.s[isp], deg = 2)
        self.bvsin = - fit[1]/(2.*fit[0])
        inew = N.linspace(self.vsin.min(),self.vsin.max(),3*self.vsin.shape[0]) 
        return inew, N.polyval(fit, inew)
        
    def plotCHI2(self, cmax = None):
        if cmax != None:
            cmax = self.chi < cmax
        P.scatter(self.vsin[cmax], self.chi[cmax])
        P.xlim((self.vsin.min(), self.vsin.max()))
        P.ylim((self.chi.min(), self.chi.max()))
        P.xlabel(r'$v \sin(i)$')
        P.ylabel(r'$\chi^2$')
        P.tight_layout()
        P.show(False)
