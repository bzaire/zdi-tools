def zdipot_in(filename, outputname, incl, vsin, **kwargs):
    """ Creates the zdipot.in with the arguments informed.
        Some parameters could be variables to be passed in
        the bash routine.
        Parameters:
        incl, vsin, vrad, beta, gamma
    """
    f = InputNml(filename, incl, vsin, **kwargs)
    with open(outputname, 'w') as file:
        file.write('#!/bin/bash'+'\n')
        file.write('zdipot -b << aa'+'\n')
        # Do you want to make a magnetic image (y/n)
        file.write(f.mag_img+'\n')
        # Do you want to read an old image file (y/n)
        file.write(f.old_mag_img+'\n')
        # Input ngrid, incl, vsin(i) and vrad
        file.write(f.ngrid+' '+f.incl+' '+f.vsin+' '+f.vrad+'\n')
        # Select type of image reconstruction (any combination of QMB)
        file.write(f.QMB+'\n')
        # Type of reconstruction
        file.write(f.mag_type+'\n')
        # Order of spherical harmonics expansion
        file.write(f.lmax+'\n')
        # Input default spot brightness
        file.write(f.spot_brightness+'\n')
        # What is a typical value for the magnetic field strength
        file.write(f.tipic_B_strenght+'\n')
        # Beta and gamma for differential rotation
        file.write(f.beta+' '+f.gamma+'\n')
        # Name of file to read spectral data from
        file.write(f.filename+'\n')
        # Are profiles affected by moon pollution (y/n)
        file.write(f.moon_pol+'\n')
        # Do you want to use the Mean observations (y/n)
        file.write(f.mean_obs+'\n')
        # Select Stokes parameters to use (any combination of IV)
        file.write(f.stokes+'\n')
        # What is the largest phase smearing allowed
        file.write(f.phase+'\n')
        # What is the oversampling factor
        file.write(f.oversampling+'\n')
        # Do you want to weight the spherical harmonics (y/n)
        file.write(f.weight+'\n')
        # Enter value for L_fac
        file.write(f.lfac+'\n')
        # Input values for Caim and Maxit
        file.write(f.caim+' '+f.max_int+'\n')
        # Name of file to save brightness data to
        file.write(f.filebrightness+'\n')
        # Do you want to save the synthetic spectra (y/n)
        file.write(f.save_syntetic_spec+'\n')
        # Name of file to save spectral data to
        file.write(f.filesyntspec+'\n')
        if len(f.QMB)!=1 :
            # Do you want to save the SHaDe coefficients (y/n)
            file.write(f.save_SHaDeCoeff+'\n')
            # Name of file to save mode coefficients to
            file.write(f.fileSHaDeCoeff+'\n')
        file.write('aa')

class InputNml(object):
    "Set the input namelist for the ZDI run"
    def __init__(self, filename, incl, vsin, vrad=0, QMB='Q', beta=0, gamma=0, stokes='IV', save_syntetic_spec = 'y'):
        # Define the basic configuration, recieve additional information
        self.mag_img = 'y'
        self.old_mag_img = 'n'
        self.ngrid = '1000'
        self.incl = str(incl)
        self.vsin = str(vsin)
        self.vrad = str(vrad)
        self.QMB = QMB
        self.mag_type = '-4'
        self.lmax = '15'
        self.spot_brightness = '1'
        self.tipic_B_strenght = '1000'
        self.beta = str(beta)
        self.gamma = str(gamma)
        self.filename = filename
        self.moon_pol = 'n'
        self.mean_obs = 'y'
        self.stokes = stokes
        self.phase = '1'
        self.oversampling = '1'
        self.weight = 'y'
        self.lfac = '0'
        self.caim = '1'
        self.max_int = '100'
        self.filebrightness = self.filename[:-2]+"m1"
        self.save_syntetic_spec = str(save_syntetic_spec)
        self.filesyntspec = filename[:-2]+'s1'
        if len(self.QMB)!=1 :
            self.save_SHaDeCoeff = 'y'
            self.fileSHaDeCoeff = self.filename[:-2]+"c1"

        # Creates a zdipit.in file with the input configuration
#        self.output_filename = "zdipot_"+self.filename[:-3]
