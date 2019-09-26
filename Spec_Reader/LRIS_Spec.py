#!/usr/bin/env python 

import numpy as np
#from specutils.io.read_fits import read_fits_spectrum1d
from spectrum1d import read_fits_spectrum1d
from astropy.io import fits
import astropy.units as u
from astropy.constants import h,c

from Spec import Spec

#As there are too many different things to keep in minds, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

class LRIS_Spec(Spec):

    def __init__(self,_name,_zspec,_fits_files,_line_center=None,
                 blue=False,red=False):
        super(LRIS_Spec,self).__init__(_name,_zspec,_fits_files,_line_center)
        self.RT   = 5.0*u.m #Telescope radius.
        self.instrument = "LRIS"
        self.blue = blue
        self.red = red
        self.__flam

    @property
    def __flam(self):
        
        spec_b = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_files[0],
                                      dispersion_unit=u.AA, 
                                      flux_unit = u.erg/u.cm**2/u.s/u.Hz)
        spec_r = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_files[1],
                                      dispersion_unit=u.AA, 
                                      flux_unit = u.erg/u.cm**2/u.s/u.Hz)

        #Line that we want to fit unless side has been decided already.
        if not self.blue and not self.red:
            lam_targ = self.line_center*(1.+self.zspec)
            if lam_targ>spec_b[0].dispersion.min() and \
               lam_targ<spec_b[0].dispersion.max():
                self.blue = True
            elif lam_targ>spec_r[0].dispersion.min() and \
                 lam_targ<spec_r[0].dispersion.max():
                self.red = True
            else:
                print("Line not within spectral ranges")
                return
        
        if self.blue:
            self.lam_obs = spec_b[0].dispersion
            fnu = spec_b[0].data*spec_b[0].unit
            ff = fits.open(self.data_prefix+"/"+self.fits_files[0])
            self.spec_err_name = "error."+self.fits_files[0]

            #https://www2.keck.hawaii.edu/inst/lris/detectors.html
            #Use mean of amplifiers.
            self.RON = 3.82

        elif self.red:
            self.lam_obs = spec_r[0].dispersion
            fnu = spec_r[0].data*spec_r[0].unit
            ff = fits.open(self.data_prefix+"/"+self.fits_files[1])
            self.spec_err_name = "error."+self.fits_files[1]
            #https://www2.keck.hawaii.edu/inst/lris/detectors.html
            #Use mean of amplifiers.
            self.RON = 4.64

        self.dlam = np.mean(self.lam_obs[1:]-self.lam_obs[:-1])#Mean lambda bin.
        self.texp = float(ff[0].header['EXPTIME'])*u.s

        self.flam = (fnu*c/self.lam_obs**2).to(u.erg/u.cm**2/u.s/u.AA)
        return
