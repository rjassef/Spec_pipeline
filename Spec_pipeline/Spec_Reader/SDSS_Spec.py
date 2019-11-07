import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.constants import h,c
import re
import os

from .spectrum1d import read_fits_spectrum1d
from .Spec import Spec
from .rebin_spec import rebin_spec

#As there are too many different things to keep in minds, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

class SDSS_Spec(Spec):

    def __init__(self,_name,_zspec,_fits_files):
        super(SDSS_Spec,self).__init__(_name,_zspec,_fits_files)
        self.RT   = 1.25*u.m #Telescope radius.
        self.instrument = "SDSS"
        self.__flam

    @property
    def __flam(self):

        s = fits.open(self.data_prefix+"/"+self.fits_files[0])
        self.lam_obs   = 10.**(s[1].data['loglam']) * u.AA
        self.flam      = s[1].data['flux'] * 1e-17 * u.erg/(u.cm**2*u.s*u.AA)
        self._flam_err = s[1].data['ivar']**-0.5 * 1e-17 * u.erg/(u.cm**2*u.s*u.AA)
        return

    @property
    def flam_err(self):
        return self._flam_err.to(u.erg/(u.cm**2*u.s*u.AA))
