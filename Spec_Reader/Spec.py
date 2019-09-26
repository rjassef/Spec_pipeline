#As there are too many different things to keep in minds, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

import numpy as np
import astropy.units as u
from astropy.constants import h,c
from obtain_error_spectrum import get_error_spec

class Spec(object):

    def __init__(self,_name,_zspec,_fits_files,_line_center):
        self.RT   = None
        self.instrument = None
        self.name  = _name
        self.zspec = _zspec
        self.fits_files = _fits_files
        self.line_center = _line_center
        self.lam_obs  = None
        self.dlam = None
        self.texp = None
        self.RON  = None
        self.spec_err_name = None
        self.K1 = None
        self.K2 = None
        self.data_prefix = "data/"

    @property
    def lam_rest(self):
        try:
            return self._lam_rest
        except AttributeError:
            pass
        self._lam_rest = self.lam_obs/(1.+self.zspec)
        return self._lam_rest
    
    @property
    def eps(self):
        try:
            return self._eps
        except AttributeError:
            pass
        self._eps = (np.pi*self.RT**2) * self.dlam * self.texp/(h*c)
        self._eps = self._eps.to(u.s*u.cm**2/u.erg)
        return self._eps

    @property
    def flam_err(self):
        try:
            return self._flam_err
        except AttributeError:
            pass
        #Try to open the error spectrum. If not found, generate it.
        try:
            cat = open(self.spec_err_name,"r")
            self._flam_err = np.loadtxt(cat,usecols=[1])
            self._flam_err = self._flam_err * u.erg/u.s/u.cm**2/u.AA
        except IOError:
            self._flam_err, self.K1, self.K2 = \
                                               get_error_spec(self,wd=15)
            np.savetxt(self.data_prefix+"/"+self.spec_err_name,
                       np.array([self.lam_obs,
                                 self._flam_err]).T)
        return self._flam_err


