#!/usr/bin/env python 

import numpy as np
#from specutils.io.read_fits import read_fits_spectrum1d
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

class GMOS_Spec(Spec):

    def __init__(self,_name,_zspec,_fits_files,_line_center=None):
        super(GMOS_Spec,self).__init__(_name,_zspec,_fits_files,_line_center)
        self.RT   = 4.0*u.m #Telescope radius.
        self.instrument = "GMOS"
        self.__flam
        self.__flam_sky
        self.__sens

    @property
    def __flam(self):

        ff = fits.open(self.data_prefix+"/"+self.fits_files[0])
        spec = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_files[0],
                                    dispersion_unit=u.AA, 
                                    flux_unit = u.erg/u.cm**2/u.s/u.Hz)
        self.lam_obs = spec[0].dispersion
        fnu = spec[0].data * spec[0].unit

        self.dlam = np.mean(self.lam_obs[1:]-self.lam_obs[:-1])#Mean lambda bin.
        self.texp = ff[0].header['EXPTIME']*u.s
        self.RON  = ff[0].header['RDNOISE']
        self.spec_err_name = "error."+self.fits_files[0]
        self.spec_err_name = re.sub(".fits",".txt",self.spec_err_name)

        self.flam = (fnu*c/self.lam_obs**2).to(u.erg/u.cm**2/u.s/u.AA)
        return 

    @property
    def __flam_sky(self):
        
        #Read the template
        sky_temp = np.loadtxt(os.environ['SPEC_PIPE_LOC']+\
                              "/Spec_pipeline/Sky_Templates/template_sky_GMOS.dat")
        lam_sky = sky_temp[:,0]*u.AA
        flam_sky_orig = sky_temp[:,1]*u.erg/u.s/u.cm**2/u.AA

        #Rebin the template to the object spectrum.
        self.flam_sky = rebin_spec(lam_sky, flam_sky_orig, self.lam_obs)

        return

    #There does not seem to be a way to recover the GRISM used for the
    #GMOS spectra. We'll assume the B600 for all until we figure
    #something out.
    @property
    def __sens(self):

        #Read the sensitivity curve.
        grname =  "B600"
        sens_temp = np.loadtxt(os.environ['SPEC_PIPE_LOC']+\
                               "/Spec_pipeline/Sensitivity_Files/"+
                               "Sens_GMOS_"+grname+".txt")
        lam_sens = sens_temp[:,0]*u.AA
        sens_orig = sens_temp[:,1]*u.dimensionless_unscaled
        
        #Rebin the template to the object spectrum.
        self.sens = rebin_spec(lam_sens, sens_orig, self.lam_obs)

        return

        
