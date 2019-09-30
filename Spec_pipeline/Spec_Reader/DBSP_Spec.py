#!/usr/bin/env python 

import numpy as np
#from specutils.io.read_fits import read_fits_spectrum1d
from astropy.io import fits
import astropy.units as u
from astropy.constants import h,c
import os

from .spectrum1d import read_fits_spectrum1d
from .Spec import Spec
from .rebin_spec import rebin_spec

#As there are too many different things to keep in minds, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

class DBSP_Spec(Spec):

    def __init__(self,_name,_zspec,_fits_files,_line_center=None,
                 blue=False,red=False):
        super(DBSP_Spec,self).__init__(_name,_zspec,_fits_files,_line_center)
        self.RT   = 2.5*u.m #Telescope radius.
        self.instrument = "DBSP"
        self.blue = blue
        self.red  = red
        self.__flam
        self.__flam_sky

    @property
    def __flam(self):
        
        spec_b = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_files[0],
                                      dispersion_unit=u.AA, 
                                      flux_unit = u.erg/u.cm**2/u.s/u.Hz)
        spec_r = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_files[1],
                                      dispersion_unit=u.AA, 
                                      flux_unit = u.erg/u.cm**2/u.s/u.Hz)

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
        elif self.red:
            self.lam_obs = spec_r[0].dispersion
            fnu = spec_r[0].data*spec_r[0].unit
            ff = fits.open(self.data_prefix+"/"+self.fits_files[1])
            self.spec_err_name = "error."+self.fits_files[1]

        self.dlam = np.mean(self.lam_obs[1:]-self.lam_obs[:-1])#Mean lambda bin.
        self.texp = float(ff[0].header['EXPTIME'])*u.s
        self.RON  = float(ff[0].header['RON'])

        self.flam = (fnu*c/self.lam_obs**2).to(u.erg/u.cm**2/u.s/u.AA)
        return

    @property
    def __flam_sky(self):

        #Figure out the spectrograph arm.
        if self.blue:
            sky_temp_fname = os.environ['SPEC_PIPE_LOC']+\
                "/Spec_pipeline/Sky_Templates/template_sky_DBSP_b.dat"
        elif self.red:
            sky_temp_fname = os.environ['SPEC_PIPE_LOC']+\
                "/Spec_pipeline/Sky_Templates/template_sky_DBSP_r.dat"
        else:
            #print("Cannot find spectrograph arm flag")
            return
            
        #Read the template
        sky_temp = np.loadtxt(sky_temp_fname)
        lam_sky = sky_temp[:,0]*u.AA
        flam_sky_orig = sky_temp[:,1]*u.erg/u.s/u.cm**2/u.AA

        #Rebin the template to the object spectrum.
        self.flam_sky = rebin_spec(lam_sky, flam_sky_orig, self.lam_obs)

        return

    @property
    def __sens(self):

        #Read the sensitivity curve.
        if self.blue:
            grname = "B"
        else:
            grname = "R"
        sens_temp = np.loadtxt(os.environ['SPEC_PIPE_LOC']+\
                               "/Spec_pipeline/Sensitivity_Files/"+
                               "Sens_DBSP_"+grname+".txt")
        lam_sens = sens_temp[:,0]*u.AA
        sens_orig = sens_temp[:,1]*u.dimensionless_unscaled
        
        #Rebin the template to the object spectrum.
        self.sens = rebin_spec(lam_sens, sens_orig, self.lam_obs)

        return
