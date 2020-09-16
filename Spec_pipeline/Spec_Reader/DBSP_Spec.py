#!/usr/bin/env python

import numpy as np
#from specutils.io.read_fits import read_fits_spectrum1d
from astropy.io import fits
import astropy.units as u
from astropy.constants import h,c
import re
import os
from scipy.interpolate import interp1d

#from .spectrum1d import read_fits_spectrum1d
from .iraf_spectrum1d import read_fits_spectrum1d
from .Spec import Spec
from .rebin_spec import rebin_spec

#As there are too many different things to keep in minds, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

class DBSP_Spec(Spec):

    """
Module that reads a DBSP spectrum and returns a spec object.

Args:
   _name (string)      : Object name or ID.

   _zspec (float)      : Spectroscopic redshift.

   _fits_files (list)  : Spectrum file name. Has to be a one element list.

   _line_center (float): Optional. Grating used for the observations.
                         Needs to have astropy units of AA.

   blue (boolean)      : Optional. If wavenlength of interest is not provided,
                         it should be indicated whether the blue side or the
                         red side spectrum should be loaded.

   red (boolean)       : Optional. If wavenlength of interest is not provided,
                         it should be indicated whether the blue side or the
                         red side spectrum should be loaded.

   show_err_plot (boolean)    : Optional. True if error-fit plot is to be
                                displayed.

   local_sky_files (list)     : Optional. List of sky files if the default ones
                                are not to be used.

   local_sens_files (list)    : Optional. List of sensitivity files if the
                                default ones are not to be used.

   dichroic (string)          : Optional.

   grating (string)           : Optional.

   grating_dispersion (float) : Optional. In astropy units of AA/mm

   detector (string)          : Optional.

   plate_scale (float)        : Optional. In astropy units of arcsec.

   pixel_size (float)         : Optional. In astropy units of microns.

   slit_width (float)         : Optional. In astropy units of arcsec.

   """

    def __init__(self,_name,_zspec,_fits_files,_line_center=None,
                 blue=False,red=False,show_err_plot=False,local_sky_files=None,local_sens_files=None, dichroic=None, grating=None, grating_dispersion=None, detector=None, plate_scale=None, pixel_size=None, slit_width=None):
        super(DBSP_Spec,self).__init__(_name,_zspec,_fits_files,_line_center,show_err_plot=show_err_plot)
        self.RT   = 2.5*u.m #Telescope radius.
        self.instrument = "DBSP"
        self.dual_spec = True
        self.blue = blue
        self.red  = red
        self.edge_drop = 75.*u.AA
        self.local_sky_files = local_sky_files
        self.local_sens_files = local_sens_files

        self.dichroic = dichroic
        self.grating = grating
        self.grating_dispersion = grating_dispersion
        self.detector = detector
        self.plate_scale = plate_scale
        self.pixel_size = pixel_size
        self.slit_width = slit_width

        self._sigma_res = None

        self.__flam
        self.__flam_sky
        self.__sens

    @property
    def __flam(self):

        spec_b = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_files[0],
                                      dispersion_unit=u.AA,
                                      flux_unit = u.erg/(u.cm**2*u.s*u.Hz))
        spec_r = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_files[1],
                                      dispersion_unit=u.AA,
                                      flux_unit = u.erg/(u.cm**2*u.s*u.Hz))

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
            #Assign the blue spectrum to be used
            spec_use = spec_b
            #Set the channel
            self.channel = 'b'
            #Finally, assign the error name file.
            self.spec_err_name = "error."+self.fits_files[0]

        elif self.red:
            #Assign the red spectrum to be used
            spec_use = spec_r
            #Set the channel
            self.channel = 'r'
            #Finally, assign the error name file.
            self.spec_err_name = "error."+self.fits_files[1]

        #Find some important aspects of the observations.
        #Dichroic
        keywords_to_load = {
            "dichroic"  : "DICHROIC",
            "detector"  : "DETNAM",
            "grating"   : "GRATING",
            "slit_width": "APERTURE",
            "apsize_pix": "apsize_pix",
            "texp"      : "EXPTIME",
            "RON"       : "RON",
            "GAIN"      : "GAIN"
        }
        self.load_keyword_headers(spec_use, keywords_to_load)

        #Dichroic
        if self.dichroic is not None:
            self.dichroic = re.sub("-","",self.dichroic)
            self.dichroic_wave = float(self.dichroic[1:])*100.*u.AA


        #Grating - there is a grating called sometimes 316/7150 and sometimes 316/7500, but it is the same grating with the same blaze as confirmed by the Palomar Observatory staff. We will alwayd use 316/7500.
        if self.grating is not None:
            self.grating  = re.sub("7150","7500",self.grating)
            self.grating  = re.sub("/","-",self.grating)
            if self.grating == "3167500":
                self.grating = "316-7500"

        #Slit width
        if self.slit_width is not None:
            self.slit_width = float(self.slit_width) * u.arcsec

        #Finish the setup
        self.run_setup(spec_use)

        return

    @property
    def __flam_sky(self):

        #Read the template
        sky_temp = np.loadtxt(self.sky_temp_fname)
        lam_sky = sky_temp[:,0]*u.AA
        flam_sky_orig = sky_temp[:,1]*u.erg/(u.s*u.cm**2*u.AA)

        #Rebin the template to the object spectrum.
        self.flam_sky = rebin_spec(lam_sky, flam_sky_orig, self.lam_obs)

        return

    @property
    def __sens(self):

        sens_temp = np.loadtxt(self.sens_temp_fname)
        lam_sens = sens_temp[:,0]*u.AA
        sens_orig = sens_temp[:,1]*u.dimensionless_unscaled

        #Rebin the template to the object spectrum.
        #self.sens = rebin_spec(lam_sens, sens_orig, self.lam_obs)

        #Interpolate the sensitivity template to the object spectrum. Extrapolate if needed, which is OK as it is a smooth function of wavelegnth for the most part.
        if np.min(self.lam_obs)<np.min(lam_sens) or np.max(self.lam_obs)>np.max(lam_sens):
            print("Warning: Extrapolating sensitivity curve to match spectral range {0:s}".format(self.name))
            print("Spec-range: {0:.1f} - {1:.2f}".format( np.min(self.lam_obs),np.max(self.lam_obs)))
            print("Sens-range: {0:.1f} - {1:.2f}".format(np.min(lam_sens),np.max(lam_sens)))
        f = interp1d(lam_sens, sens_orig, kind='linear', fill_value='extrapolate')
        self.sens = f(self.lam_obs)

        return

    #We use the numbers here: https://www.astro.caltech.edu/palomar/observer/200inchResources/dbspoverview.html#grating . We'll assume that the dispersion numbers quoted are the FWHM (like for LRIS) and we'll assume a 1" slit, the same as we did for LRIS, since really what matters here is the minimum between the seeing and the slit size, and we want to be conservative.
    @property
    def sigma_res(self):

        if self._sigma_res is not None:
            return self._sigma_res

        slit_size = 1.0*u.arcsec

        if self.grating_dispersion is not None:
            res = self.grating_dispersion
        else:
             print("Grating dispersion not set.")
             print("Using minimum of 1 AA/mm ")
             res = 1.0*u.AA/u.mm

        FWHM_res = (slit_size/self.plate_scale)*self.pixel_size * res
        self._sigma_res = (FWHM_res/(2.*(2.*np.log(2.))**0.5)).to(u.AA)
        return self._sigma_res
