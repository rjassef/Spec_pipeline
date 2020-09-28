#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.constants import h,c
import re
import os

#from .spectrum1d import read_fits_spectrum1d
from .iraf_spectrum1d import read_fits_spectrum1d
from .Spec import Spec
from .rebin_spec import rebin_spec

#As there are too many different things to keep in minds, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

class GMOS_Spec(Spec):

    """
Module that read a GMOS spectrum and returns a spec object.

Args:
   name (string)    : Object name or ID.

   zspec (float)    : Spectroscopic redshift.

   fits_files (list): Spectrum file name. Has to be a one element list.

   local_sky_files (list)     : Optional. List of sky files if the default ones
                                are not to be used.

   local_sens_files (list)    : Optional. List of sensitivity files if the
                                default ones are not to be used.

   inst_conf                  : Optional. Configurations dictionary.

   header_kws                 : Optional. Dictionary. Defaults are:
                                "apsize_pix": "apsize_pix",
                                "texp"      : "EXPTIME",
                                "RON"       : "RDNOISE",
                                "GAIN"      : "GAINMULT"

   """

    def __init__(self,name,zspec,fits_files,show_err_plot=False,local_sky_files=None,local_sens_files=None, inst_conf=None, header_kws=None):

        super(GMOS_Spec,self).__init__(name,zspec,fits_files,show_err_plot=show_err_plot,  local_sky_files=local_sky_files, local_sens_files=local_sens_files)

        self.dual_spec = False
        self.RT   = 4.0*u.m #Telescope radius.
        self.instrument = "GMOS"

        #Force the addition of the blue exponential to the error.
        self.force_blue_exp_err = True

        #self.grname = _grname
        if inst_conf is not None:
            for kw in inst_conf.keys():
                kwuse = kw
                if kwuse[:4]=='blue' and self.blue:
                    kwuse = kw[5:]
                setattr(self,kwuse,inst_conf[kw])

        self.__flam
        self._Spec__flam_sky
        self._Spec__sens

        return

    @property
    def __flam(self):

        spec = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_files[0],
                                    dispersion_unit=u.AA,
                                    flux_unit = u.erg/(u.cm**2*u.s*u.Hz))
        self.spec_err_name = "error."+self.fits_files[0]

        #Load attributes from header keywords. Set the default ones, and overwrite them with, or add to them, the ones set by the user.
        keywords_to_load = {
            "apsize_pix": "apsize_pix",
            "texp"      : "EXPTIME",
            "RON"       : "RDNOISE",
            "GAIN"      : "GAINMULT"
        }
        if header_kws is not None:
            for kw in header_kws.keys():
                keywords_to_load[kw] = header_kws[kw]
        self.load_keyword_headers(spec, keywords_to_load)

        #Slit width
        if self.slit_width is not None:
            try:
                self.slit_width.unit
            except AttributeError:
                self.slit_width = float(self.slit_width) * u.arcsec

        #Finish the setup
        self.run_setup(spec)
