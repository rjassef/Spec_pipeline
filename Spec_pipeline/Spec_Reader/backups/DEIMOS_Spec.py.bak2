#!/usr/bin/env python

import numpy as np
#from specutils.io.read_fits import read_fits_spectrum1d
from astropy.io import fits
import astropy.units as u
from astropy.constants import h,c
import re
import os

#from .spectrum1d import read_fits_spectrum1d
from .iraf_spectrum1d import read_fits_spectrum1d
from .Spec import Spec
from .rebin_spec import rebin_spec

#As there are too many different things to keep in mind, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

class DEIMOS_Spec(Spec):

    """
Module that read an DEIMOS spectrum and returns a spec object.

Args:
   name (string)      : Object name or ID.

   zspec (float)      : Spectroscopic redshift.

   fits_files (list)  : Spectrum file name. Has to be a one element list.

   line_center (float): Optional. Grating used for the observations.
                         Needs to have astropy units of AA.

   blue (boolean)      : Optional. If wavenlength of interest is not provided,
                         it should be indicated whether the blue side or the
                         red side spectrum should be loaded.

   red (boolean)       : Optional. If wavenlength of interest is not provided,
                         it should be indicated whether the blue side or the
                         red side spectrum should be loaded.

   show_err_plot (boolean) : Show plot of the error fit if it is carried out.

   local_sky_files (list)     : Optional. List of sky files if the default ones
                                are not to be used.

   local_sens_files (list)    : Optional. List of sensitivity files if the
                                default ones are not to be used.

   inst_conf                  : Optional. Configurations dictionary.

   header_kws                 : Optional. Dictionary. Defaults are:
                                "detector"  : "CCDGEOM",
                                "apsize_pix": "apsize_pix",
                                "texp"      : "EXPTIME",

   """

    def __init__(self,name,zspec,fits_files,line_center=None,
                 blue=False,red=False,show_err_plot=False,local_sky_files=None,local_sens_files=None, inst_conf=None, header_kws=None):

        super(DEIMOS_Spec,self).__init__(name, zspec, fits_files, line_center, show_err_plot=show_err_plot, local_sky_files=local_sky_files, local_sens_files=local_sens_files)

        self.RT   = 5.0*u.m #Telescope radius.
        self.instrument = "DEIMOS"
        self.blue = blue
        self.red = red

        #Deimos is not a dual spectrograph really, but we treat it as such since Dan reduced them with a dual spectrograph scheme.
        self.dual_spec = True
        self.edge_drop = 0.*u.AA
        self.dichroic  = "None"
        self.dichroic_wave = None

        #Parse the configuration.
        if inst_conf is not None:
            for kw in inst_conf.keys():
                kwuse = kw
                if kwuse[:4]=='blue' and self.blue:
                    kwuse = kw[5:]
                elif kwuse[:3]=='red' and self.red:
                    kwuse = kw[4:]
                setattr(self,kwuse,inst_conf[kw])

        self.__flam
        self._Spec__flam_sky
        self._Spec__sens


    @property
    def __flam(self):

        spec_b = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_files[0],
                                      dispersion_unit=u.AA,
                                      flux_unit = u.erg/(u.cm**2*u.s*u.Hz))
        spec_r = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_files[1],
                                      dispersion_unit=u.AA,
                                      flux_unit = u.erg/(u.cm**2*u.s*u.Hz))

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


        #Load attributes from header keywords. Set the default one, and overwrite them with, or add to them, the ones set by the user.
        keywords_to_load = {
            "detector"  : "CCDGEOM",
            "apsize_pix": "apsize_pix",
            "texp"      : "EXPTIME",
        }
        if header_kws is not None:
            for kw in header_kws.keys():
                keywords_to_load[kw] = header_kws[kw]
        self.load_keyword_headers(spec_use, keywords_to_load)

        #Detectors
        if self.detector is not None:
            if re.search("MIT/LL",self.detector):
                self.detector = "MIT_LL"
            else:
                pass

        #Slit width
        if self.slit_width is not None:
            try:
                self.slit_width.unit
            except AttributeError:
                m = re.match("long_(.*)",self.slit_width)
                self.slit_width = float(m.group(1)) * u.arcsec

        #Since DEIMOS has the pretty uncommon feature of using filters to control the wavelength range, we should include the filter in the sensitivity curve. Since this scapes the name convention in Spec.py, we'll overload it as a local sens file if none have been given.
        if self.local_sens_files is None:
            sens_temp_fname = os.environ['SPEC_PIPE_LOC'] + "/Spec_pipeline/Sensitivity_Files/" + "Sens_{0:s}_{1:s}_{2:s}_{3:s}.txt".format(self.instrument, self.detector, self.grating, self.filter)
            self.local_sens_files = [sens_temp_fname, sens_temp_fname]

        #Finish the setup
        self.run_setup(spec_use)

        return
