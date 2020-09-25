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

#As there are too many different things to keep in mind, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

class LRIS_Spec(Spec):

    """
Module that read an LRIS spectrum and returns a spec object.

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

   show_err_plot (boolean) : Show plot of the error fit if it is carried out.

   local_sky_files (list)     : Optional. List of sky files if the default ones
                                are not to be used.

   local_sens_files (list)    : Optional. List of sensitivity files if the
                                default ones are not to be used.

   inst_conf                  : Optional. Configurations dictionary.
   """

    def __init__(self,_name,_zspec,_fits_files,_line_center=None,
                 blue=False,red=False,show_err_plot=False,local_sky_files=None,local_sens_files=None, inst_conf=None):

        super(LRIS_Spec,self).__init__(_name,_zspec,_fits_files,_line_center,show_err_plot=show_err_plot, local_sky_files=local_sky_files, local_sens_files=local_sens_files)

        self.RT   = 5.0*u.m #Telescope radius.
        self.instrument = "LRIS"
        self.dual_spec = True
        self.blue = blue
        self.red = red
        self.edge_drop = 75.*u.AA

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

    return

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

        #Find some important aspects of the observations.
        #Dichroic
        keywords_to_load = {
            "dichroic"  : "DICHNAME",
            "detector"  : "CCDGEOM",
            "slit_width": "SLITNAME",
            "apsize_pix": "apsize_pix",
            "texp"      : "EXPTIME",
        }
        if self.blue:
            keywords_to_load['grating'] = "GRISNAME"
            keywords_to_load['detector'] = "CCDGEOM"
        else:
            keywords_to_load['grating'] = "GRANAME"
            keywords_to_load['detector'] = "DETECTOR"
        self.load_keyword_headers(spec_use, keywords_to_load)

        #Detectors
        if self.detector is not None:
            if re.search("e2v",self.detector):
                self.detector = "e2v"
            elif re.search("LBNL", self.detector):
                self.detector = "LBNL"
            elif re.search("Mark 2", self.detector):
                self.detector = "Mark2"
            else:
                pass

        #Dichroic
        if self.dichroic is not None:
            self.dichroic_wave = float(self.dichroic)*u.nm

        #Grating
        if self.grating is not None:
            self.grating  = re.sub("/","-",self.grating)

        #Slit width
        if self.slit_width is not None:
            try:
                self.slit_width.unit
            except AttributeError:
                m = re.match("long_(.*)",self.slit_width)
                self.slit_width = float(m.group(1)) * u.arcsec

        #Finish the setup
        self.run_setup(spec_use)

        return
