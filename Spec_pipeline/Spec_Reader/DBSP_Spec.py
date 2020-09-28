#!/usr/bin/env python

import numpy as np
import astropy.units as u
import re

from .iraf_spectrum1d import read_fits_spectrum1d
from .Spec import Spec

#As there are too many different things to keep in minds, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

class DBSP_Spec(Spec):

    """
Module that reads a DBSP spectrum and returns a spec object.

Args:
   name (string)              : Object name or ID.

   zspec (float)              : Spectroscopic redshift.

   fits_file (string)         : Spectrum file name.

   blue (boolean)             : Optional. Indicated the provided spectrum is
                                from the blue arm of the spectrograph. Either
                                blue or red must be provided.

   red (boolean)              : Optional. Indicated the provided spectrum is
                                from the red arm of the spectrograph. Either
                                blue or red must be provided.

   show_err_plot (boolean)    : Optional. True if error-fit plot is to be
                                displayed.

   local_sky_file (string)    : Optional. Sky file if the default ones
                                are not to be used.

   local_sens_file (string)   : Optional. Sensitivity file if the
                                default ones are not to be used.

   inst_conf (dict)           : Optional. Configurations dictionary.

   header_kws (dict)          : Optional. Dictionary.

   """

    def __init__(self, name, zspec, fits_file, blue=False, red=False, show_err_plot=False, local_sky_file=None, local_sens_file=None, inst_conf=None, header_kws=None):

        #Set basic instrument properties.
        RT   = 2.5*u.m #Telescope radius.
        instrument = "DBSP"
        self.dual_spec = True
        self.blue = blue
        self.red  = red
        self.edge_drop = 75.*u.AA

        #Load the super class.
        iinit = super(DBSP_Spec,self).__init__(name, zspec, fits_file, show_err_plot=show_err_plot, local_sky_file=local_sky_file, local_sens_file=local_sens_file, inst_conf=inst_conf, header_kws=header_kws, RT=RT, instrument=instrument)

        #Problem loading super class. Do not continue.
        if iinit == 1:
            return

        #Load the spectra
        self.__flam
        self._Spec__flam_sky
        self._Spec__sens

        return

    @property
    def __flam(self):

        #Load the spectrum.
        spec = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_file, dispersion_unit=u.AA, flux_unit = u.erg/(u.cm**2*u.s*u.Hz))

        #Set the name of the error file.
        self.spec_err_name = "error."+self.fits_file

        #Set the channel.
        if self.blue:
            self.channel = 'b'
        elif self.red:
            self.channel = 'r'

        #Load attributes from header keywords.
        self.load_keyword_headers(spec, self.keywords_to_load)

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
            try:
                self.slit_width.unit
            except AttributeError:
                self.slit_width = float(self.slit_width) * u.arcsec

        #Finish the setup
        self.run_setup(spec)

        return
