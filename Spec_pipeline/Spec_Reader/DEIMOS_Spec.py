#!/usr/bin/env python

import numpy as np
import astropy.units as u
import re
import os

from .iraf_spectrum1d import read_fits_spectrum1d
from .Spec import Spec

#As there are too many different things to keep in mind, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

class DEIMOS_Spec(Spec):

    """
Module that read an DEIMOS spectrum and returns a spec object.

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

   header_kws                 : Optional. Dictionary.

   """

    def __init__(self, name, zspec, fits_file, blue=False, red=False, show_err_plot=False, local_sky_file=None, local_sens_file=None, inst_conf=None, header_kws=None):

        #Set basic instrument properties.
        RT   = 5.0*u.m #Telescope radius.
        instrument = "DEIMOS"
        self.blue = blue
        self.red = red
        self.dual_spec = True

        iinit = super(DEIMOS_Spec,self).__init__(name, zspec, fits_file, show_err_plot=show_err_plot, local_sky_file=local_sky_file, local_sens_file=local_sens_file, inst_conf=inst_conf, header_kws=header_kws, RT=RT, instrument=instrument)

        #Problem loading super class. Do not continue.
        if iinit==1:
            return

        #Deimos is not a dual spectrograph really, but we treat it as such since Dan reduced them with a dual spectrograph scheme.
        self.edge_drop = 0.*u.AA
        self.dichroic  = "None"
        self.dichroic_wave = None

        #Load the spectra
        self.__flam
        self._Spec__flam_sky
        self._Spec__sens


    @property
    def __flam(self):

        #Load the spectrum
        spec = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_file, dispersion_unit=u.AA, flux_unit = u.erg/(u.cm**2*u.s*u.Hz))

        #Finally, assign the error name file.
        self.spec_err_name = "error."+self.fits_file

        #Set the channel.
        if self.blue:
            self.channel = 'b'
        elif self.red:
            self.channel = 'r'

        #Load attributes from header keywords.
        self.load_keyword_headers(spec, self.keywords_to_load)

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
        if self.local_sens_file is None:
            self.local_sens_file = os.environ['SPEC_PIPE_LOC'] + "/Spec_pipeline/Sensitivity_Files/" + "Sens_{0:s}_{1:s}_{2:s}_{3:s}.txt".format(self.instrument, self.detector, self.grating, self.filter)

        #Finish the setup
        self.run_setup(spec)

        return
