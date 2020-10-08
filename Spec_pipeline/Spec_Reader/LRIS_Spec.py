#!/usr/bin/env python

import numpy as np
import astropy.units as u
import re

from .iraf_spectrum1d import read_fits_spectrum1d
from .Spec import Spec


#As there are too many different things to keep in mind, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

class LRIS_Spec(Spec):
    """

    Module that reads an LRIS spectrum and returns a spec object.

    Parameters
    ----------
    name : string
        Object name or ID.

    zspec : float
        Spectroscopic redshift.

    fits_file : string
        Spectrum file name.

    blue : boolean, optional
        Indicated the provided spectrum is from the blue arm of the spectrograph. Either blue or red must be provided.

    red : boolean, optional
        Indicated the provided spectrum is from the red arm of the spectrograph. Either blue or red must be provided.

    show_err_plot : boolean, optional
        True if error-fit plot is to be displayed.

    local_sky_file : string, optional
        Sky file if the default ones are not to be used.

    local_sens_file : string, optional
        Sensitivity file if the default ones are not to be used.

    inst_conf : dictionary, optional
        Configurations dictionary.

    header_kws : dictionary, optional
        Dictionary with header keywords to use. Have precedence over default header keywords.

    """

    def __init__(self, name, zspec, fits_file, blue=False, red=False, show_err_plot=False, local_sky_file=None, local_sens_file=None, inst_conf=None, header_kws=None):

        #Load basic instrument properties.
        RT   = 5.0*u.m #Telescope radius.
        instrument = "LRIS"
        self.dual_spec = True
        self.blue = blue
        self.red = red
        self.edge_drop = 75.*u.AA

        super(LRIS_Spec,self).__init__(name, zspec, fits_file, show_err_plot=show_err_plot, local_sky_file=local_sky_file, local_sens_file=local_sens_file, inst_conf=inst_conf, header_kws=header_kws, RT=RT, instrument=instrument)

        #Finally, load the spectra.
        self.__flam
        self._Spec__flam_sky
        self._Spec__sens

        return

    @property
    def __flam(self):

        #Load the spectrum
        spec = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_file, dispersion_unit=u.AA, flux_unit = u.erg/(u.cm**2*u.s*u.Hz))

        #Assign the error name file.
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
        self.run_setup(spec)

        return
