#!/usr/bin/env python

import numpy as np
import astropy.units as u

from .iraf_spectrum1d import read_fits_spectrum1d
from .Spec import Spec

#As there are too many different things to keep in minds, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

class GMOS_Spec(Spec):

    """
    Module that read a GMOS spectrum and returns a spec object.

    Parameters
    ----------

    name : string
        Object name or ID.

    zspec : float
        Spectroscopic redshift.

    fits_file : string
        Spectrum file name.

    show_err_plot : boolean, optional
        True if error-fit plot is to be displayed.

    local_sky_file : string, optional
        Sky file if the default ones are not to be used.

    local_sens_file : string, optional
        Sensitivity file if the default ones are not to be used.

    inst_conf : dictionary, optional
        Configurations dictionary.

    header_kws: dictionary, optional
        Dictionary with header keywords to use. Have precedence over default header keywords.

    """

    def __init__(self, name, zspec, fits_file, show_err_plot=False, local_sky_file=None, local_sens_file=None, inst_conf=None, header_kws=None):

        #Load basic properties.
        RT   = 4.0*u.m #Telescope radius.
        instrument = "GMOS"

        iinit = super(GMOS_Spec,self).__init__(name, zspec, fits_file, show_err_plot=show_err_plot, local_sky_file=local_sky_file, local_sens_file=local_sens_file, inst_conf=inst_conf, header_kws=header_kws, RT=RT, instrument=instrument)

        #Problem loading super class. Do not continue.
        if iinit == 1:
            return

        #Force the addition of the blue exponential to the error.
        self.force_blue_exp_err = True

        #Load the spectrum
        self.__flam
        self._Spec__flam_sky
        self._Spec__sens

        return

    @property
    def __flam(self):

        spec = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_file,
                                    dispersion_unit=u.AA,
                                    flux_unit = u.erg/(u.cm**2*u.s*u.Hz))
        self.spec_err_name = "error."+self.fits_file

        #Load attributes from header keywords.
        self.load_keyword_headers(spec, self.keywords_to_load)

        #Slit width
        if self.slit_width is not None:
            try:
                self.slit_width.unit
            except AttributeError:
                self.slit_width = float(self.slit_width) * u.arcsec

        #Finish the setup
        self.run_setup(spec)
