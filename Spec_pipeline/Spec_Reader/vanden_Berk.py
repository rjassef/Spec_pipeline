import numpy as np
from astropy.table import Table
import astropy.units as u
from astropy.constants import h,c
import re
import os

from .Spec import Spec
from .rebin_spec import rebin_spec

class vandenBerk_Spec(Spec):

    """
Module that reads the vanden Berk composite spectrum and returns a spec object.
The composite is loaded as a spectrum at redshift _zspec.

Args:

   _zspec (float)      : Spectroscopic redshift.

   _fits_files (list)  : Spectrum file name. Has to be a one element list.

   """

    def __init__(self,_zspec,):
        super(vandenBerk_Spec,self).__init__("vB01_template",_zspec)
        self.RT   = 1.25*u.m #Telescope radius.
        self.instrument = "SDSS"
        self.__flam

    @property
    def __flam(self):

        sfile = os.environ['SPEC_PIPE_LOC']+\
                "/Spec_pipeline/Spec_Reader/vandenberk_composite.txt"
        qso = Table.read(sfile,format='ascii.cds')

        lam_rest_full  = qso['Wave'].to(u.AA)
        self.lam_obs = lam_rest_full*(1.+self.zspec)

        flux_unit      = 1e-17 * u.erg/(u.cm**2*u.s*u.AA)
        self.flam      = qso['FluxD'] * flux_unit
        self._flam_err = qso['e_FluxD'] * flux_unit

        return

    @property
    def flam_err(self):
        return self._flam_err.to(u.erg/(u.cm**2*u.s*u.AA))
