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
   _name (string)    : Object name or ID.

   _zspec (float)    : Spectroscopic redshift.

   _fits_files (list): Spectrum file name. Has to be a one element list.

   _grname (string)  : Grating used for the observations.

   """

    def __init__(self,_name,_zspec,_fits_files,_grname,show_err_plot=False):
        super(GMOS_Spec,self).__init__(_name,_zspec,_fits_files,show_err_plot=show_err_plot)
        self.grname = _grname
        self.RT   = 4.0*u.m #Telescope radius.
        self.instrument = "GMOS"
        self.__flam
        self.__flam_sky
        self.__sens


    @property
    def __flam(self):

        ff = fits.open(self.data_prefix+"/"+self.fits_files[0])
        spec = read_fits_spectrum1d(self.data_prefix+"/"+self.fits_files[0],
                                    dispersion_unit=u.AA,
                                    flux_unit = u.erg/(u.cm**2*u.s*u.Hz))

        #Find the grism and remove data outside the edges of the sensitivity curves.
        sens_temp = np.loadtxt(os.environ['SPEC_PIPE_LOC']+\
                            "/Spec_pipeline/Sensitivity_Files/"+
                            "Sens_GMOS_"+self.grname+".txt")
        lam_sens = sens_temp[:,0]*u.AA
        kuse = (spec[0].dispersion>np.min(lam_sens)) & \
                (spec[0].dispersion<np.max(lam_sens))

        #Finally, figure out the sky template edges and trim the spectrum to that limit.
        sky_temp = np.loadtxt(os.environ['SPEC_PIPE_LOC']+\
                              "/Spec_pipeline/Sky_Templates/template_sky_GMOS.dat")
        lam_sky = sky_temp[:,0]*u.AA
        kuse = (kuse) & (spec[0].dispersion>np.min(lam_sky)) & \
                (spec[0].dispersion<np.max(lam_sky))

        #Now, assign the wavelength and flux to the object.
        self.lam_obs = spec[0].dispersion[kuse]
        fnu = spec[0].data[kuse]*spec[0].unit

        # self.lam_obs = spec[0].dispersion
        # fnu = spec[0].data * spec[0].unit

        #Pixel scale with 2x2 binning EEV detector: https://www.gemini.edu/instrumentation/gmos/capability#Imaging
        self.PIXSIZE = 2  * 0.073 * u.arcsec

        #No slit information in the headers. We'll assume the default 1.25" size from the spec class. This should not matter as the extraction aperture should be there in all the GMOS headers.

        #If no apsize_pix read from headers, assume the slit size for the extraction aperture.
        if 'apsize_pix' in spec[0].header:
            self.apsize_pix = spec[0].header['apsize_pix']
        else:
            self.apsize_pix = (self.slit_width/self.PIXSIZE).to(1.).value


        self.dlam = np.mean(self.lam_obs[1:]-self.lam_obs[:-1])#Mean lambda bin.
        self.texp = ff[0].header['EXPTIME']*u.s
        self.RON  = ff[0].header['RDNOISE']
        self.GAIN = ff[0].header['GAINMULT']
        self.spec_err_name = "error."+self.fits_files[0]
        self.spec_err_name = re.sub(".fits",".txt",self.spec_err_name)

        self.flam = (fnu*c/self.lam_obs**2).to(u.erg/(u.cm**2*u.s*u.AA))
        return

    @property
    def __flam_sky(self):

        #Read the template
        sky_temp = np.loadtxt(os.environ['SPEC_PIPE_LOC']+\
                              "/Spec_pipeline/Sky_Templates/template_sky_GMOS.dat")
        lam_sky = sky_temp[:,0]*u.AA
        flam_sky_orig = sky_temp[:,1]*u.erg/(u.s*u.cm**2*u.AA)

        #Rebin the template to the object spectrum.
        self.flam_sky = rebin_spec(lam_sky, flam_sky_orig, self.lam_obs)

        return

    #There does not seem to be a way to recover the GRISM used for the
    #GMOS spectra. Parameter must be provided in object declaration.
    @property
    def __sens(self):

        #Read the sensitivity curve.
        sens_temp = np.loadtxt(os.environ['SPEC_PIPE_LOC']+\
                               "/Spec_pipeline/Sensitivity_Files/"+
                               "Sens_GMOS_"+self.grname+".txt")
        lam_sens = sens_temp[:,0]*u.AA
        sens_orig = sens_temp[:,1]*u.dimensionless_unscaled

        #Rebin the template to the object spectrum.
        self.sens = rebin_spec(lam_sens, sens_orig, self.lam_obs)

        return

    #Values taken from https://www.gemini.edu/instrumentation/gmos/components#Gratings. We will assume the resolutions at the blaze. Note that those resolving powers are for a 0.5" slit. As with other instruments, we will assume a 1" slit (so we just multiply by 2). As with DBSP, we'll assume that it is the FWHM
    @property
    def sigma_res(self):

        if self.grname=="B600":
            FWHM_res = 4160./1688. * u.AA * 2.
        elif self.grname=="R400":
            FWHM_res = 7640./1918. * u.AA * 2.
        else:
            return None

        sigma_res = FWHM_res/(2.*(2.*np.log(2.))**0.5)
        return sigma_res
