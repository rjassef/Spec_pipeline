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

   """

    def __init__(self,_name,_zspec,_fits_files,_line_center=None,
                 blue=False,red=False,show_err_plot=False):
        super(LRIS_Spec,self).__init__(_name,_zspec,_fits_files,_line_center,show_err_plot=show_err_plot)
        self.RT   = 5.0*u.m #Telescope radius.
        self.instrument = "LRIS"
        self.blue = blue
        self.red = red
        #self.edge_drop = 50.*u.AA
        self.edge_drop = 75.*u.AA
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

            #Open the fits file for the headers.
            ff = fits.open(self.data_prefix+"/"+self.fits_files[0])

            #Figure out the limits on which we can use the spectra.
            dichroic_wave = float(ff[0].header['DICHNAME'])*u.nm
            kuse = (spec_use[0].dispersion<dichroic_wave-self.edge_drop)

            #https://www2.keck.hawaii.edu/inst/lris/detectors.html
            #Use mean of amplifiers.
            self.RON = 3.82
            self.GAIN = 1.61
            self.PIXSIZE = 0.135*u.arcsec

            #Find the grism
            grism_aux = re.search("^(.*?)/.*$",spec_b[0].header['GRISNAME'])
            self.grism = "B"+grism_aux[1]
            grname = self.grism

            #Set the sky template. Remove data outisde the edges of the sky template.
            self.sky_temp_fname = os.environ['SPEC_PIPE_LOC']+\
            "/Spec_pipeline/Sky_Templates/template_sky_LRIS_b.dat"

            #Finally, assign the error name file.
            self.spec_err_name = "error."+self.fits_files[0]

        elif self.red:

            #Assign the red spectrum to be used
            spec_use = spec_r

            #Open the fits file for the headers.
            ff = fits.open(self.data_prefix+"/"+self.fits_files[1])

            #Figure out the limits on which we can use the spectra.
            dichroic_wave = float(ff[0].header['DICHNAME'])*u.nm
            kuse = (spec_use[0].dispersion>dichroic_wave+self.edge_drop)

            #https://www2.keck.hawaii.edu/inst/lris/detectors.html
            #Use mean of amplifiers.
            self.RON = 4.64
            self.GAIN = 1.197
            self.PIXSIZE = 0.135*u.arcsec

            #Find the Grating and remove data outside the edges of the sensitivity curves.
            grating_aux = re.search("^(.*?)/.*$",spec_r[0].header['GRANAME'])
            self.grating = "R"+grating_aux[1]
            grname = self.grating

            #Set the sky template. Remove data outisde the edges of the sky template.
            self.sky_temp_fname = os.environ['SPEC_PIPE_LOC']+\
            "/Spec_pipeline/Sky_Templates/template_sky_LRIS_r.dat"

            #Finally, assign the error name file.
            self.spec_err_name = "error."+self.fits_files[1]


        #Read the slit width if available in the headers. Otherwise, we'll assume the default 1.25" size from the spec class.
        if 'SLITNAME' in spec_use[0].header:
            m = re.match("long_(.*)",spec_use[0].header['SLITNAME'])
            self.slit_width = float(m.group(1)) * u.arcsec

        #If no apsize_pix read from headers, assume the slit size for the extraction aperture.
        if 'apsize_pix' in spec_use[0].header:
            self.apsize_pix = spec_use[0].['apsize_pix']
        else:
            self.apsize_pix = (self.slit_width/self.PIXSIZE).to(1.).value

        #Find the grism and remove data outside the edges of the sensitivity curves.
        sens_temp = np.loadtxt(os.environ['SPEC_PIPE_LOC']+\
                            "/Spec_pipeline/Sensitivity_Files/"+
                            "Sens_LRIS_"+grname+".txt")
        lam_sens = sens_temp[:,0]*u.AA
        kuse = (kuse) & (spec_use[0].dispersion>np.min(lam_sens)) & \
                (spec_use[0].dispersion<np.max(lam_sens))

        #Finally, figure out the sky template edges and trim the spectrum to that limit.
        sky_temp = np.loadtxt(self.sky_temp_fname)
        lam_sky = sky_temp[:,0]*u.AA
        kuse = (kuse) & (spec_use[0].dispersion>np.min(lam_sky)) & \
                (spec_use[0].dispersion<np.max(lam_sky))

        #Now, assign the wavelength and flux to the object.
        self.lam_obs = spec_use[0].dispersion[kuse]
        fnu = spec_use[0].data[kuse]*spec_use[0].unit

        #In the error name, replace fits for txt, as we will write it in ASCII.
        self.spec_err_name = re.sub(".fits",".txt",self.spec_err_name)

        #Mean bin size and exposure time. Useful for error estimation.
        self.dlam = np.mean(self.lam_obs[1:]-self.lam_obs[:-1])
        self.texp = float(ff[0].header['EXPTIME'])*u.s

        #Convert flam to fnu
        self.flam = (fnu*c/self.lam_obs**2).to(u.erg/(u.cm**2*u.s*u.AA))

        #Close the fits file.
        ff.close()

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

        #Read the sensitivity curve.
        if self.blue:
            grname = self.grism
        elif self.red:
            grname = self.grating
        else:
            return
        sens_temp = np.loadtxt(os.environ['SPEC_PIPE_LOC']+\
                               "/Spec_pipeline/Sensitivity_Files/"+
                               "Sens_LRIS_"+grname+".txt")
        lam_sens = sens_temp[:,0]*u.AA
        sens_orig = sens_temp[:,1]*u.dimensionless_unscaled

        #Rebin the template to the object spectrum.
        self.sens = rebin_spec(lam_sens, sens_orig, self.lam_obs)

        return


    #Resolutions, taken from https://www2.keck.hawaii.edu/inst/lris/dispersive_elements.html. We'll assume the value for a 1" slit. We typically used larger ones, but the seeing was probably not much larger than 1". We will also take the minimum value, as this is meant to put a minimum limit in the velocity  widths.
    @property
    def sigma_res(self):

        if self.blue:
            if self.grism=="B300":
                FWHM_res = 8.4*u.AA
            elif self.grism=="B400":
                FWHM_res = 6.5*u.AA
            elif self.grism=="B600":
                FWHM_res = 3.8*u.AA
            elif self.grism=="B1200":
                FWHM_res = 1.56*u.AA
            else:
                return None

        elif self.red:
            if self.grating=="R300":
                FWHM_res = 9.18*u.AA
            elif self.grating=="R400":
                FWHM_res = 6.9*u.AA
            elif self.grating=="R600":
                FWHM_res = 4.7*u.AA
            else:
                return None
        else:
            return None

        sigma_res = FWHM_res/(2.*(2.*np.log(2.))**0.5)
        return sigma_res
