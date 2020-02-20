#!/usr/bin/env python

import numpy as np
#from specutils.io.read_fits import read_fits_spectrum1d
from astropy.io import fits
import astropy.units as u
from astropy.constants import h,c
import re
import os

from .spectrum1d import read_fits_spectrum1d
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
                 blue=False,red=False):
        super(LRIS_Spec,self).__init__(_name,_zspec,_fits_files,_line_center)
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

            #Figure out the dichroic and cut the wavelength range.
            ff = fits.open(self.data_prefix+"/"+self.fits_files[0])
            dichroic_wave = float(ff[0].header['DICHNAME'])*u.nm
            kuse = (spec_b[0].dispersion<dichroic_wave-self.edge_drop)

            #https://www2.keck.hawaii.edu/inst/lris/detectors.html
            #Use mean of amplifiers.
            self.RON = 3.82

            #Find the grism and remove data outside the edges of the sensitivity curves.
            grism_aux = re.search("^(.*?)/.*$",spec_b[0].header['GRISNAME'])
            self.grism = "B"+grism_aux[1]
            sens_temp = np.loadtxt(os.environ['SPEC_PIPE_LOC']+\
                               "/Spec_pipeline/Sensitivity_Files/"+
                               "Sens_LRIS_"+self.grism+".txt")
            lam_sens = sens_temp[:,0]*u.AA
            kuse = (kuse) & (spec_b[0].dispersion>np.min(lam_sens)) & \
                (spec_b[0].dispersion<np.max(lam_sens))

            #Set the sky template. Remove data outisde the edges of the sky template.
            self.sky_temp_fname = os.environ['SPEC_PIPE_LOC']+\
            "/Spec_pipeline/Sky_Templates/template_sky_LRIS_b.dat"
            sky_temp = np.loadtxt(self.sky_temp_fname)
            lam_sky = sky_temp[:,0]*u.AA
            kuse = (kuse) & (spec_b[0].dispersion>np.min(lam_sky)) & \
                (spec_b[0].dispersion<np.max(lam_sky))

            #Finally, assign the wavelength and flux to the object.
            self.lam_obs = spec_b[0].dispersion[kuse]
            fnu = spec_b[0].data[kuse]*spec_b[0].unit
            self.spec_err_name = "error."+self.fits_files[0]

        elif self.red:
            ff = fits.open(self.data_prefix+"/"+self.fits_files[1])
            dichroic_wave = float(ff[0].header['DICHNAME'])*u.nm
            kuse = (spec_r[0].dispersion>dichroic_wave+self.edge_drop)

            #https://www2.keck.hawaii.edu/inst/lris/detectors.html
            #Use mean of amplifiers.
            self.RON = 4.64

            #Find the Grating and remove data outside the edges of the sensitivity curves.
            grating_aux = re.search("^(.*?)/.*$",spec_r[0].header['GRANAME'])
            self.grating = "R"+grating_aux[1]
            sens_temp = np.loadtxt(os.environ['SPEC_PIPE_LOC']+\
                               "/Spec_pipeline/Sensitivity_Files/"+
                               "Sens_LRIS_"+self.grating+".txt")
            lam_sens = sens_temp[:,0]*u.AA
            kuse = (kuse) & (spec_r[0].dispersion>np.min(lam_sens)) & \
                (spec_r[0].dispersion<np.max(lam_sens))

            #Set the sky template. Remove data outisde the edges of the sky template.
            self.sky_temp_fname = os.environ['SPEC_PIPE_LOC']+\
            "/Spec_pipeline/Sky_Templates/template_sky_LRIS_r.dat"
            sky_temp = np.loadtxt(self.sky_temp_fname)
            lam_sky = sky_temp[:,0]*u.AA
            kuse = (kuse) & (spec_r[0].dispersion>np.min(lam_sky)) & \
                (spec_r[0].dispersion<np.max(lam_sky))

            #Finally, assign the wavelength and flux to the object.
            self.lam_obs = spec_r[0].dispersion[kuse]
            fnu = spec_r[0].data[kuse]*spec_r[0].unit
            self.spec_err_name = "error."+self.fits_files[1]


        lam_sky_min = np.min(sky_temp[:,0])*u.AA
        lam_sky_max = np.max(sky_temp[:,0])*u.AA

        self.spec_err_name = re.sub(".fits",".txt",self.spec_err_name)

        self.dlam = np.mean(self.lam_obs[1:]-self.lam_obs[:-1])#Mean lambda bin.
        self.texp = float(ff[0].header['EXPTIME'])*u.s

        self.flam = (fnu*c/self.lam_obs**2).to(u.erg/(u.cm**2*u.s*u.AA))
        ff.close()
        return

    @property
    def __flam_sky(self):

        '''#Figure out the spectrograph arm.
        if self.blue:
            sky_temp_fname = os.environ['SPEC_PIPE_LOC']+\
                "/Spec_pipeline/Sky_Templates/template_sky_LRIS_b.dat"
        elif self.red:
            sky_temp_fname = os.environ['SPEC_PIPE_LOC']+\
                "/Spec_pipeline/Sky_Templates/template_sky_LRIS_r.dat"
        else:
            #print("Cannot find spectrograph arm flag")
            return'''

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
