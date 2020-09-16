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

#As there are too many different things to keep in minds, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

class DBSP_Spec(Spec):

    """
Module that reads a DBSP spectrum and returns a spec object.

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

   show_err_plot (boolean)    : Optional. True if error-fit plot is to be
                                displayed.

   local_sky_files (list)     : Optional. List of sky files if the default ones
                                are not to be used.

   local_sens_files (list)    : Optional. List of sensitivity files if the
                                default ones are not to be used.

   dichroic (string)          : Optional.

   grating (string)           : Optional.

   grating_dispersion (float) : Optional. In astropy units of AA/mm

   detector (string)          : Optional.

   plate_scale (float)        : Optional. In astropy units of arcsec.

   pixel_size (float)         : Optional. In astropy units of microns.

   slit_width (float)         : Optional. In astropy units of arcsec.

   """

    def __init__(self,_name,_zspec,_fits_files,_line_center=None,
                 blue=False,red=False,show_err_plot=False,local_sky_files=None,local_sens_files=None, dichroic=None, grating=None, grating_dispersion=None, detector=None, plate_scale=None, pixel_size=None, slit_width=None):
        super(DBSP_Spec,self).__init__(_name,_zspec,_fits_files,_line_center,show_err_plot=show_err_plot)
        self.RT   = 2.5*u.m #Telescope radius.
        self.instrument = "DBSP"
        self.dual_spec = True
        self.blue = blue
        self.red  = red
        self.edge_drop = 75.*u.AA
        self.local_sky_files = local_sky_files
        self.local_sens_files = local_sens_files

        self.dichroic = dichroic
        self.grating = grating
        self.grating_dispersion = grating_dispersion
        self.detector = detector
        self.plate_scale = plate_scale
        self.pixel_size = pixel_size
        self.slit_width = slit_width

        self._sigma_res = None

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
        if self.dichroic is None:
            if 'DICHROIC' in spec_use[0].header.keys():
                self.dichroic = spec_use[0].header['DICHROIC']
                self.dichroic = re.sub("-","",self.dichroic)
            else:
                print("Warning: DICHROIC keyword not found. No dichroic set")

        if self.dichroic is not None:
            dichroic_wave = float(self.dichroic[1:])*100.*u.AA
            if self.blue:
                kuse = (spec_use[0].dispersion<dichroic_wave-self.edge_drop)
            else:
                kuse = (spec_use[0].dispersion>dichroic_wave+self.edge_drop)
        else:
            kuse = (spec_use[0].dispersion>0*u.AA)


        #Detector
        if self.detector is None:
            if 'DETNAM' in spec_use[0].header.keys():
                self.detector = spec_use[0].header['DETNAM']
            else:
                print("Warning: DETNAM keyword not found. No detector set")

        if self.detector is not None:
            self.load_detector_properties()

        #Grating - there is a grating called sometimes 316/7150 and sometimes 316/7500, but it is the same grating with the same blaze as confirmed by the Palomar Observatory staff. We will alwayd use 316/7500.
        if self.grating is None:
            if 'GRATING' in spec_use[0].header.keys():
                self.grating  = spec_use[0].header['GRATING']
                self.grating  = re.sub("7150","7500",self.grating)
                self.grating  = re.sub("/","-",self.grating)
                if self.grating == "3167500":
                    self.grating = "316-7500"
            else:
                print("Warning: GRATING keyword not found. No grating set")

        if self.grating is not None:
            self.load_grating_properties()

        #Slit width
        if self.slit_width is None:
            if 'APERTURE' in spec_use[0].header.keys():
                self.slit_width = float(spec_use[0].header['APERTURE']) * u.arcsec
            else:
                print("Warning: APERTURE keyword not found. No slit width set")


        #If no apsize_pix read from headers, assume the slit size for the extraction aperture.
        if 'apsize_pix' in spec_use[0].header:
            self.apsize_pix = spec_use[0].header['apsize_pix']
        else:
            self.apsize_pix = (self.slit_width/self.plate_scale).to(1.).value

        #Set the sensitivity template.
        if self.local_sens_files is None:
            self.sens_temp_fname = os.environ['SPEC_PIPE_LOC'] + "/Spec_pipeline/Sensitivity_Files/" + "Sens_DBSP_{0:s}_{1:s}_{2:s}_{3:s}.txt".format(self.detector, self.grating, self.dichroic, self.channel)
        else:
            if self.blue:
                self.sens_temp_fname = self.local_sens_files[0]
            else:
                self.sens_temp_fname = self.local_sens_files[1]

        #Set the sky template to use.
        if self.local_sky_files is None:
            self.sky_temp_fname = os.environ['SPEC_PIPE_LOC'] + "/Spec_pipeline/Sky_Templates/" + "template_sky_DBSP_{0:s}_{1:.2f}arcsec_{2:s}.txt".format(self.grating,self.slit_width.to(u.arcsec).value,self.channel)
        else:
            if self.blue:
                self.sky_temp_fname = self.local_sky_files[0]
            else:
                self.sky_temp_fname = self.local_sky_files[1]

        #Finally, figure out the sky template edges and trim the spectrum to that limit.
        try:
            sky_temp = np.loadtxt(self.sky_temp_fname)
        except IOError:
            print("Could not open sky file ",self.sky_temp_fname)
            return
        lam_sky = sky_temp[:,0]*u.AA
        kuse_sky = (spec_use[0].dispersion>np.min(lam_sky)) & \
                (spec_use[0].dispersion<np.max(lam_sky))

        #Display a warning if we are missing any range because of the sensitivity curve or the sky template.
        lam = spec_use[0].dispersion[kuse]
        if np.min(lam)<np.min(lam_sky) or np.max(lam)>np.max(lam_sky):
            print("Wavelength range for object {0:s} limited because of sky template".format(self.name))
            print("Spec-range: {0:.1f} - {1:.2f}".format(np.min(lam),np.max(lam)))
            print("Sky-range: {0:.1f} - {1:.2f}".format(np.min(lam_sky),np.max(lam_sky)))
        kuse = (kuse) & (kuse_sky)

        #Now, assign the wavelength and flux to the object.
        self.lam_obs = spec_use[0].dispersion[kuse]
        fnu = spec_use[0].data[kuse]*spec_use[0].unit

        #Change the .fits for .txt in the error file name as it will be saved in ASCII
        self.spec_err_name = re.sub(".fits",".txt",self.spec_err_name)

        #Mean bin size, exposure time, RON and GAIN. Useful for error estimation.
        self.dlam = np.mean(self.lam_obs[1:]-self.lam_obs[:-1])
        self.texp = float(spec_use[0].header['EXPTIME'])*u.s
        self.RON  = float(spec_use[0].header['RON'])
        self.GAIN = float(spec_use[0].header['GAIN'])

        #Convert fnu to flambda.
        self.flam = (fnu*c/self.lam_obs**2).to(u.erg/(u.cm**2*u.s*u.AA))

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

        sens_temp = np.loadtxt(self.sens_temp_fname)
        lam_sens = sens_temp[:,0]*u.AA
        sens_orig = sens_temp[:,1]*u.dimensionless_unscaled

        #Rebin the template to the object spectrum.
        #self.sens = rebin_spec(lam_sens, sens_orig, self.lam_obs)

        #Interpolate the sensitivity template to the object spectrum. Extrapolate if needed, which is OK as it is a smooth function of wavelegnth for the most part.
        if np.min(self.lam_obs)<np.min(lam_sens) or np.max(self.lam_obs)>np.max(lam_sens):
            print("Warning: Extrapolating sensitivity curve to match spectral range")
            print("Spec-range: {0:.1f} - {1:.2f}".format( np.min(self.lam_obs),np.max(self.lam_obs)))
            print("Sens-range: {0:.1f} - {1:.2f}".format(np.min(lam_sens),np.max(lam_sens)))
        f = interp1d(lam_sens, sens_orig, kind='linear', fill_value='extrapolate')
        self.sens = f(self.lam_obs)

        return

    #We use the numbers here: https://www.astro.caltech.edu/palomar/observer/200inchResources/dbspoverview.html#grating . We'll assume that the dispersion numbers quoted are the FWHM (like for LRIS) and we'll assume a 1" slit, the same as we did for LRIS, since really what matters here is the minimum between the seeing and the slit size, and we want to be conservative.
    @property
    def sigma_res(self):

        if self._sigma_res is not None:
            return self._sigma_res

        slit_size = 1.0*u.arcsec

        if self.grating_dispersion is not None:
            res = self.grating_dispersion
        else:
             print("Grating dispersion not set.")
             print("Using minimum of 1 AA/mm ")
             res = 1.0*u.AA/u.mm

        FWHM_res = (slit_size/self.plate_scale)*self.pixel_size * res
        self._sigma_res = (FWHM_res/(2.*(2.*np.log(2.))**0.5)).to(u.AA)
        return self._sigma_res
