#As there are too many different things to keep in minds, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

import numpy as np
import astropy.units as u
from astropy.constants import h,c
import os
import re

#from .obtain_error_spectrum import get_error_spec
from .obtain_error_spectrum_with_extra_poly import get_error_spec

class Spec(object):

    def __init__(self,_name,_zspec,_fits_files=None,_line_center=None,show_err_plot=False,error_fit_blue_exp=True):
        self.RT   = None
        self.instrument = None
        self.name  = _name
        self.zspec = _zspec
        self.fits_files = _fits_files
        self.line_center = _line_center
        self.lam_obs  = None
        self.dlam = None
        self.texp = None
        self.RON  = None
        self.spec_err_name = None
        self.Kp = None
        self.data_prefix = "data/"
        self.save_err = True
        self.show_err_plot=show_err_plot
        self.print_err_plot=False
        self.dual_spec=False
        self.red=False
        self.blue=False

        #If no slit width is given, assume 1.25" as discussed on telecon from 05/26/2020
        self.slit_width = 1.25 * u.arcsec

    @property
    def lam_rest(self):
        try:
            return self._lam_rest
        except AttributeError:
            pass
        if self.lam_obs is not None:
            self._lam_rest = self.lam_obs/(1.+self.zspec)
        else:
            self._lam_rest = None
        return self._lam_rest

    def eps(self,sens=None,lam_obs=None,dlam=None):
        if sens is None:
            sens = self.sens
        if lam_obs is None:
            lam_obs = self.lam_obs
        if dlam is None:
            dlam = self.dlam
        self._eps = sens * (np.pi*self.RT**2) * dlam * self.texp/\
                    (h*c/lam_obs)
        self._eps = self._eps.to(u.s*u.cm**2*u.AA/u.erg)
        return self._eps

    @property
    def flam_err(self):
        try:
            return self._flam_err
        except AttributeError:
            pass
        #Try to open the error spectrum. If not found, generate it.
        try:
            cat = open(self.data_prefix+"/"+self.spec_err_name,"r")
            self._flam_err = np.loadtxt(cat,usecols=[1])
            self._flam_err = self._flam_err * u.erg/(u.s*u.cm**2*u.AA)
            cat.close()
        except IOError:
            self._flam_err, self.Kp = get_error_spec(self, wd=15, show_plot=self.show_err_plot)
            if self.save_err:
                np.savetxt(self.data_prefix+"/"+self.spec_err_name,
                           np.array([self.lam_obs,
                                     self._flam_err]).T)
        return self._flam_err

    #####

    def load_detector_properties(self):
        #Load the detectors file.
        dets = np.genfromtxt(os.environ['SPEC_PIPE_LOC']+"/Spec_pipeline/Configurations/{0:s}_detectors.txt".format(self.instrument),dtype=str)

        if self.detector in dets[:,0]:
            det_use = dets[dets[:,0]==self.detector,:]
            kws = ["plate_scale", "pixel_size", "RON", "GAIN"]
            units_kws = ["arcsec", "micron", "", ""]
            for k, kw in enumerate(kws):
                if (not hasattr(self,kw)) or (getattr(self,kw) is None):
                    try:
                        setattr(self, kw, float(det_use[0,k+1])*u.Unit(units_kws[k]))
                    except ValueError:
                        setattr(self, kw, None)
            #self.plate_scale = float(det_use[0,1])*u.arcsec
            #self.pixel_size  = float(det_use[0,2])*u.micron
        else:
            print("Warning: Detector {0:s} not found.".format(self.detector))

        return

    def load_grating_properties(self):

        #Load the detectors file.
        if self.dual_spec:
            grts = np.genfromtxt(os.environ['SPEC_PIPE_LOC']+"/Spec_pipeline/Configurations/{0:s}_gratings_{1:s}.txt".format(self.instrument,self.channel),dtype=str)
        else:
            grts = np.genfromtxt(os.environ['SPEC_PIPE_LOC']+"/Spec_pipeline/Configurations/{0:s}_gratings.txt".format(self.instrument),dtype=str)

        if self.grating in grts[:,0]:
            grt_use = grts[grts[:,0]==self.grating,:]
            kws = ["grating_dispersion", "FWHM_res"]
            units_kws = ["Angstrom/mm", "Angstrom"]
            for k,kw in enumerate(kws):
                if (not hasattr(self,kw)) or (getattr(self,kw) is None):
                    try:
                        setattr(self, kw, float(grt_use[0,k+1])*u.Unit(units_kws[k]))
                    except ValueError:
                        setattr(self, kw, None)
            #self.grating_dispersion = float(grt_use[0,1])*u.AA/u.mm
        else:
            print("Warning: Grating {0:s} not found.".format(self.grating))

        return

    def load_keyword_headers(self, spec, keywords):

        for key in keywords.keys():
            if (not hasattr(self,key)) or (getattr(self,key) is None):
                if keywords[key] in spec[0].header.keys():
                    setattr(self,key, spec[0].header[keywords[key]])
                else:
                    print("Warning: {0:s} keyword not found.".format(keywords[key]))
                    setattr(self,key, None)

    def run_setup(self, spec_use):

        #Cut the edges close to the dichroic.
        if (hasattr(self,'dichroic')) and (self.dichroic is not None):
            if self.blue:
                kuse = (spec_use[0].dispersion<self.dichroic_wave-self.edge_drop)
            else:
                kuse = (spec_use[0].dispersion>self.dichroic_wave+self.edge_drop)
        else:
            kuse = (spec_use[0].dispersion>0*u.AA)

        #Load the detector properties.
        if (hasattr(self,'detector')) and (self.detector is not None):
            self.load_detector_properties()

        #Load the grating properties.
        if (hasattr(self,'grating')) and (self.grating is not None):
            self.load_grating_properties()

        #If no apsize_pix read from headers, assume the slit size for the extraction aperture.
        if self.apsize_pix is None:
            try:
                self.apsize_pix = (self.slit_width/self.plate_scale).to(1.).value
            except TypeError:
                pass

        #Set the sensitivity template.
        if self.local_sens_files is None:
            self.sens_temp_fname = os.environ['SPEC_PIPE_LOC'] + "/Spec_pipeline/Sensitivity_Files/" + "Sens_{0:s}_{1:s}_{2:s}_{3:s}_{4:s}.txt".format(self.instrument, self.detector, self.grating, self.dichroic, self.channel)
        else:
            if self.blue:
                self.sens_temp_fname = self.local_sens_files[0]
            else:
                self.sens_temp_fname = self.local_sens_files[1]

        #Set the sky template to use.
        if self.local_sky_files is None:
            self.sky_temp_fname = os.environ['SPEC_PIPE_LOC'] + "/Spec_pipeline/Sky_Templates/" + "template_sky_{0:s}_{1:s}_{2:.2f}arcsec_{3:s}.txt".format( self.instrument, self.grating, self.slit_width.to(u.arcsec).value, self.channel)
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
        self.texp = np.float(self.texp) * u.s

        #Convert fnu to flambda.
        self.flam = (fnu*c/self.lam_obs**2).to(u.erg/(u.cm**2*u.s*u.AA))

        return
