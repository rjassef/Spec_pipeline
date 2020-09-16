#As there are too many different things to keep in minds, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

import numpy as np
import astropy.units as u
from astropy.constants import h,c
import os

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


    def load_detector_properties(self):
        #Load the detectors file.
        dets = np.genfromtxt(os.environ['SPEC_PIPE_LOC']+"/Spec_pipeline/Configurations/{0:s}_detectors.txt".format(self.instrument),dtype=str)

        if self.detector in dets[:,0]:
            det_use = dets[dets[:,0]==self.detector,:]
            self.plate_scale = float(det_use[0,1])*u.arcsec
            self.pixel_size  = float(det_use[0,2])*u.micron
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
            self.grating_dispersion = float(grt_use[0,1])*u.AA/u.mm
        else:
            print("Warning: Grating {0:s} not found.".format(self.grating))

        return
