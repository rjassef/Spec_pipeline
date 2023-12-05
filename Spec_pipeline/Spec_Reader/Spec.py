#As there are too many different things to keep in minds, we'll load
#the spectra as objects, so we can load the appropiate sensitivity
#curves and sky spectra without having to think too much about it
#during the code execution.

import numpy as np
import astropy.units as u
from astropy.constants import h,c
import os
import re
from scipy.interpolate import interp1d

#from .obtain_error_spectrum import get_error_spec
from .obtain_error_spectrum_with_extra_poly import get_error_spec
from .rebin_spec import rebin_spec

class Spec(object):

    """
    Superclass for spectrum objects. Subclasses for different instruments have been created, and it is recommended to use those as a template for setting up a new instrument rather than trying to create directly as spec object.

    Parameters
    ----------
    name : string
        Object ID

    zspec : float
        Spectroscopic redshift.

    fits_file : string
        Spectrum file name.

    show_err_plot : boolean, optional
        True if error-fit plot is to be displayed.

    error_fit_blue_exp : boolean, optional
        If True, forces the error estimations to use an exponential to model extra uncertainty in the blue edge. Empirically seen to be needed by most spectra.

    local_sky_file : string, optional
        Sky file if the default ones are not to be used.

    local_sens_file : string, optional
        Sensitivity file if the default ones are not to be used.

    inst_conf : dictionary, optional
        Configurations dictionary.

    header_kws : dictionary, optional
        Dictionary with header keywords to use. Have precedence over default header keywords.

    RT : float with astropy.units of meters, optional
        Aperture radius of the telescope. Default is 1.0 m.

    instrument : string, optional
        Instrument's name. Should be consistent with the name used in the default configuration files.

    """

    def __init__(self, name, zspec, fits_file=None, show_err_plot=False, error_fit_blue_exp=True, local_sky_file=None, local_sens_file=None, inst_conf=None, header_kws=None, RT=1.0*u.m, instrument=None):

        self.name  = name
        self.zspec = zspec
        self.fits_file = fits_file
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

        self.RT = RT
        self.instrument = instrument
        if self.instrument is None:
            print("Must provide instrument's name.")
            return 1

        if not hasattr(self,"dual_spec"):
            self.dual_spec = False
            self.red=False
            self.blue=False

        if self.dual_spec and (not self.blue and not self.red):
            print("Must indicate whether this a blue or red arm spectrum.")
            return 1

        self.local_sky_file = local_sky_file
        self.local_sens_file = local_sens_file

        self.dichroic = None
        self.grating = None
        self.grating_dispersion = None
        self.detector = None
        self.plate_scale = None
        self.pixel_size = None
        self.slit_width = None

        self._sigma_res = None

        #If no slit width is given, assume 1.25" as discussed on telecon from 05/26/2020
        #self.slit_width = 1.25 * u.arcsec

        #Load the user provided configurations.
        if inst_conf is not None:
            for kw in inst_conf.keys():
                kwuse = kw
                if kwuse[:4]=='blue' and self.blue:
                    kwuse = kw[5:]
                elif kwuse[:3]=='red' and self.red:
                    kwuse = kw[4:]
                setattr(self,kwuse,inst_conf[kw])



        #Set the default header keywords and overwrite them with, or add to them, the ones set by the user.
        self.keywords_to_load = dict()
        try:
            cat = open(os.environ['SPEC_PIPE_LOC']+"/Spec_pipeline/Configurations/{0:s}_Header_Keywords.txt".format(self.instrument))
            for line in cat:
                if line[0]=='#':
                    continue
                x = line.split()
                if len(x)<2:
                    continue
                if len(x)>2:
                    if x[2]=='blue' and not self.blue:
                        continue
                    elif x[2]=='red' and not self.red:
                        continue
                self.keywords_to_load[x[0]] = x[1]
        except IOError:
            pass

        #Add spec aperture keyword.
        self.keywords_to_load["apsize_pix"] = "apsize_pix"

        #Add the user provided ones.
        if header_kws is not None:
            for kw in header_kws.keys():
                self.keywords_to_load[kw] = header_kws[kw]

        return 0

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
        if (hasattr(self,'dichroic_wave')) and (self.dichroic_wave is not None):
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
        if self.local_sens_file is None:
            #self.sens_temp_fname = os.environ['SPEC_PIPE_LOC'] + "/Spec_pipeline/Sensitivity_Files/" + "Sens_{0:s}_{1:s}_{2:s}_{3:s}_{4:s}.txt".format(self.instrument, self.detector, self.grating, self.dichroic, self.channel)
            self.sens_temp_fname = os.environ['SPEC_PIPE_LOC'] + "/Spec_pipeline/Sensitivity_Files/" + "Sens_{0:s}_{1:s}_{2:s}".format(self.instrument, self.detector, self.grating)
            if self.dual_spec:
                self.sens_temp_fname += "_{0:s}_{1:s}".format(self.dichroic, self.channel)
            self.sens_temp_fname += ".txt"
        else:
            self.sens_temp_fname = self.local_sens_file
            # if self.blue:
            #     self.sens_temp_fname = self.local_sens_files[0]
            # else:
            #     self.sens_temp_fname = self.local_sens_files[1]

        #Set the sky template to use.
        if self.local_sky_file is None:
            #self.sky_temp_fname = os.environ['SPEC_PIPE_LOC'] + "/Spec_pipeline/Sky_Templates/" + "template_sky_{0:s}_{1:s}_{2:.2f}arcsec_{3:s}.txt".format( self.instrument, self.grating, self.slit_width.to(u.arcsec).value, self.channel)
            self.sky_temp_fname = os.environ['SPEC_PIPE_LOC'] + "/Spec_pipeline/Sky_Templates/" + "template_sky_{0:s}_{1:s}_{2:.2f}arcsec".format( self.instrument, self.grating, self.slit_width.to(u.arcsec).value)
            if self.dual_spec:
                self.sky_temp_fname += "_{0:s}".format(self.channel)
            self.sky_temp_fname += ".txt"
        else:
            self.sky_temp_fname = self.local_sky_file
            # if self.blue:
            #     self.sky_temp_fname = self.local_sky_files[0]
            # else:
            #     self.sky_temp_fname = self.local_sky_files[1]

        #Finally, figure out the sky template edges and trim the spectrum to that limit.
        try:
            sky_temp = np.loadtxt(self.sky_temp_fname)
            lam_sky = sky_temp[:,0]*u.AA
            lam_sky_min = np.min(lam_sky)
            lam_sky_max = np.max(lam_sky)
        except (IOError, OSError):
            print("Could not open sky file ",self.sky_temp_fname)
            lam_sky_min = 0.*u.AA
            lam_sky_max = 1.e6*u.AA

        kuse_sky = (spec_use[0].dispersion>lam_sky_min) & \
                (spec_use[0].dispersion<lam_sky_max)

        #Display a warning if we are missing any range because of the sensitivity curve or the sky template.
        lam = spec_use[0].dispersion[kuse]
        if np.min(lam)<lam_sky_min or np.max(lam)>lam_sky_max:
            print("Wavelength range for object {0:s} limited because of sky template".format(self.name))
            print("Spec-range: {0:.1f} - {1:.2f}".format(np.min(lam),np.max(lam)))
            print("Sky-range: {0:.1f} - {1:.2f}".format(lam_sky_min,lam_sky_max))
        kuse = (kuse) & (kuse_sky)

        #Now, assign the wavelength and flux to the object.
        self.lam_obs = spec_use[0].dispersion[kuse]
        fnu = spec_use[0].data[kuse]*spec_use[0].unit

        #Change the .fits for .txt in the error file name as it will be saved in ASCII
        self.spec_err_name = re.sub(".fits",".txt",self.spec_err_name)

        #Mean bin size, exposure time, RON and GAIN. Useful for error estimation.
        self.dlam = np.mean(self.lam_obs[1:]-self.lam_obs[:-1])
        try:
            self.texp.unit
        except AttributeError:
            self.texp = float(self.texp) * u.s

        #Convert fnu to flambda if in fnu units. If not, raise a warning.
        try:
            (fnu/(u.erg/u.s/u.cm**2/u.Hz)).to(u.dimensionless_unscaled)
            self.flam = (fnu*c/self.lam_obs**2).to(u.erg/(u.cm**2*u.s*u.AA))
        except u.UnitConversionError:
            print("Could not convert fnu to flam due to units issue. Leaving spectrum with original units.")
            self.flam = fnu

        return

    @property
    def __sens(self):

        try:
            sens_temp = np.loadtxt(self.sens_temp_fname)
        except (IOError,OSError):
            print("Could not open file {0:s}".format(self.sens_temp_fname))
            return
        lam_sens = sens_temp[:,0]*u.AA
        sens_orig = sens_temp[:,1]*u.dimensionless_unscaled

        #Rebin the template to the object spectrum.
        #self.sens = rebin_spec(lam_sens, sens_orig, self.lam_obs)

        #Interpolate the sensitivity template to the object spectrum. Extrapolate if needed, which is OK as it is a smooth function of wavelegnth for the most part.
        if np.min(self.lam_obs)<np.min(lam_sens) or np.max(self.lam_obs)>np.max(lam_sens):
            print("Warning: Extrapolating sensitivity curve to match spectral range {0:s}".format(self.name))
            print("Spec-range: {0:.1f} - {1:.2f}".format( np.min(self.lam_obs),np.max(self.lam_obs)))
            print("Sens-range: {0:.1f} - {1:.2f}".format(np.min(lam_sens),np.max(lam_sens)))
        f = interp1d(lam_sens, sens_orig, kind='linear', fill_value='extrapolate')
        self.sens = f(self.lam_obs)

        return

    @property
    def __flam_sky(self):

        #Read the template
        try:
            sky_temp = np.loadtxt(self.sky_temp_fname)
        except (IOError, OSError):
            print("Could not open file {0:s}".format(self.sky_temp_fname))
            return

        lam_sky = sky_temp[:,0]*u.AA
        flam_sky_orig = sky_temp[:,1]*u.erg/(u.s*u.cm**2*u.AA)

        #Rebin the template to the object spectrum.
        self.flam_sky = rebin_spec(lam_sky, flam_sky_orig, self.lam_obs)

        return

    @property
    def sigma_res(self):

        if self._sigma_res is not None:
            return self._sigma_res

        if self.FWHM_res is None:
            slit_size = 1.0*u.arcsec
            if self.grating_dispersion is not None:
                res = self.grating_dispersion
            else:
                print("Grating dispersion not set.")
                print("Using minimum of 1 AA/mm ")
                res = 1.0*u.AA/u.mm
            self.FWHM_res = (slit_size/self.plate_scale)*self.pixel_size * res

        self._sigma_res = (self.FWHM_res/(2.*(2.*np.log(2.))**0.5)).to(u.AA)
        return self._sigma_res
