import numpy as np
import astropy.units as u
from astropy.constants import c,h
from astropy.table import Table

from Spec_pipeline.Spec_Reader.rebin_spec import rebin_spec

class Fake_Spec(object):

    def __init__(self,z=2.5,flux_sky_scale=1.0):

        self.z    = z    #Redshift
        self.texp  = 15.*u.minute #Exposure time.
        self.flux_sky_scale = flux_sky_scale

        self.waveunit = u.AA
        self.flamunit = u.erg/u.s/u.cm**2/u.AA

        self.__flam()
        self.__flam_sky()
        self.__sens()
        self.__flam_err()

        return


    def __flam(self):
        qso = Table.read("vandenberk_composite.txt",format='ascii.cds')
        self.flam_orig_full = qso['FluxD'] * self.flamunit
        self.lam_rest_full  = qso['Wave'].to(self.waveunit)

        self.lam_obs_full = self.lam_rest_full*(1.+self.z)
        kuse = np.argwhere((self.lam_obs_full>=self.lam_spec_min) & 
                           (self.lam_obs_full<=self.lam_spec_max))
        kuse = kuse.flatten()

        self.lam_rest  = self.lam_rest_full[kuse]
        self.lam_obs   = self.lam_rest*(1.+self.z)

        self.flam      = self.flam_orig_full[kuse]*self.norm

        return

    def __flam_sky(self):
        sky_temp = np.loadtxt(self.sky_fname)
        lam_sky_temp  = sky_temp[:,0]*self.waveunit
        flam_sky_temp = sky_temp[:,1]*self.flamunit#*self.flux_sky_scale
        
        self.flam_sky = rebin_spec(lam_sky_temp, flam_sky_temp, self.lam_obs)
        return

    def __sens(self):
        sens_temp = np.loadtxt(self.sens_fname)
        lam_sens_orig  = sens_temp[:,0]*self.waveunit
        sens_orig      = sens_temp[:,1]*u.dimensionless_unscaled
        self.sens = rebin_spec(lam_sens_orig, sens_orig, self.lam_obs)
        return

    def __flam_err(self):
        
        #Transform to counts.
        self.dlam = np.mean(self.lam_obs[1:]-self.lam_obs[:-1])
        eps = self.sens * self.texp * np.pi*self.RT**2 * \
              self.dlam / (h*c/self.lam_obs)
        eps = eps.to(u.s*u.cm**2*u.AA/u.erg)
        counts     = self.flam     * eps
        counts_sky = self.flam_sky * eps
        self.flam_err = (counts + counts_sky + self.RON**2)**0.5 / eps
        self.flam_err = self.flam_err.to(self.flamunit)
        return


    #####

    @property
    def norm(self):
        fnu     = 3631.*10.**(-0.4*self.mag_norm)*u.Jy
        flam    = (fnu*c/self.lam_eff_norm**2).to(self.flamunit)
        k = np.argmin(np.abs(self.lam_obs_full-self.lam_eff_norm))
        _norm   = flam/self.flam_orig_full[k]
        return _norm


    def gen_obs_spec(self):
        self.flam.to(self.flamunit)
        self.flam_err.to(self.flamunit)
        return np.random.normal(self.flam.value,
                                self.flam_err.value)*self.flam.unit
        
