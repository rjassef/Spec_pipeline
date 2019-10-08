import numpy as np
import astropy.units as u
from astropy.constants import c,h
from astropy.table import Table

from Spec_pipeline.Spec_Reader.rebin_spec import rebin_spec

class Fake_LRIS_Spec(object):

    def __init__(self,_side,z=2.5,Vmag=24.0):

        self.z    = z    #Redshift
        self.Vmag = Vmag #Normalization magnitude
        self.side = _side
        if self.side!='blue' and self.side!='red':
            print("Please declare red or blue side")
            return
        
        self.RT    = 5.*u.m       #Telescope radius
        self.texp  = 15.*u.minute #Exposure time.
        
        if self.side=='blue':
            self.RON = 3.82
        else:
            self.RON = 4.64

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

        if self.side=='blue':
            lam_min = 3200.*u.AA
            lam_max = 5500.*u.AA
        else:
            lam_min =  5500.*u.AA
            lam_max = 10000.*u.AA

        self.lam_obs_full = self.lam_rest_full*(1.+self.z)
        kuse = np.argwhere((self.lam_obs_full>=lam_min) & 
                           (self.lam_obs_full<=lam_max))
        kuse = kuse.flatten()

        self.lam_rest  = self.lam_rest_full[kuse]
        self.lam_obs   = self.lam_rest*(1.+self.z)

        self.flam      = self.flam_orig_full[kuse]*self.norm

        return

    def __flam_sky(self):
        
        if self.side=='blue':
            fname = "../Spec_pipeline/Sky_Templates/template_sky_LRIS_b.dat"
        else:
            fname = "../Spec_pipeline/Sky_Templates/template_sky_LRIS_r.dat"

        sky_temp = np.loadtxt(fname)
        lam_sky_temp  = sky_temp[:,0]*self.waveunit
        flam_sky_temp = sky_temp[:,1]*self.flamunit*1e-3
        self.flam_sky = rebin_spec(lam_sky_temp, flam_sky_temp, self.lam_obs)
        return

    def __sens(self):
        
        if self.side=='blue':
            fname = "../Spec_pipeline/Sensitivity_Files/Sens_LRIS_B600.txt"
        else:
            fname = "../Spec_pipeline/Sensitivity_Files/Sens_LRIS_R400.txt"
        sens_temp = np.loadtxt(fname)
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
        lam_eff = 5500.*u.AA
        fnu_V   = 3631.*10.**(-0.4*self.Vmag)*u.Jy
        flam_V  = (fnu_V*c/lam_eff**2).to(self.flamunit)
        k = np.argmin(np.abs(self.lam_obs_full-lam_eff))
        _norm   = flam_V/self.flam_orig_full[k]
        return _norm


    def gen_obs_spec(self):
        self.flam.to(self.flamunit)
        self.flam_err.to(self.flamunit)
        return np.random.normal(self.flam.value,
                                self.flam_err.value)*self.flam.unit
        
