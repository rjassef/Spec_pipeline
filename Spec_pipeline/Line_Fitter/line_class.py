# Implementation of the line fitting routines. This is the main class,
# we will have a separate class with different defaults for different
# emission lines.

import numpy as np
import astropy.units as u
from astropy.constants import c
from scipy.special import betainc

from . import MC_errors as MC
from . import fit_chi2 as fit
from . import plot_fit

class Line_fit(object):
    
    def __init__(self,_line_name,_line_center,_line_velocity_region,
                 _continuum_regions):
        self.line_name   = _line_name
        self.line_center = _line_center
        self.line_velocity_region = _line_velocity_region
        self.continuum_regions = _continuum_regions

        #Values that will come from the fitting.
        self.lam_cen_fit       = None
        self.flam_line_cen_fit = None
        self.sigma_v_fit       = None
        self.a                 = None
        self.b                 = None
        self.flam_mod          = None

        #Values that will come from the MC
        self.lam_cen_low = None
        self.lam_cen_hig = None
        self.flam_line_cen_low = None
        self.flam_line_cen_hig = None
        self.sigma_v_low = None
        self.sigma_v_hig = None

        #Values from Ftest
        self.chi2 = None
        self.chi2_no_line = None
        self.F = None
        self.p = None

        #Set the default constraints for the fitting.
        self.sigma_v_min =  100.*u.km/u.s
        self.sigma_v_max = 5000.*u.km/u.s
        self.delta_lam_cen_max = 2.*u.AA
        #Minimum peak height fraction compared to initial guess.
        self.frac_flam_line_cen_min = 1e-1

    @property
    def FWHM_v_low(self):
        if self.sigma_v_low is None:
            print("First run the MC")
            return
        return self.sigma_v_low*2.*(2.*np.log(2.))**0.5

    @property
    def FWHM_v_hig(self):
        if self.sigma_v_hig is None:
            print("First run the MC")
            return
        return self.sigma_v_hig*2.*(2.*np.log(2.))**0.5

    @property
    def FWHM_v(self):
        if self.sigma_v_fit is None:
            print("First fit the line")
            return
        return self.sigma_v_fit*2.*(2.*np.log(2.))**0.5

    def run_fit(self, spec,lam_cen_0=None,
                sigma_v_0=3000.*u.km/u.s):
        if lam_cen_0 is None:
            lam_cen_0 = self.line_center
        if spec.lam_rest is None or \
           lam_cen_0<np.min(spec.lam_rest) or \
           lam_cen_0>np.max(spec.lam_rest):
            print("Line not within spectrum")
            return
        self.lam_cen_fit, self.flam_line_cen_fit, self.sigma_v_fit, \
            self.a, self.b, self.flam_mod = \
                                            fit.fit(spec,self, 
                                                    lam_cen_0, sigma_v_0)
        return

    def run_MC(self,spec,nrep,Ncpu=None,save_chain=None):
        if self.sigma_v_fit is None:
            print("First fit the line")
            return
        [self.lam_cen_low, self.lam_cen_hig,\
         self.flam_line_cen_low, self.flam_line_cen_hig,\
         self.sigma_v_low, self.sigma_v_hig] = MC.MC_errors(
             nrep, spec, self, 
             Ncpu=Ncpu, save_chain=save_chain)
        return

    def plot(self,spec,plot_fname=None):
        plot_fit.plot_fit(spec.lam_rest, spec.flam, self.flam_mod, 
                          self.lam_cen_fit,
                          self.continuum_regions, 
                          self.line_velocity_region,
                          spec.name, self.FWHM_v,plot_fname=plot_fname)

    def run_Ftest(self,spec):
        if self.sigma_v_fit is None:
            print("First fit the line")
            return
        self.chi2 = fit.chi2_fit([self.lam_cen_fit.value,
                                  self.flam_line_cen_fit.value,
                                  self.sigma_v_fit.value],
                                 spec, self, self.a, self.b,
                                 self.lam_cen_fit, 
                                 self.flam_line_cen_fit,
                                 check_constraints=False)
        self.chi2_no_line = fit.chi2_fit([self.lam_cen_fit.value,
                                          0.,
                                          self.sigma_v_fit.value],
                                         spec, self, self.a, self.b,
                                         self.lam_cen_fit, 
                                         self.flam_line_cen_fit,
                                         check_constraints=False)
            
        #Get the degrees of freedom.
        dv = (c*(spec.lam_rest/self.lam_cen_fit-1.)).to(u.km/u.s)
        dvabs = np.abs(dv)
        n_datapoints = len(spec.flam[dvabs<self.line_velocity_region])
        nu = n_datapoints - 2 - 3
        nu_no_line = n_datapoints
        chi2_nu = self.chi2/float(nu)
        chi2_no_line_nu = self.chi2_no_line/float(nu_no_line)

        self.F = ((self.chi2_no_line-self.chi2)/float(nu_no_line-nu)) / \
                 (self.chi2/float(nu))
        self.F = self.F.value
        if self.F<0.:
            self.p = 1.0
        else:
            nnu1 = float(nu_no_line-nu)
            nnu2 = float(nu)
            w = nnu1*self.F/(nnu1*self.F+nnu2)
            self.p = 1.-betainc(nnu1/2.,nnu2/2.,w)
        return
