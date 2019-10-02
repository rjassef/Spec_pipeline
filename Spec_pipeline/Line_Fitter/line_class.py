# Implementation of the line fitting routines. This is the main class,
# we will have a separate class with different defaults for different
# emission lines.

import numpy as np
import astropy.units as u
from astropy.constants import c
from scipy.special import betainc

#from . import MC_errors as MC
from . import MC_errors_general as MC
#from . import fit_chi2 as fit
from . import fit_general as fit
from . import plot_fit

class Line_fit(object):
    
    def __init__(self,_line_name,_line_center,_line_velocity_region,
                 _continuum_regions):
        self.line_name   = _line_name
        self.line_center = _line_center
        self.line_velocity_region = _line_velocity_region
        self.continuum_regions = _continuum_regions

        self.x0_line = None
        self.x0_cont = None
        self.xopt_line = None
        self.xopt_cont = None
        
        #Values from Ftest
        self.chi2 = None
        self.chi2_no_line = None
        self.F = None
        self.p = None


    def run_fit(self, spec):

        #Check that the fit can be run.
        if not self.can_fit_be_run(spec):
            return

        #If no initial guesses have been set, set the default ones
        #here.
        if self.x0_line is None or self.x0_cont is None:
            self.set_initial_fit_values(spec)
        
        #Run the fit.
        self.xopt_line, self.xopt_cont = fit.fit(spec, self)
        self.set_cont_pars(self.xopt_cont)
        self.set_line_pars(self.xopt_line)
        return

    def run_MC(self,spec,nrep,Ncpu=None,save_chain=None):
        if self.xopt_line is None:
            print("First run the line fit.")
            return
        MC.MC_errors(nrep, spec, self, Ncpu=Ncpu, save_chain=save_chain)
        return

    def plot(self,spec,plot_fname=None,chain=None):
        plot_fit.plot_fit(spec,self,chain=chain,plot_fname=plot_fname)

    def run_Ftest(self,spec):
        if self.xopt_line is None:
            print("First fit the line")
            return

        #Only consider regions within a certain velocity range of the
        #canonical emission line center.
        v = (c*(spec.lam_rest/self.line_center-1.)).to(u.km/u.s)
        vabs = np.abs(v)
        iuse = np.argwhere(vabs<self.line_velocity_region)

        
        #With emission line.
        self.chi2 = fit.chi2_line_fit(self.xopt_line,
                                      spec, self,
                                      iuse, self.xopt_cont,
                                      check_constraints=False)

        #Without emission line.
        x_noline = np.copy(self.xopt_line)
        x_noline[1] = 0.
        self.chi2_no_line = fit.chi2_line_fit(x_noline,
                                              spec, self,
                                              iuse, self.xopt_cont,
                                              check_constraints=False)
        
        
        #Get the degrees of freedom.
        n_datapoints = len(iuse)
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
