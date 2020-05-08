# Implementation of the line fitting routines. This is the main class,
# we will have a separate class with different defaults for different
# emission lines.

import numpy as np
from scipy.special import betainc
from astropy.constants import c

from . import MC_errors_general as MC
from . import new_fit_general as fit
from . import plot_fit

class Line_fit(object):

    def __init__(self,_line_name,_default_spec=None):

        self.line_name = _line_name
        self.default_spec = _default_spec

        #Fit arrays
        self.x0 = None
        self.xopt = None

        self.nlines = None

        #Values from Ftest
        self.chi2 = None
        self.chi2_no_line = None
        self.F = None
        self.p = None

        #MC Output chain
        self.MC_chain = None

        return

    def get_spec_use(self,spec):
        if spec is None:
            if self.default_spec is None:
                print("Need to provide an spectrum")
                return None
            else:
                return self.default_spec
        return spec

    def zline(self, spec_use=None):
        spec = self.get_spec_use(spec_use)
        if spec is None:
            return
        #return spec.zspec+self.dv_fit/c
        return (1+spec.zspec)*(1+self.dv_fit/c)-1

    def run_fit(self, spec_use=None):

        spec = self.get_spec_use(spec_use)
        if spec is None:
            return

        #Check that the fit can be run.
        if not self.can_fit_be_run(spec):
            return

        #If no initial guesses have been set, set the default ones
        #here.
        if self.x0 is None:
            self.set_initial_fit_values(spec)

        #Run the fit.
        self.xopt = fit.fit(spec, self)
        self.set_pars(self.xopt)
        return

    def run_MC(self,nrep,spec_use=None,Ncpu=None,save_chain=None):
        spec = self.get_spec_use(spec_use)
        if spec is None:
            return
        if self.xopt is None:
            print("First run the line fit.")
            return
        MC.MC_errors(nrep, spec, self, Ncpu=Ncpu, save_chain=save_chain)
        return

    def plot(self,spec_use=None,plot_fname=None,chain_file=None,chain=None):
        spec = self.get_spec_use(spec_use)
        if spec is None:
            return
        if chain is None:
            chain = self.MC_chain
        plot_fit.plot_fit(spec,self,chain=chain,chain_file=chain_file,plot_fname=plot_fname)

    def run_Ftest(self,spec_use=None):

        spec = self.get_spec_use(spec_use)
        if spec is None:
            return

        if self.xopt is None:
            print("First fit the line")
            return

        #Get the indices of the line fitting region.
        iuse = self.get_i_line(spec)

        #With emission line.
        self.chi2 = fit.chi2_fit(self.xopt, spec, self, iuse, check_constraints=False)

        #Without emission line.
        x_noline = np.copy(self.xopt)
        x_noline[:self.npar_line] = 0
        self.chi2_no_line = fit.chi2_line_fit(x_noline, spec, self, iuse, check_constraints=False)


        #Get the degrees of freedom.
        n_datapoints = len(iuse)
        nu = n_datapoints - self.npar_line
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

        ###
        # In order for a subclass to work, a significant number of functions need to be defined within it. This function just checks that the subclass is ready to be used. Should only be called for debugging purposes.
    def check_subclass_ready(self):

        list_of_attibutes = [
            'dv_fit',
            'can_fit_be_run',
            'set_initial_fit_values',
            'set_pars',
            'get_i_fit',
            'get_i_line',
            'get_i_cont',
            'npar_fit',
            'npar_line',
            'meet_constraints',
            'flam_model',
            'flam_cont_model',
            'parse_chain_output',
            'ncont_reg',
            'continuum_regions',
            'flamunit'
            ]

        for att in list_of_attibutes:
            try:
                eval("self."+att)
            except AttributeError as err:
                print(err)
        return
