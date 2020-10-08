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

        """
        This function returns an updated redshift for each line being fit according to its best-fit velocity offset.

        Parameters
        ----------
        spec_use : Spec object, optional
            Spectrum being fit.

        """

        spec = self.get_spec_use(spec_use)
        if spec is None:
            return
        #return spec.zspec+self.dv_fit/c
        return (1+spec.zspec)*(1+self.dv_fit/c)-1

    def run_fit(self, spec_use=None):

        """
        This function runs the fit.

        Parameters
        ----------
        spec_use : Spec object, optional
            Spectrum that will be fit.

        """

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

        """
        This function runs the MC. Needs to be run after run_fit.

        Parameters
        ----------
        nrep : int
            Number of MC resamples.

        spec_use : Spec object, optional
            Spectrum that will be fit.

        Ncpu : int, optional
            Number of CPU cores to be used. Default is to use all available.

        save_chain : file path, optional
            File name to save the MC chain. If none given, chain is not saved to a file.

        """

        spec = self.get_spec_use(spec_use)
        if spec is None:
            return
        if self.xopt is None:
            print("First run the line fit.")
            return
        MC.MC_errors(nrep, spec, self, Ncpu=Ncpu, save_chain=save_chain)
        return

    def plot(self,spec_use=None,plot_fname=None,chain_file=None,chain=None):

        """
        This function plots the best-fit model to the spectrum within the fitting ranges.

        Parameters
        ----------
        spec_use : Spec object, optional
            Spectrum that will be fit.

        plot_fname : file path, optional
            Filename to save a hard copy of the plot. If none provided, the plot is shown to screen.

        chain_file : file path, optional
            MC chain file. Overrides the default and the current chain obtained from run_MC if any.

        chain : 2D array of shape nrep x npar, optional
            If provided, uses this one instead of that calculated by run_MC for this plot.

        """
        spec = self.get_spec_use(spec_use)
        if spec is None:
            return
        if chain is None:
            chain = self.MC_chain
        plot_fit.plot_fit(spec,self,chain=chain,chain_file=chain_file,plot_fname=plot_fname)

    def run_Ftest(self,spec_use=None):

        """
        This function runs an F-test. Can only be run after the fit. Results are saved in self.p and self.F .

        Parameters
        ----------
        spec_use : Spec object, optional
            Spectrum that will be fit.

        """

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
        self.chi2_no_line = np.zeros(self.nlines)
        flam_line_fit_back = np.copy(self.flam_line_fit)
        for i in range(self.nlines):
            self.flam_line_fit[i] = 0*self.flamunit
            flam_mod = self.flam_model(spec.lam_rest[iuse])
            diff = spec.flam[iuse]-flam_mod
            flam_err_use = spec.flam_err[iuse]
            self.chi2_no_line[i] = np.sum((diff/flam_err_use)**2)
            self.flam_line_fit[i] = flam_line_fit_back[i]

        #x_noline = np.copy(self.xopt)
        #x_noline[:self.npar_line] = 0
        #self.chi2_no_line = fit.chi2_fit(x_noline, spec, self, iuse, check_constraints=False)


        #Get the degrees of freedom.
        n_datapoints = len(iuse)
        nu = n_datapoints - self.npar_fit
        nu_no_line = n_datapoints - (self.npar_fit-self.npar_line)
        chi2_nu = self.chi2/float(nu)
        chi2_no_line_nu = self.chi2_no_line/float(nu_no_line)

        self.F = ((self.chi2_no_line-self.chi2)/float(nu_no_line-nu)) / \
                 (self.chi2/float(nu))
        self.F = self.F.value
        self.p = np.zeros(self.nlines)
        for i in range(self.nlines):
            if self.F[i]<0.:
                self.p[i] = 1.0
            else:
                nnu1 = float(nu_no_line-nu)
                nnu2 = float(nu)
                w = nnu1*self.F[i]/(nnu1*self.F[i]+nnu2)
                self.p[i] = 1.-betainc(nnu1/2.,nnu2/2.,w)
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
            'flamunit',
            'nlines'
            ]

        for att in list_of_attibutes:
            try:
                eval("self."+att)
            except AttributeError as err:
                print(err)
        return
