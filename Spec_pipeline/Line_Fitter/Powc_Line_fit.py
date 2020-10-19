import numpy as np
import astropy.units as u
from astropy.constants import c
from scipy.special import betainc
import matplotlib.pyplot as plt
import re

from .complex_line import Complex_Line_fit
from .new_fit_general as fit

class Powc_Line_fit(Complex_Line_fit):

    def __init__(self, line_name=None, spec=None):
        super(Powc_Line_fit,self).__init__(line_name, spec)
        return

    #####
    # Continuum Model

    #For now, we define a linear continuum.
    def flam_cont_model(self,lam,x_cont=None):
        a, b = self.cont_par_parser(x_cont)
        return b * (lam/(5000.*u.AA))**a

    def cont_par_parser(self,x_cont):
        if x_cont is None:
            a = self.xopt[-2]
            b = self.xopt[-1]*self.flamunit
        else:
            a = x_cont[0]
            b = x_cont[1]*self.flamunit
        return a, b

    def cont_initial_fit_values(self,spec):
        i_cont = self.get_i_cont(spec)
        mean_cont = np.mean(spec.flam[i_cont]).to(self.flamunit)
        mean_lam  = np.mean(spec.lam_rest[i_cont]).to(self.waveunit)
        a0 = 1.0
        b0 = mean_cont.value * (mean_lam/(5000.*u.AA)).to(1.).value
        x0_cont = [a0,b0]
        return x0_cont

    @property
    def npar_cont(self):
        return 2

    def meet_cont_constraints(self,x_cont):
        a = x_cont[0]
        b = x_cont[1]*self.flamunit
        cont_fluxes = (b*(self._Complex_Line_fit__lam_fit/(5000.*u.AA))**a).to(self.flamunit).value
        if len(cont_fluxes[cont_fluxes<0])>0:
            return False
        return True


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
            print("First fit the lines")
            return

        for line in self.multi_line:

            #Get the indices of the line fitting region.
            iuse = line.get_i_line(spec)

            #With emission line.
            chi2 = fit.chi2_fit(self.xopt, spec, self, iuse, check_constraints=False)

            #Without emission line.
            chi2_no_line = np.zeros(line.nlines)
            flam_line_fit_back = np.copy(line.flam_line_fit)
            for i in range(line.nlines):
                line.flam_line_fit[i] = 0*self.flamunit
                flam_mod = self.flam_cont_model(spec.lam_rest[iuse])
                diff = spec.flam[iuse]-flam_mod
                flam_err_use = spec.flam_err[iuse]
                chi2_no_line[i] = np.sum((diff/flam_err_use)**2)
                line.flam_line_fit[i] = flam_line_fit_back[i]

            #Get the degrees of freedom.
            n_datapoints = len(iuse)
            nu = n_datapoints - line.npar_fit
            nu_no_line = n_datapoints - (line.npar_fit-line.npar_line)
            chi2_nu = chi2/float(nu)
            chi2_no_line_nu = chi2_no_line/float(nu_no_line)

            line.F = ((chi2_no_line-chi2)/float(nu_no_line-nu)) / (chi2/float(nu))
            line.F = line.F.value
            line.p = np.zeros(line.nlines)
            for i in range(line.nlines):
                if line.F[i]<0.:
                    line.p[i] = 1.0
                else:
                    nnu1 = float(nu_no_line-nu)
                    nnu2 = float(nu)
                    w = nnu1*line.F[i]/(nnu1*line.F[i]+nnu2)
                    line.p[i] = 1.-betainc(nnu1/2.,nnu2/2.,w)
        return

    #This is a fast fit doing linear squares in natural log space to make the Gaussians and the power-law continuum linear.
    def run_fast_fit(self):

        spec = self.default_spec

        #Fit the power-law continuum first.
        i_cont = self.get_i_cont(spec)
        lam_c_all  = spec.lam_rest[i_cont]
        flam_c_all = spec.flam[i_cont]
        lam_c  = lam_c_all[flam_c_all>0.*self.flamunit]
        flam_c = flam_c_all[flam_c_all>0.*self.flamunit]

        y = np.log((flam_c/self.flamunit).to(1.).value)
        x = np.log((lam_c/(5000.*u.AA)).to(1.).value)
        A = np.zeros((2,2))
        B = np.zeros(2)
        for i in range(2):
            for j in range(2):
                A[i,j] = np.sum(x**(i+j))
            B[i] = np.sum(x**i * y)
        par_fit = np.linalg.solve(A,B)
        xopt_cont = np.array([par_fit[1],np.exp(par_fit[0])])

        #Now, fit the each of the Gaussian lines.
        xopt_lines = np.zeros(self.ntot_lines*3)
        n = 0
        for line in self.multi_line:
            i_line = line.get_i_line(spec)
            lam_l_all = spec.lam_rest[i_line]
            flam_l_all = spec.flam[i_line]
            lam_l  = lam_l_all[flam_l_all>0.*self.flamunit]
            flam_l = flam_l_all[flam_l_all>0.*self.flamunit]
            flam_l -= self.flam_cont_model(lam_l,xopt_cont)

            for k in range(line.nlines):
                v = c*(lam_l/line.line_center[k]-1.)
                i_line_use = np.argwhere((np.abs(v)<500.*u.km/u.s) & (flam_l>0.*self.flamunit))
                v_use = v[i_line_use]
                lam_l_use = lam_l[i_line_use]
                flam_l_use = flam_l[i_line_use]
                if len(v_use)==0:
                    continue

                y = np.log((flam_l_use/self.flamunit).to(1.).value)
                x = v_use.to(u.km/u.s).value

                A = np.zeros((3,3))
                B = np.zeros(3)
                for i in range(3):
                    for j in range(3):
                        A[i,j] = np.sum(x**(i+j))
                    B[i] = np.sum(x**i * y)
                try:
                    par_fit = np.linalg.solve(A,B)
                except np.linalg.LinAlgError:
                    xopt_lines[n*3:(n+1)*3] = 0.
                    continue
                if par_fit[2]>=0.:
                    xopt_lines[n*3:(n+1)*3] = 0.
                    continue

                dv_fit        = -par_fit[1]/(2.*par_fit[2])
                flam_line_fit = np.exp(par_fit[0]-par_fit[1]**2/(4*par_fit[2]))
                sigma_v_fit   = (-1./(2.*par_fit[2]))**0.5

                #Do not accept lines narrower than the instrumental resolution.
                if sigma_v_fit < line.sigma_v_min[k].to(line.vunit).value:
                    continue

                #Save the best-fit parameters.
                xopt_lines[n*3  ] = dv_fit
                xopt_lines[n*3+1] = flam_line_fit
                xopt_lines[n*3+2] = sigma_v_fit

                #Finally, subtract the emission line from the line flux so a nearby line does not identify it again by mistake.
                dv_fit        = dv_fit        * line.vunit
                flam_line_fit = flam_line_fit * line.flamunit
                sigma_v_fit   = sigma_v_fit   * line.vunit
                flam_l -= flam_line_fit * np.exp((v-dv_fit)/sigma_v_fit)

                n+=1


        self.xopt = np.zeros(len(xopt_cont)+len(xopt_lines))
        self.xopt[:self.npar_line] = xopt_lines
        self.xopt[self.npar_line:] = xopt_cont

        self.set_line_pars(self.xopt[:self.npar_line])

        return

    #Because of the fact that we are fitting each line in a completely independent manner, we actually need to redefine flam_line_model to make sure parameters are assigned appropriately.
    def flam_line_model(self,lam,x_line=None):

        if x_line is None:
            x_line = self.xopt[:self.npar_line]

        #Note that n_multilines counts the number of Multi_Line_fit objects, but each of those can be composed of more than 1 line.
        flam_line_model = np.zeros((self.ntot_lines,len(lam)))*self.flamunit
        k = 0
        for line in self.multi_line:
            for i in range(line.nlines):
                dv        = x_line[k*3]   * line.vunit
                flam_line = x_line[k*3+1] * line.flamunit
                sigma_v   = x_line[k*3+2] * line.vunit

                v = c*(lam/line.line_center[i]-1.)
                flam_line_model[k] = flam_line * np.exp(-0.5*((v-dv)/sigma_v)**2)
                k+=1

        return np.sum(flam_line_model,axis=0)

    @property
    def npar_line(self):
        return self.ntot_lines*3

    def set_line_pars(self,x_line):

        self.dv_fit        = self.xopt[0:self.npar_line:3]*self.vunit
        self.flam_line_fit = self.xopt[1:self.npar_line:3]*self.flamunit
        self.sigma_v_fit   = self.xopt[2:self.npar_line:3]*self.vunit

        k = 0
        for line in self.multi_line:
            line.dv_fit        = np.zeros(line.nlines)*line.vunit
            line.flam_line_fit = np.zeros(line.nlines)*line.flamunit
            line.sigma_v_fit   = np.zeros(line.nlines)*line.vunit
            for i in range(line.nlines):
                line.dv_fit[i]        = self.dv_fit[k]
                line.flam_line_fit[i] = self.flam_line_fit[k]
                line.sigma_v_fit[i]   = self.sigma_v_fit[k]
                k+=1

        return
