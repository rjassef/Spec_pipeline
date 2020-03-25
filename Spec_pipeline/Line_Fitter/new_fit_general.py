import numpy as np
import astropy.units as u
from astropy.constants import c
from scipy.optimize import fmin

#####
# This version of fit_general fits simultaneously the continuum and the emission line instead of first fitting the continuum and then the emission lines. This should give more estability to the process.
#####


def chi2_fit(x,spec,line_fitter,iuse,check_constraints=True):

    x_line = x[:line_fitter.npar_line]
    x_cont = x[line_fitter.npar_line:]

    if check_constraints:
        if not line_fitter.meet_cont_constraints(x_cont):
            return np.inf
        if not line_fitter.meet_line_constraints(x_line):
            return np.inf

    #Construct the model
    flam_mod = line_fitter.flam_model(spec.lam_rest[iuse],x_line,x_cont)

    #Get the chi2
    diff = spec.flam[iuse]-flam_mod
    flam_err_use = spec.flam_err[iuse]
    chi2 = np.sum((diff/flam_err_use)**2)

    return chi2


#Main fit module.

def fit(spec, line_fitter, x0_cont=None, x0_line=None):

    if x0_cont is None:
        x0_cont = line_fitter.x0_cont
    if x0_line is None:
        x0_line = line_fitter.x0_line

    #Here we'll actually fit the whole thing together.
    i_cont = line_fitter.get_i_cont(spec)
    i_line = line_fitter.get_i_line(spec)
    i_all  = np.unique(np.concatenate((i_cont,i_line)))

    x0 = np.concatenate((x0_line,x0_cont))
    xopt = fmin(chi2_fit, x0  , args=(spec,line_fitter,i_all),disp=False)
    xopt = fmin(chi2_fit, xopt, args=(spec,line_fitter,i_all),disp=False)
    xopt_line = xopt[:line_fitter.npar_line]
    xopt_cont = xopt[line_fitter.npar_line:]

    return xopt_line, xopt_cont
