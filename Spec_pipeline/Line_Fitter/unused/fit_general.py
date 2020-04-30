import numpy as np
import astropy.units as u
from astropy.constants import c
from scipy.optimize import fmin


def chi2_line_fit(x,spec,line_fitter,iuse,x_cont,
                  check_constraints=True):

    #Check the strict constraints are met.
    if check_constraints:
        if not line_fitter.meet_line_constraints(x):
            return np.inf

    #Construct the model.
    flam_mod = line_fitter.flam_model(spec.lam_rest[iuse],x,x_cont)

    #Get the chi2
    diff = spec.flam[iuse]-flam_mod
    flam_err_use = spec.flam_err[iuse]
    chi2 = np.sum((diff/flam_err_use)**2)

    return chi2

def chi2_cont_fit(x,spec,line_fitter,iuse,check_constraints=True):

    #Check the strict constraints are met.
    if check_constraints:
        if not line_fitter.meet_cont_constraints(x):
            return np.inf

    #Construct the model.
    flam_mod = line_fitter.flam_cont_model(spec.lam_rest[iuse],x)

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
    
    #################
    # Continuum fit.
    #################

    #Define the indices of the continuum regions.
    i_cont = line_fitter.get_i_cont(spec)

    #Run the fit
    xopt_cont = fmin(chi2_cont_fit, line_fitter.x0_cont,
                     args=(spec,line_fitter,i_cont),disp=False)
    xopt_cont = fmin(chi2_cont_fit, xopt_cont,
                     args=(spec,line_fitter,i_cont),disp=False)

    #################
    # Line fit.
    #################
    
    #Determine the indices of the line fitting region.
    i_line = line_fitter.get_i_line(spec)

    #Run the line fit.
    xopt_line = fmin(chi2_line_fit, line_fitter.x0_line  ,
                     args=(spec,line_fitter,i_line,xopt_cont),disp=False)
    xopt_line = fmin(chi2_line_fit, xopt_line,
                     args=(spec,line_fitter,i_line,xopt_cont),disp=False)

    return xopt_line, xopt_cont
