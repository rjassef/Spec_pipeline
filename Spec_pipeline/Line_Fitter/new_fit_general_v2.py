import numpy as np
import astropy.units as u
from astropy.constants import c
from scipy.optimize import fmin
from scipy.optimize import minimize, NonlinearConstraint

#####
# This version of fit_general fits simultaneously the continuum and the emission line instead of first fitting the continuum and then the emission lines. This should give more estability to the process.
#####


def chi2_fit(x,spec,line_fitter,iuse,check_constraints=True):

    if check_constraints:
        if not line_fitter.meet_constraints(x):
            return np.inf

    #Construct the model
    flam_mod = line_fitter.flam_model(spec.lam_rest[iuse],x)

    #Get the chi2
    diff = spec.flam[iuse]-flam_mod
    flam_err_use = spec.flam_err[iuse]
    chi2 = np.sum((diff/flam_err_use)**2)

    return chi2


#Main fit module.

def fit(spec, line_fitter, x0=None):

    if x0 is None:
        x0 = line_fitter.x0

    #Here we'll actually fit the whole thing together.
    i_all = line_fitter.get_i_fit(spec)

    #We'll try with only setting a non-linear constraint, not a bound. 
    con = lambda x: int(not line_fitter.meet_constraints(x))
    nlc = NonlinearConstraint(con, -np.inf, -0.5)
    #print(nlc)

    print(x0)
    print(con(x0))

    #Fit
    res = minimize(chi2_fit, x0  , args=(spec,line_fitter,i_all,False), constraints=[{'type':'eq', 'fun':con}], method='SLSQP')
    print(res)
    xopt = res.x
    print(xopt, chi2_fit(xopt,spec,line_fitter,i_all, True))

    xopt = fmin(chi2_fit, x0  , args=(spec,line_fitter,i_all,True),disp=False)
    xopt = fmin(chi2_fit, xopt, args=(spec,line_fitter,i_all,True),disp=False)
    print(xopt, chi2_fit(xopt,spec,line_fitter,i_all))

    return xopt
