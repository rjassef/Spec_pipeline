import numpy as np
import astropy.units as u
from astropy.constants import c
from scipy.optimize import fmin

###

def join_units(a,b):
    aux = np.concatenate((a.value,b.value))
    aux = aux*a.unit
    return aux

###

def flam_model(x,a,b,lam,line_fitter):
    
    lam_cen       = x[0] * u.AA
    flam_line_cen = x[1] * u.erg/u.s/u.cm**2/u.AA
    sigma_v       = x[2] * u.km/u.s

    flam_mod  = line_fitter.flam_cont_model(lam,a,b)
    flam_mod += line_fitter.flam_line_model(lam,
                                            lam_cen,flam_line_cen,sigma_v)
    return flam_mod

def meet_constraints(x,line_fitter,lam_cen_0,flam_line_cen_0):
    
    lam_cen       = x[0] * u.AA
    flam_line_cen = x[1] * u.erg/u.s/u.cm**2/u.AA
    sigma_v       = x[2] * u.km/u.s
   

    if sigma_v<line_fitter.sigma_v_min or sigma_v>line_fitter.sigma_v_max:
        return False

    if np.abs(lam_cen-lam_cen_0)>line_fitter.delta_lam_cen_max:
        return False

    #Do not allow absorption lines.
    if flam_line_cen<0*u.erg/u.s/u.cm**2/u.AA:
        return False

    #Do not allow the peak of the emission line to drift much below
    #the guess value. There is a failure mode on which the lines
    #are fit at any central wavelength with any width but with a flux
    #of effectively 0.
    if flam_line_cen_0>0 and \
            (flam_line_cen/flam_line_cen_0).to(1.)<\
            line_fitter.frac_flam_line_cen_min:
        return False

    return True


def chi2_fit(x,spec,line_fitter,a,b,lam_cen_0,flam_line_cen_0,
             check_constraints=True):

    lam_cen       = x[0] * u.AA
    flam_line_cen = x[1] * u.erg/u.s/u.cm**2/u.AA
    sigma_v       = x[2] * u.km/u.s

    #Check the strict constraints are met.
    if check_constraints:
        if not meet_constraints(x,line_fitter,lam_cen_0,flam_line_cen_0):
            return np.inf

    #Construct the model.
    flam_mod = flam_model(x,a,b,spec.lam_rest,line_fitter)
    flam_mod = flam_mod.to(u.erg/u.s/u.cm**2/u.AA)

    #Only consider regions within a certain velocity range of the
    #emission line.
    v = (c*(spec.lam_rest/lam_cen-1.)).to(u.km/u.s)
    vabs = np.abs(v)
    diff = spec.flam[vabs<line_fitter.line_velocity_region] - \
           flam_mod[vabs<line_fitter.line_velocity_region]

    #Get the chi2
    flam_err_use = spec.flam_err[vabs<line_fitter.line_velocity_region]
    chi2 = np.sum((diff/flam_err_use)**2)

    return chi2

###

'''def Continuum_fit(spec,continuum_regions):

    lam_rest = spec.lam_rest
    flam = spec.flam
    flam_err = spec.flam_err

    #Define the continuum regions.
    for i in range(len(continuum_regions[:][0])):
        lam1 = continuum_regions[i][0]
        lam2 = continuum_regions[i][1]
        lam_cont_fit_aux = lam_rest[(lam_rest>=lam1) & (lam_rest<=lam2)]
        flam_cont_fit_aux = flam[(lam_rest>=lam1) & (lam_rest<=lam2)]
        flam_err_cont_fit_aux = flam_err[(lam_rest>=lam1) & (lam_rest<=lam2)]
        if i==0:
            lam_cont_fit      = lam_cont_fit_aux
            flam_cont_fit     = flam_cont_fit_aux
            flam_err_cont_fit = flam_err_cont_fit_aux
        else:
            lam_cont_fit      = join_units(lam_cont_fit, lam_cont_fit_aux)
            flam_cont_fit     = join_units(flam_cont_fit, flam_cont_fit_aux)
            flam_err_cont_fit = join_units(flam_err_cont_fit, 
                                           flam_err_cont_fit_aux)

    #Do a least squares fit.
    flame2 = flam_err_cont_fit**2
    o  = np.sum(1./flame2)
    f  = np.sum(flam_cont_fit/flame2)
    l  = np.sum(lam_cont_fit/flame2)
    l2 = np.sum(lam_cont_fit**2/flame2)
    fl = np.sum(flam_cont_fit*lam_cont_fit/flame2)
    a  = (fl*o-f*l)/(l2*o-l**2)
    b  = (f - a*l)/o

    return a,b'''

####

def Continuum_fit(spec,line_fitter):

    #Define the indices of the continuum regions.
    i_cont = np.argwhere(
        ((spec.lam_rest>=line_fitter.continuum_regions[0][0]) &
         (spec.lam_rest<=line_fitter.continuum_regions[0][1])) |
        ((spec.lam_rest>=line_fitter.continuum_regions[1][0]) &
         (spec.lam_rest<=line_fitter.continuum_regions[1][1])))
    lam_cont_fit = spec.lam_rest[i_cont]
    flam_cont_fit = spec.flam[i_cont]
    flam_err_cont_fit = spec.flam_err[i_cont]
    
    #Do a least squares fit.
    flame2 = flam_err_cont_fit**2
    o  = np.sum(1./flame2)
    f  = np.sum(flam_cont_fit/flame2)
    l  = np.sum(lam_cont_fit/flame2)
    l2 = np.sum(lam_cont_fit**2/flame2)
    fl = np.sum(flam_cont_fit*lam_cont_fit/flame2)
    a  = (fl*o-f*l)/(l2*o-l**2)
    b  = (f - a*l)/o

    return a,b

####

#Main fit module. 

def fit(spec, line_fitter, lam_cen_0, sigma_v_0):

    ################
    #Chi2 Continuum fit
    ################
    a,b = Continuum_fit(spec,line_fitter)

    ###########
    # Line fit
    ###########
    #Initial guesses
    flam_line_cen_0 = \
        (np.max(spec.flam[np.abs(spec.lam_rest-lam_cen_0)<3*u.AA]) - \
        line_fitter.flam_cont_model(lam_cen_0, a, b))*0.5

    xini = [lam_cen_0.to(u.AA).value,
            flam_line_cen_0.to(u.erg/u.s/u.cm**2/u.AA).value,
            sigma_v_0.to(u.km/u.s).value]
    xopt = fmin(chi2_fit,xini,args=(spec,line_fitter,a,b,
                                    lam_cen_0,
                                    flam_line_cen_0), 
                disp=False)
    xopt = fmin(chi2_fit,xopt,args=(spec,line_fitter,a,b,
                                    lam_cen_0,
                                    flam_line_cen_0), 
                disp=False)

    lam_cen       = xopt[0] * u.AA
    flam_line_cen = xopt[1] * u.erg/u.s/u.cm**2/u.AA
    sigma_v       = xopt[2] * u.km/u.s

    flam_mod = flam_model(xopt,a,b,spec.lam_rest,line_fitter)    

    return lam_cen, flam_line_cen, sigma_v, a, b, flam_mod

