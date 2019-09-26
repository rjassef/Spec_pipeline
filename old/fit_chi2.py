import numpy as np
import astropy.units as u
from astropy.constants import c
from scipy.optimize import fmin

###

#Set the constraints. We'll set the constraints for CIV now, but can
#be changed as required.
sigma_v_min = 100.*u.km/u.s
sigma_v_max = 5000.*u.km/u.s

delta_lam_cen_max = 2.*u.AA

frac_flam_line_cen_min = 1e-1 #Minimum peak height fraction compared
                              #to initial guess.

###

def join_units(a,b):
    aux = np.concatenate((a.value,b.value))
    aux = aux*a.unit
    return aux

###

def flam_line_model(lam_cen, flam_line_cen, sigma_v, lam):
    v  = c*(lam/lam_cen-1.)
    return flam_line_cen * np.exp(-0.5*(v/sigma_v)**2)
    
def cont_model(a, b, lam):
    return a*lam + b
    
def flam_model(x,a,b,lam):
    
    lam_cen       = x[0] * u.AA
    flam_line_cen = x[1] * u.erg/u.s/u.cm**2/u.AA
    sigma_v       = x[2] * u.km/u.s

    flam_mod  = cont_model(a,b,lam)
    flam_mod += flam_line_model(lam_cen,flam_line_cen,sigma_v,lam)
    return flam_mod

def chi2_fit(x,flam,flam_err,lam,a,b,vmax,lam_cen_0,flam_line_cen_0):

    lam_cen       = x[0] * u.AA
    flam_line_cen = x[1] * u.erg/u.s/u.cm**2/u.AA
    sigma_v       = x[2] * u.km/u.s

    #Apply constrains.
    if sigma_v<sigma_v_min or sigma_v>sigma_v_max:
        return np.inf

    if np.abs(lam_cen-lam_cen_0)>delta_lam_cen_max:
        return np.inf

    #Do not allow absorption lines.
    if flam_line_cen<0*u.erg/u.s/u.cm**2/u.AA:
        return np.inf

    #Do not allow the peak of the emission line to drift much below
    #the guess value. There is a failure mode on which the lines
    #are fit at any central wavelength with any width but with a flux
    #of effectively 0.
    if flam_line_cen_0>0 and \
            (flam_line_cen/flam_line_cen_0).to(1.)<frac_flam_line_cen_min:
        return np.inf
    

    #Construct the model.
    flam_mod = flam_model(x,a,b,lam).to(u.erg/u.s/u.cm**2/u.AA)

    #Only consider regions within a certain velocity range of the
    #emission line.
    v = (c*(lam/lam_cen-1.)).to(u.km/u.s)
    vabs = np.abs(v)
    diff = flam[vabs<vmax] - flam_mod[vabs<vmax]

    #Get the chi2
    chi2 = np.sum((diff/flam_err[vabs<vmax])**2)

    return chi2

###

def Continuum_fit(lam_rest,flam,flam_err,continuum_regions):

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

    return a,b

####

#Main fit module. 

def fit(lam_rest, flam, flam_err, continuum_regions, line_velocity_region, 
        lam_cen_0, sigma_v_0):

    ################
    #Chi2 Continuum fit
    ################
    a,b = Continuum_fit(lam_rest,flam,flam_err,continuum_regions)

    ###########
    # Line fit
    ###########
    #Initial guesses
    flam_line_cen_0 = (np.max(flam[np.abs(lam_rest-lam_cen_0)<3*u.AA]) - \
        cont_model(a, b, lam_cen_0))*0.5

    xini = [lam_cen_0.to(u.AA).value,
            flam_line_cen_0.to(u.erg/u.s/u.cm**2/u.AA).value,
            sigma_v_0.to(u.km/u.s).value]
    xopt = fmin(chi2_fit,xini,args=(flam,flam_err,lam_rest,a,b,
                                    line_velocity_region,lam_cen_0,
                                    flam_line_cen_0), 
                disp=False)
    xopt = fmin(chi2_fit,xopt,args=(flam,flam_err,lam_rest,a,b,
                                    line_velocity_region,lam_cen_0,
                                    flam_line_cen_0), 
                disp=False)

    lam_cen       = xopt[0] * u.AA
    flam_line_cen = xopt[1] * u.erg/u.s/u.cm**2/u.AA
    sigma_v       = xopt[2] * u.km/u.s

    flam_mod = flam_model(xopt,a,b,lam_rest)    

    return lam_cen, flam_line_cen, sigma_v, a, b, flam_mod

