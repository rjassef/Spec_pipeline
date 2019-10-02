#!/usr/bin/env python 

import numpy as np
import astropy.units as u
from scipy.optimize import fmin
from pysynphot import observation
from pysynphot import spectrum
from astropy.constants import h,c

###

def rolling_linear_regression(lam,flam,window):

    lam_lr      = np.zeros(len(lam)-window)
    flam_lr     = np.zeros(len(lam)-window)
    flam_std_lr = np.zeros(len(lam)-window)
    for k in range(len(lam)-window):
        lam_use  = lam[k:k+window]
        flam_use = flam[k:k+window]

        f  = flam_use.mean()
        l  = lam_use.mean()
        l2 = (lam_use*lam_use).mean()
        fl = (flam_use*lam_use).mean()

        a = (fl-f*l)/(l2-l**2)
        b = f-a*l

        flam_mod = a*lam_use + b
        
        kk = int((window-1)/2)
        lam_lr[k]  = lam_use[kk].value
        flam_lr[k] = flam_mod[kk].value
        flam_std_lr[k] = (flam_use-flam_mod).std().value

    lam_lr *= lam.unit
    flam_lr *= flam.unit
    flam_std_lr *= flam.unit

    return lam_lr, flam_lr, flam_std_lr

###

def S_func(x,lam_obs,flam,SN_lam,flam_sky,eps,RON):

    K1 = np.abs(x[0])
    K2 = np.abs(x[1])

    #if K1<0.1 or K1>100.:
    #    return 1e32
    #if K2<0.1 or K2>100.:
    #if K2<0.:
    #    return 1e32
    
    top = K1*eps*np.abs(flam)*lam_obs
    bot = K1*eps*(np.abs(flam)+flam_sky*K2)*lam_obs + RON**2
    SN_lam_mod = top/(bot**0.5)

    S = np.sum((SN_lam-SN_lam_mod)**2)
    return S
    

###

def get_error_pars(lam_obs,flam,SN_lam,flam_sky,eps, RON):

    #Both K1 and K2 should be around 1.0
    K1_0 = 1.0
    K2_0 = 0.0
    x0 = np.array([K1_0, K2_0])
    xopt = fmin(S_func,x0  ,args=(lam_obs,flam,SN_lam,flam_sky,eps,RON),
                                  disp=False)
    xopt = fmin(S_func,xopt,args=(lam_obs,flam,SN_lam,flam_sky,eps,RON),
                                  disp=False)

    K1 = np.abs(xopt[0])
    K2 = np.abs(xopt[1])
   
    return K1, K2

###

def get_error_spec(spec, wd=15):

    #Calculate the rolling mean and std of the spectrum in a window of
    #wd pixels.
    lam_mean, flam_mean, flam_std = rolling_linear_regression(
        spec.lam_obs,spec.flam,wd)

    #Repeat with the sky spectrum and the eps factor.
    lam_sky_mean, flam_sky_mean, flam_sky_std = rolling_linear_regression(
        spec.lam_obs,spec.flam_sky,wd)
    if np.isscalar(spec.eps.value):
        eps_use = spec.eps
    else:
        lam_eps_use, eps_use, eps_std = rolling_linear_regression(
            spec.lam_obs,spec.eps,wd)
        
        
    #Get the zero-th order SN array
    SN_lam = flam_mean/flam_std

    #Fit the SN array to the error parameters.
    K1, K2 = get_error_pars(lam_mean,flam_mean,SN_lam,flam_sky_mean,
                            eps_use,spec.RON)

    #Get the error spectrum
    flam_err = np.sqrt(K1*spec.eps*(np.abs(spec.flam)+
                                    K2*spec.flam_sky)*spec.lam_obs + 
                       spec.RON**2)/(K1*spec.eps*spec.lam_obs)

    return flam_err.to(u.erg/u.s/u.cm**2/u.AA), K1, K2


