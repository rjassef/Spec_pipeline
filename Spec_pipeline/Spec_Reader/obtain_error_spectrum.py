#!/usr/bin/env python 

import numpy as np
import astropy.units as u
from scipy.optimize import fmin
from scipy import interpolate
from pysynphot import observation
from pysynphot import spectrum
from astropy.constants import h,c
import matplotlib.pyplot as plt

###

def rolling_linear_regression(lam,flam,window):

    lam_lr      = np.zeros(len(lam)-window)*lam.unit
    flam_lr     = np.zeros(len(lam)-window)*flam.unit
    flam_std_lr = np.zeros(len(lam)-window)*flam.unit
    for k in range(len(lam)-window):
        lam_use  = lam[k:k+window]
        flam_use = flam[k:k+window]

        ndeg = 3
        p = np.polyfit(lam_use.value,flam_use.value,ndeg)
        flam_mod = 0
        for i in range(ndeg+1):
            flam_mod += p[ndeg-i]*(lam_use.value)**i
        flam_mod = flam_mod*flam_use.unit

        kk = int((window-1)/2)
        lam_lr[k]  = lam_use[kk]
        flam_lr[k] = flam_mod[kk]
        flam_std_lr[k] = np.std(flam_use-flam_mod)

    return lam_lr, flam_lr, flam_std_lr

###

#https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html

def smooth(x,wd):
    w = np.ones(wd,'d')
    s = np.r_[x[wd-1:0:-1],x,x[-2:-wd-1:-1]]
    y = np.convolve(w/w.sum(),s,mode='valid')
    return y

def smooth_errors(lam,flam,flam_sky,sens,wd):

    #lam_lr, flam_lr, flam_std_lr = rolling_linear_regression(lam,flam,wd)
    #aux1  , flam_sky_lr, aux2    = rolling_linear_regression(lam,flam_sky,wd)
    #aux1  , sens_lr    , aux2    = rolling_linear_regression(lam,sens,wd)
    #return lam_lr, flam_lr, flam_std_lr, flam_sky_lr, sens_lr

    lam_sm      = smooth(lam.value,wd)
    flam_sm     = smooth(flam.value,wd)
    flam_sky_sm = smooth(flam_sky.value,wd)
    sens_sm     = smooth(sens,wd)

    flam_mod    = np.interp(lam.value,lam_sm,flam_sm)*flam.unit
    flam_sky_mod= np.interp(lam.value,lam_sm,flam_sky_sm)*flam_sky.unit
    sens_mod    = np.interp(lam.value,lam_sm,sens_sm)

    lam_lr      = np.zeros(len(lam)-wd)*lam.unit
    flam_lr     = np.zeros(len(lam)-wd)*flam.unit
    flam_std_lr = np.zeros(len(lam)-wd)*flam.unit
    flam_sky_lr = np.zeros(len(lam)-wd)*flam.unit
    sens_lr     = np.zeros(len(lam)-wd)
    for k in range(len(lam)-wd):
        lam_use  = lam[k:k+wd]
        flam_use = flam[k:k+wd]
        flam_mod_use = flam_mod[k:k+wd]
        flam_sky_mod_use = flam_sky_mod[k:k+wd]
        sens_mod_use = sens_mod[k:k+wd]

        kk = int((wd-1)/2)
        lam_lr[k]  = lam_use[kk]
        flam_lr[k] = flam_mod_use[kk]
        flam_std_lr[k] = np.std(flam_use-flam_mod_use)
        flam_sky_lr[k] = flam_sky_mod_use[kk]
        sens_lr[k]     = sens_mod_use[kk]

    return lam_lr, flam_lr, flam_std_lr, flam_sky_lr, sens_lr

###

def S_func(x,flam,flam_std,flam_sky,eps,RON):

    K1 = x[0]
    K2 = x[1]

    if K1<0. or K2<0.:
        return np.inf
    
    flam_std_mod = np.sqrt( K1*eps*(flam+K2*flam_sky) + RON**2 ) / (K1*eps)

    S = np.sum((flam_std-flam_std_mod)**2).value
    return S
    

###

def get_error_pars(flam,SN_lam,flam_sky,eps,RON):

    #Both K1 and K2 should be around 1.0
    #K1_0 = 1e-3
    #K2_0 = 1e-3
    K1_0 = 1.0
    K2_0 = 1e-3
    x0 = np.array([K1_0, K2_0])
    xopt = fmin(S_func,x0  ,args=(flam,SN_lam,flam_sky,eps,RON),
                                  disp=False)
    xopt = fmin(S_func,xopt,args=(flam,SN_lam,flam_sky,eps,RON),
                                  disp=False)

    K1 = xopt[0]
    K2 = xopt[1]
   
    return K1, K2

###

def get_error_spec(spec, wd=15):

    #Calculate the rolling mean and std of the spectrum in a window of
    #wd pixels.
    lam_meanx, flam_meanx, flam_stdx, \
        flam_sky_meanx, sens_usex = \
                                    smooth_errors(
                                        spec.lam_obs,spec.flam,
                                        spec.flam_sky,spec.sens,wd)

    
    #For this exercise, we need to not consider the Lyalpha forrest
    #region. Real IGM absorption appears as noise and throws
    #everything off the board.
    kuse = np.where(lam_meanx>1300.*u.AA*(1.+spec.zspec))
    lam_mean = lam_meanx[kuse]
    flam_mean = flam_meanx[kuse]
    flam_std = flam_stdx[kuse]
    flam_sky_mean = flam_sky_meanx[kuse]
    sens_use = sens_usex[kuse]

    #Estimate the size of the new wavelength bin.
    dlam_mean = np.mean(lam_mean[1:]-lam_mean[:-1])

    #Get the eps factor to use.
    eps_use = spec.eps(sens=sens_use,lam_obs=lam_mean,dlam=dlam_mean)

    #Fit the std array to the error parameters.
    K1, K2 = get_error_pars(flam_mean,flam_std,flam_sky_mean,
                            eps_use,spec.RON)

    #Get the error spectrum
    flam_err = np.sqrt(K1*spec.eps()*(np.abs(spec.flam)+
                                    K2*spec.flam_sky) + 
                       spec.RON**2)/(K1*spec.eps())

    return flam_err.to(u.erg/(u.s*u.cm**2*u.AA)), K1, K2


