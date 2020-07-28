#!/usr/bin/env python

import numpy as np
import astropy.units as u
from scipy.optimize import fmin
from scipy import interpolate
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

#Taken from Section 3 of http://www.ucolick.org/~bolte/AY257/s_n.pdf . We will neglect the Dark Noise. Note that the RON used here is an effective RON:
#
# RON_eff = ( (RON**2 + (GAIN/2.)**2) * extraction_aperture_pixels )**(0.5)
def flam_std_func(Kp,flam,flam_sky,lam,RON_eff,eps):

    flam_blue_exp = Kp[2]*np.exp(Kp[3]*(lam-5000.*u.AA)/( 1000.*u.AA)) / eps
    flam_red_exp  = Kp[4]*np.exp(Kp[5]*(lam-5000.*u.AA)/(10000.*u.AA)) / eps
    flam_sky_use  = Kp[1]*flam_sky + flam_blue_exp + flam_red_exp
    return np.sqrt( Kp[0]*eps*(np.abs(flam)+flam_sky_use) + RON_eff**2 ) / (Kp[0]*eps)

    #flam_poly = Ap*np.exp(Bp*(lam-5000.*u.AA)/(1000.*u.AA)) / eps
    #flam_sky_plus_poly = K2*flam_sky + flam_poly
    #return np.sqrt( K1*eps*(np.abs(flam)+flam_sky_plus_poly) + RON_eff**2 ) / (K1*eps)

###

def S_func(x,flam,flam_std,flam_sky,lam,eps,RON_eff,fit_blue_exp,fit_red_exp):

    Kp = np.zeros(6)
    Kp[0] = x[0] #K1
    Kp[1] = x[1] #K2
    kk = 0
    if fit_blue_exp:
        Kp[2] = x[2] #Ab
        Kp[3] = x[3] #Bb
        kk+=2
    if fit_red_exp:
        Kp[4] = x[2+kk] #Ar
        Kp[5] = x[3+kk] #Br

    if Kp[0]<0. or Kp[1]<0. or Kp[2]<0. or Kp[3]>0. or Kp[4]<0. or Kp[5]<0.:
        return np.inf

    flam_std_mod = flam_std_func(Kp,flam,flam_sky,lam,RON_eff,eps)
    flam_std_mod.to(flam_std.unit)

    S = np.sum((flam_std-flam_std_mod)**2).value
    return S


###

def get_error_pars(flam,flam_std,flam_sky,lam,eps,RON_eff,fit_blue_exp,fit_red_exp):

    npar = 2
    if fit_blue_exp:
        npar += 2
    if fit_red_exp:
        npar += 2
    Kp_0 = np.zeros(npar)

    Kp_0[0] = 1.0    #K1
    Kp_0[1] = 1.e-3  #K2
    kk = 0
    if fit_blue_exp:
        Kp_0[2] =  1.0 #Ab
        Kp_0[3] = -10. #Bb
        kk += 2
    if fit_red_exp:
        Kp_0[2+kk] = 1.0 #Ar
        Kp_0[3+kk] = 10. #Br

    Kp_opt = fmin(S_func,Kp_0  , args=(flam,flam_std,flam_sky,lam,eps,RON_eff,fit_blue_exp,fit_red_exp),
                                  disp=False)
    Kp_opt = fmin(S_func,Kp_opt, args=(flam,flam_std,flam_sky,lam,eps,RON_eff,fit_blue_exp,fit_red_exp),
                                  disp=False)

    Kp = np.zeros(6)
    Kp[0] = Kp_opt[0] #K1
    Kp[1] = Kp_opt[1] #K2
    kk = 0
    if fit_blue_exp:
        Kp[2] = Kp_opt[2] #Ab
        Kp[3] = Kp_opt[3] #Bb
        kk+=2
    if fit_red_exp:
        Kp[4] = Kp_opt[2+kk] #Ar
        Kp[5] = Kp_opt[3+kk] #Br

    return Kp

###

def get_error_spec(spec, wd=15, show_plot=False):

    #Calculate the rolling mean and std of the spectrum in a window of
    #wd pixels.
    lam_meanx, flam_meanx, flam_stdx, \
        flam_sky_meanx, sens_usex = \
                                    smooth_errors(
                                        spec.lam_obs,spec.flam,
                                        spec.flam_sky,spec.sens,wd)

    #Default behaviour is to use the blue exponential extra error if we are in the blue side of a dual spectrograph and the red exponential if we are in the red side of the red exponential. This can be overridden using the force_blue_exp_err or force_red_exp_err.
    fit_blue_exp = False
    fit_red_exp  = False
    if spec.dual_spec:
        if spec.blue:
            fit_blue_exp = True
        else:
            fit_red_exp = True
    try:
        fit_blue_exp = spec.force_blue_exp_err
    except AttributeError:
        pass
    try:
        fit_red_exp = spec.force_red_exp_err
    except AttributeError:
        pass
    #Unless we are consering the extra blue exponential component, for this exercise we need to not consider the Lyalpha forest region. Real IGM absorption appears as noise and throws everything off the board in comparison to the sky template.
    if fit_blue_exp or spec.red:
        kuse = np.argwhere(lam_meanx>0.*u.AA)
    else:
        kuse = np.argwhere(lam_meanx>1300.*u.AA*(1.+spec.zspec))
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
    RON_eff = ((spec.RON**2 + (spec.GAIN/2.)**2)*spec.apsize_pix)**0.5
    Kp = get_error_pars(flam_mean,flam_std,flam_sky_mean,lam_mean,eps_use,RON_eff, fit_blue_exp, fit_red_exp)

    #Get the error spectrum
    flam_err = flam_std_func(Kp,spec.flam,spec.flam_sky,spec.lam_obs,RON_eff,spec.eps())

    #If we are fitting the error with the blue exponential and we are making diagnostic plots, we also want to get the error without that term so that we can compare them.
    if (fit_blue_exp or fit_red_exp) and show_plot:
        kuse2 = np.argwhere(lam_meanx>1300.*u.AA*(1.+spec.zspec))
        lam_mean2 = lam_meanx[kuse2]
        flam_mean2 = flam_meanx[kuse2]
        flam_std2 = flam_stdx[kuse2]
        flam_sky_mean2 = flam_sky_meanx[kuse2]
        sens_use2 = sens_usex[kuse2]

        #Estimate the size of the new wavelength bin.
        dlam_mean2 = np.mean(lam_mean2[1:]-lam_mean2[:-1])

        #Get the eps factor to use.
        eps_use2 = spec.eps(sens=sens_use2,lam_obs=lam_mean2,dlam=dlam_mean2)

        #Fit the std array to the error parameters.
        RON_eff = ((spec.RON**2 + (spec.GAIN/2.)**2)*spec.apsize_pix)**0.5
        Kp_2 = get_error_pars(flam_mean2,flam_std2,flam_sky_mean2,lam_mean2,eps_use2,RON_eff, False, False)

        flam_err_2 = flam_std_func(Kp_2,spec.flam,spec.flam_sky,spec.lam_obs,RON_eff,spec.eps())

    if show_plot:
        fig, ax = plt.subplots()
        flamunit = u.erg/u.s/u.cm**2/u.AA
        plt.plot(lam_mean,flam_std.to(flamunit),'--b',label='First error estimate')
        if fit_blue_exp or fit_red_exp:
            plt.plot(lam_mean2,flam_std2.to(flamunit),'-g', label='$>1300$\AA')
        flam_std_mod = flam_std_func(Kp,flam_mean,flam_sky_mean,lam_mean,RON_eff,eps_use)
        plt.plot(lam_mean,flam_std_mod.to(flamunit),'-k')
        plt.plot(spec.lam_obs,flam_err.to(flamunit),'-r',label='Best-fit')
        if fit_blue_exp or fit_red_exp:
            plt.plot(spec.lam_obs,flam_err_2.to(flamunit),'-m',label='Best-fit no Exp.')
        iduse = spec.name
        if spec.instrument in ["LRIS","DBSP"]:
            if spec.blue:
                iduse += ".blue"
            else:
                iduse += ".red"
        plt.title("{0:s} z={1:.3f}".format(iduse,spec.zspec))
        plt.ylabel("Error Spectrum (erg/s/cm$^2$/\AA)")
        plt.xlabel("Observed Wavelength (\AA)")
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        textbox = ""
        for i in range(len(Kp)):
            if i>0:
                textbox+="\n"
            textbox += "Kp[{0:d}] = {1:.2e}".format(i,Kp[i])
        if fit_blue_exp:
            xloc = 0.80
        else:
            xloc = 0.05
        yloc = 0.99
        plt.text(xloc, yloc, textbox, transform=ax.transAxes, fontsize=8,
            verticalalignment='top', horizontalalignment='left', bbox=props)
        plt.legend(loc='upper center')

        if spec.print_err_plot:
            plt.savefig(iduse+".err.png")
            plt.close()
        else:
            plt.show()

    return flam_err.to(u.erg/(u.s*u.cm**2*u.AA)), Kp
