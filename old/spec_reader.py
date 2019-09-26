#!/usr/bin/env python 

import numpy as np
from specutils.io.read_fits import read_fits_spectrum1d
from astropy.io import fits
import astropy.units as u
from astropy.constants import h,c
from astropy.io import fits
import re
import matplotlib.pyplot as plt

import fit_error

####

def join_units(a,b):
    aux = np.concatenate((a.value,b.value))
    aux = aux*a.unit
    return aux

####

def read_spec(fnames, z, instrument, line_center=None,spec_side=None):

    #Line that we want to fit.
    if line_center is not None:
        lam_targ = line_center*(1.+z)

    #######################################
    #GMOS - 1 spectrum only, easiest case.#
    #######################################
    if instrument=="GMOS":
        ff = fits.open(fnames[0])
        spec = read_fits_spectrum1d(fnames[0],
                                    dispersion_unit=u.AA, 
                                    flux_unit = u.erg/u.cm**2/u.s/u.Hz)
        lam = spec[0].dispersion
        fnu = spec[0].data * spec[0].unit

        RT = 4.0*u.m #Telescope radius.
        dlam = np.mean(lam[1:]-lam[:-1]) #Mean wavelength bin.
        texp = ff[0].header['EXPTIME']*u.s
        RON = ff[0].header['RDNOISE']
        
        spec_err_name = "error."+fnames[0]


    ###########################################################
    #ASPY - Astropy generated fake QSO spectrum. 1 spec only. #
    ###########################################################
    elif instrument=="ASPY":
        spec = read_fits_spectrum1d(fnames[0],
                                    dispersion_unit=u.AA, 
                                    flux_unit = u.erg/u.cm**2/u.s/u.Hz)
        lam = spec.wavelength
        fnu = spec.flux
        
        RT   = 4.0*u.m #Telescope radius.
        dlam = 1.0*u.AA
        texp = 1000.*u.s
        RON  = 1.0

        spec_err_name = "error."+fnames[0]

    ####################################################################
    #DBSP - Dual spectrum. Use redshift to figure out on which side is #
    #CIV.                                                              #
    ####################################################################
    elif instrument=="DBSP":
        spec_b = read_fits_spectrum1d(fnames[0],
                                      dispersion_unit=u.AA, 
                                      flux_unit = u.erg/u.cm**2/u.s/u.Hz)
        spec_r = read_fits_spectrum1d(fnames[1],
                                      dispersion_unit=u.AA, 
                                      flux_unit = u.erg/u.cm**2/u.s/u.Hz)

        if spec_side is None:
            if lam_targ>spec_b[0].dispersion.min() and \
                    lam_targ<spec_b[0].dispersion.max():
                spec_side = 1
            elif lam_targ>spec_r.dispersion.min() and \
                    lam_targ<spec_r.dispersion.max():
                spec_side = 2
            else:
                spec_side = -1
                

        if spec_side==1:
            lam = spec_b[0].dispersion
            fnu = spec_b[0].data*spec_b[0].unit
            ff = fits.open(fnames[0])
            spec_err_name = "error."+fnames[0]
        elif spec_side==2:
            lam = spec_r.dispersion
            fnu = spec_r.data*spec_r.unit
            ff = fits.open(fnames[1])
            spec_err_name = "error."+fnames[1]
        else:
            lam      = None
            flam     = None
            flam_err = None
            eps      = None
            RON      = None
            print line_center,"not within spectrum"
            return lam, flam, flam_err, eps, RON

        RT   = 2.5*u.m #Telescope radius.
        dlam = np.mean(lam[1:]-lam[:-1]) #Mean wavelength bin.
        texp = float(ff[0].header['EXPTIME'])*u.s
        RON  = float(ff[0].header['RON'])

    ###########################################
    #LRIS - Dual spec, use redshift as in DBSP#
    ###########################################
    elif instrument=='LRIS':
        spec_b = read_fits_spectrum1d(fnames[0],
                                      dispersion_unit=u.AA, 
                                      flux_unit = u.erg/u.cm**2/u.s/u.Hz)
        spec_r = read_fits_spectrum1d(fnames[1],
                                      dispersion_unit=u.AA, 
                                      flux_unit = u.erg/u.cm**2/u.s/u.Hz)

        if spec_side is None:
            if lam_targ>spec_b[0].dispersion.min() and \
                    lam_targ<spec_b[0].dispersion.max():
                spec_side = 1
            elif lam_targ>spec_r[0].dis persion.min() and \
                    lam_targ<spec_r[0].dispersion.max():
                spec_side = 2
            else:
                spec_side = -1
        
        if spec_side==1:
            lam = spec_b[0].dispersion
            fnu = spec_b[0].data*spec_b[0].unit
            ff = fits.open(fnames[0])
            spec_err_name = "error."+fnames[0]

            #https://www2.keck.hawaii.edu/inst/lris/detectors.html
            #Use mean of amplifiers.
            RON = 3.82

            #Load the sensitivity. 
            
            #Load the sky template. 
            sky_file = np.loadtxt("Sky_Templates/sky_keck_b.w.txt")
            sky_lam    = sky_file[:,0]*u.AA
            sky_counts = sky_file[:,1]

        elif spec_side==2:
            lam = spec_r[0].dispersion
            fnu = spec_r[0].data*spec_r[0].unit
            ff = fits.open(fnames[1])
            spec_err_name = "error."+fnames[1]
            #https://www2.keck.hawaii.edu/inst/lris/detectors.html
            #Use mean of amplifiers.
            RON = 4.64
        else:
            lam      = None
            flam     = None
            flam_err = None
            eps      = None
            RON      = None
            print line_center,"not within spectrum"
            return lam, flam, flam_err, eps, RON

        RT   = 5.*u.m #Telescope radius.
        dlam = np.mean(lam[1:]-lam[:-1]) #Mean wavelength bin.
        texp = float(ff[0].header['EXPTIME'])*u.s
        

    #Break down if instrument is not specified.
    else: 
        print instrument,". Instrument not yet configured."
        exit()

    #Failure in case line is not within the spectrum.
    if line_center is not None and \
            (lam_targ<lam.min() or lam_targ>lam.max()):
        lam  = None
        flam = None
        flam_err = None
        eps  = None
        RON  = None
        print line_center,"not within spectrum"
        return lam, flam, flam_err, eps, RON

    flam = (fnu*c/lam**2).to(u.erg/u.cm**2/u.s/u.AA)
    eps = ((np.pi*RT**2) * dlam * texp/(h*c)).to(u.s*u.cm**2/u.erg)

    #Try to open the error spectrum. If not found, generate it. Change
    #the final .fits for .txt.
    spec_err_name = re.sub(".fits",".txt",spec_err_name)
    try:
        cat = open(spec_err_name,"r")
        flam_err = np.loadtxt(cat,usecols=[1])
        flam_err *= u.erg/u.s/u.cm**2/u.AA
    except IOError:
        flam_err, K1, K2 = fit_error.get_error_spec(lam,flam,eps,RON,wd=15)
        np.savetxt(spec_err_name,np.array([lam,flam_err]).T)

    return lam, flam, flam_err, eps, RON
