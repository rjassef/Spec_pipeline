#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.constants import h,c

#####
# Gemini South
#####

def proc_geminis_sky(outname):

    #Obtained from
    #http://www.gemini.edu/sciops/telescopes-and-sites/observing-condition-constraints/optical-sky-background

    #Read the sky spectrum. Assume a 1" slit and 0.1" pixels.
    sky = np.loadtxt("skybg_50_10.GeminiS.dat")
    lam_sky  = sky[:,0]*u.nm
    flam_sky = sky[:,1] * 1./u.s/u.nm/u.arcsec**2/u.m**2

    #Convert sky spectrum to energy units.
    flam_sky *= (h*c/lam_sky) * (1.*u.arcsec * 0.1*u.arcsec)

    #Write file. Strip the units.
    lam_sky  = lam_sky.to(u.AA).value
    flam_sky = flam_sky.to(u.erg/u.s/u.cm**2/u.AA).value
    np.savetxt(outname,np.array([lam_sky, flam_sky]).T,
               fmt='%15.8e %15.8e')

#####
# LRIS
#####

def proc_LRIS_blue_sky(outname):
    proc_LRIS("sky_keck_b.w.txt",outname)

def proc_LRIS_red_sky(outname):
    proc_LRIS("sky_keck_r.w.txt",outname)

def proc_LRIS(sky_name,outname):

    #Units are not clear, although not fundamental. We'll assume they
    #are counts.
   
    #Read the sky spectrum.
    sky = np.loadtxt(sky_name)
    lam_sky = sky[:,0]*u.AA
    flam_sky = sky[:,1] #Electrons

    #Find the bin size. Assume the lam_sky is the central wavelength
    #of the bin. For most bin except the extremes, will take the bin
    #size as half the distance between the center of the previous bin
    #and the center of the next bin.
    bin_size = np.zeros(lam_sky.shape)
    bin_size[ 0] = (lam_sky[ 1]-lam_sky[ 0]).value
    bin_size[-1] = (lam_sky[-1]-lam_sky[-2]).value
    bin_size[1:-1] = 0.5*( (lam_sky[2: ]-lam_sky[0:-2]).value )
    bin_size *= u.AA

    #We'll assume a 15min exposure time (the standard for our
    #observations).
    D_t = 10.*u.m #Telescope diameter. 
    A_t = np.pi*(D_t/2.)**2
    texp = 15.*u.minute
    flam_sky = flam_sky * (h*c/lam_sky)/(texp*A_t*bin_size)

    #Write file. Strip the units.
    lam_sky  = lam_sky.to(u.AA).value
    flam_sky = flam_sky.to(u.erg/u.s/u.cm**2/u.AA).value
    np.savetxt(outname,np.array([lam_sky, flam_sky]).T,
               fmt='%15.8e %15.8e')

##

proc_geminis_sky("template_sky_GMOS.dat")
proc_LRIS_blue_sky("template_sky_LRIS_b.dat")
proc_LRIS_red_sky("template_sky_LRIS_r.dat")
