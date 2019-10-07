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
# LRIS & DBSP
#####

def proc_LRIS_blue_sky(outname):
    D_t = 10*u.m #Telescope diameter
    proc_LRIS_new("sky_keck_b.w.txt",outname)
    #proc_LRIS_DBSP("sky_keck_b.w.txt",outname,D_t)

def proc_LRIS_red_sky(outname):
    D_t = 10*u.m #Telescope diameter
    proc_LRIS_new("sky_keck_r.w.txt",outname)
    #proc_LRIS_DBSP("sky_keck_r.w.txt",outname,D_t)

def proc_LRIS_new(sky_name, outname):
    
    #Units are fnu in erg/s/cm**2/Hz
    
    sky = np.loadtxt(sky_name)
    lam_sky = sky[:,0]*u.AA
    fnu_sky = sky[:,1]*u.erg/u.s/u.cm**2/u.Hz
    
    flam_sky = (fnu_sky*c/lam_sky**2).to(u.erg/u.s/u.cm**2/u.AA)

    np.savetxt(outname,np.array([lam_sky.value, flam_sky.value]).T,
               fmt='%15.8e %15.8e')

    

def proc_DBSP_blue_sky(outname):
    D_t = 5*u.m #Telescope diameter
    proc_LRIS_DBSP("sky_palomar_b.w.txt",outname,D_t)
    
def proc_DBSP_red_sky(outname):
    D_t = 5*u.m #Telescope diameter
    proc_LRIS_DBSP("sky_palomar_r.w.txt",outname,D_t)
        
def proc_LRIS_DBSP(sky_name,outname,D_t):

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
    A_t = np.pi*(D_t/2.)**2
    texp = 15.*u.minute
    flam_sky = flam_sky * (h*c/lam_sky)/(texp*A_t*bin_size)

    #Write file. Strip the units.
    lam_sky  = lam_sky.to(u.AA).value
    flam_sky = flam_sky.to(u.erg/u.s/u.cm**2/u.AA).value
    flam_sky = np.where(flam_sky>=0,flam_sky,0.)
    np.savetxt(outname,np.array([lam_sky, flam_sky]).T,
               fmt='%15.8e %15.8e')

##

#proc_geminis_sky("template_sky_GMOS.dat")
proc_LRIS_blue_sky("template_sky_LRIS_b.dat")
proc_LRIS_red_sky("template_sky_LRIS_r.dat")
#proc_DBSP_blue_sky("template_sky_DBSP_b.dat")
#proc_DBSP_red_sky("template_sky_DBSP_r.dat")

