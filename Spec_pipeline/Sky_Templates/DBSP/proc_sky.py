#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.constants import h,c

#####
# LRIS & DBSP
#####

def proc_LRIS_new(sky_name, outname):

    #Units are fnu in erg/s/cm**2/Hz

    sky = np.loadtxt(sky_name)
    lam_sky = sky[:,0]*u.AA
    fnu_sky = sky[:,1]*u.erg/u.s/u.cm**2/u.Hz

    flam_sky = (fnu_sky*c/lam_sky**2).to(u.erg/u.s/u.cm**2/u.AA)

    np.savetxt(outname,np.array([lam_sky.value, flam_sky.value]).T,
               fmt='%15.8e %15.8e')


def proc_DBSP_blue_sky(outname,dichroic):
    D_t = 5*u.m #Telescope diameter
    if dichroic=="D55":
        proc_LRIS_DBSP("sky_palomar_{0:s}_b.w.txt".format(dichroic),outname,D_t)
    elif dichroic=="D68":
        proc_LRIS_new("sky_palomar_{0:s}_b.w.txt".format(dichroic),outname)
    else:
        print("Dichroic {0:s} not recognized".format(dichroic))

def proc_DBSP_red_sky(outname,dichroic):
    D_t = 5*u.m #Telescope diameter
    if dichroic=="D55":
        proc_LRIS_DBSP("sky_palomar_{0:s}_r.w.txt".format(dichroic),outname,D_t)
    elif dichroic=="D68":
        proc_LRIS_new("sky_palomar_{0:s}_r.w.txt".format(dichroic),outname)
    else:
        print("Dichroic {0:s} not recognized".format(dichroic))

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

    return

##

proc_DBSP_blue_sky("template_sky_DBSP_D55_b.dat","D55")
proc_DBSP_red_sky("template_sky_DBSP_D55_r.dat","D55")
proc_DBSP_blue_sky("template_sky_DBSP_D68_b.dat","D68")
proc_DBSP_red_sky("template_sky_DBSP_D68_r.dat","D68")
