#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.constants import h,c


#####
# LRIS & DBSP
#####

def proc_LRIS_blue_sky(outname):
    D_t = 10*u.m #Telescope diameter
    proc_LRIS_new("sky_keck_b.w.txt",outname)

def proc_LRIS_red_sky(outname):
    D_t = 10*u.m #Telescope diameter
    proc_LRIS_new("sky_keck_r.w.txt",outname)

def proc_LRIS_new(sky_name, outname):

    #Units are fnu in erg/s/cm**2/Hz

    sky = np.loadtxt(sky_name)
    lam_sky = sky[:,0]*u.AA
    fnu_sky = sky[:,1]*u.erg/u.s/u.cm**2/u.Hz

    flam_sky = (fnu_sky*c/lam_sky**2).to(u.erg/u.s/u.cm**2/u.AA)

    np.savetxt(outname,np.array([lam_sky.value, flam_sky.value]).T,
               fmt='%15.8e %15.8e')

##

proc_LRIS_blue_sky("template_sky_LRIS_b.dat")
proc_LRIS_red_sky("template_sky_LRIS_r.dat")
