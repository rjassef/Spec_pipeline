#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.constants import h,c

#####
# ESO
#####

def proc_eso_sky(deg,outname):

    #See README file for details.

    #Read the sky spectrum. Assume a 1" slit and 0.1" pixels.
    sky = np.loadtxt("Night_Sky_{0:d}.dat".format(deg))
    lam_sky  = sky[:,0]*u.nm
    flam_sky = sky[:,1] * 1./u.s/u.m**2/u.micron/u.arcsec**2

    #Convert sky spectrum to energy units.
    flam_sky *= (h*c/lam_sky) * (1.*u.arcsec * 0.1*u.arcsec)

    #Write file. Strip the units.
    lam_sky  = lam_sky.to(u.AA).value
    flam_sky = flam_sky.to(u.erg/u.s/u.cm**2/u.AA).value
    np.savetxt(outname,np.array([lam_sky, flam_sky]).T,
               fmt='%15.8e %15.8e')

##

proc_eso_sky(  0,"template_sky_ESO_0.dat")
proc_eso_sky( 90,"template_sky_ESO_90.dat")
proc_eso_sky(180,"template_sky_ESO_180.dat")


