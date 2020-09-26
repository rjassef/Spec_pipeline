#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.constants import h,c

def proc_DEIMOS_sky(outname):

    #Read the sky spectrum.
    sky = np.loadtxt("sky_deimos.dat")
    lam_sky = sky[:,0]*u.AA
    fnu_sky = sky[:,1]*u.erg/u.s/u.cm**2/u.Hz

    #Calculate flam.
    flam_sky = fnu_sky*c/lam_sky**2

    #Write file. Strip the units.
    lam_sky  = lam_sky.to(u.AA).value
    flam_sky = flam_sky.to(u.erg/u.s/u.cm**2/u.AA).value
    np.savetxt(outname,np.array([lam_sky, flam_sky]).T,
               fmt='%15.8e %15.8e')

    return

##

proc_DEIMOS_sky("template_sky_DEIMOS.dat")
