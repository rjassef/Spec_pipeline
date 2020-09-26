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

####

proc_geminis_sky("template_sky_GMOS.dat")
