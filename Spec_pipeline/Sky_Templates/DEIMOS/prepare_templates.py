#!/usr/bin/env python

import numpy as np
from astropy.constants import c
import astropy.units as u
import matplotlib.pyplot as plt

from Spec_pipeline.Spec_Reader.iraf_spectrum1d import read_fits_spectrum1d as spec_read

import matplotlib as mpl
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

#Both files are sky templates from the 2013 run. One is just for the blue side, and the other is the whole spectrum. The one for the blue side was produced more recently by Dan and fixes an issue that led to an overestimation of the blue sky levels. So use that for the blue side, and let's start the red side at 6500A.

#Since it from and for that run only in the end, we know that the grating was the 600ZD, the filter the GG400, and the band the slit the 1.5".

#Blue side
temp_b = spec_read("deimos_13dec_sky_b.f.fits")
lam_b  = temp_b[0].dispersion.to(u.AA)
flam_b = temp_b[0].data * temp_b[0].unit * c/lam_b**2
flam_b = flam_b.to(u.erg/u.s/u.cm**2/u.AA)
np.savetxt("template_sky_DEIMOS_600ZD_1.50arcsec_b.txt", np.array([lam_b.value, flam_b.value]).T)

plt.plot(lam_b, flam_b, '-b')
plt.ylabel(r'$F_{{\lambda}}$ ({})'.format(flam_b.unit))
plt.xlabel(r'$\lambda$ ({})'.format(lam_b.unit))
plt.savefig("template_sky_DEIMOS_600ZD_1.50arcsec_b.png")
plt.close()

#Red side
temp_r = np.loadtxt("sky_deimos.dat")
temp_r = temp_r[temp_r[:,0]>=6500.,:]
lam_r  = temp_r[:,0]*u.AA
fnu_r  = temp_r[:,1]*u.erg/u.s/u.cm**2/u.Hz
flam_r = fnu_r * c/lam_r**2
flam_r = flam_r.to(u.erg/u.s/u.cm**2/u.AA)
np.savetxt("template_sky_DEIMOS_600ZD_1.50arcsec_r.txt", np.array([lam_r.value, flam_r.value]).T)

plt.plot(lam_r, flam_r, '-b')
plt.ylabel(r'$F_{{\lambda}}$ ({})'.format(flam_r.unit))
plt.xlabel(r'$\lambda$ ({})'.format(lam_r.unit))
plt.savefig("template_sky_DEIMOS_600ZD_1.50arcsec_r.png")
plt.close()
