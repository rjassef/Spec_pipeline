#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.constants import c
import matplotlib.pyplot as plt

from Spec_pipeline.Spec_Reader.iraf_spectrum1d import read_fits_spectrum1d

import matplotlib as mpl
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

#B600
spec = read_fits_spectrum1d("gmosS_12dec_sky_b.f.fits")
fnu = spec[0].data * spec[0].unit
lam = spec[0].dispersion
flam = (fnu * c/lam**2).to(u.erg/u.s/u.cm**2/u.AA)
np.savetxt("template_sky_GMOS_B600_1.50arcsec.txt", np.array([lam.value, flam.value]).T)

plt.plot(lam, flam, '-b')
plt.ylabel(r'$F_{{\lambda}}$ ({})'.format(flam.unit))
plt.xlabel(r'$\lambda$ ({})'.format(lam.unit))
plt.savefig("template_sky_GMOS_B600_1.50arcsec.png")
plt.close()


#R400
spec = read_fits_spectrum1d("gemini_r400_sky.f.fits")
fnu = spec[0].data * spec[0].unit
lam = spec[0].dispersion
flam = (fnu * c/lam**2).to(u.erg/u.s/u.cm**2/u.AA)
np.savetxt("template_sky_GMOS_R400_1.50arcsec.txt", np.array([lam.value, flam.value]).T)

plt.plot(lam, flam, '-b')
plt.ylabel(r'$F_{{\lambda}}$ ({})'.format(flam.unit))
plt.xlabel(r'$\lambda$ ({})'.format(lam.unit))
plt.savefig("template_sky_GMOS_R400_1.50arcsec.png")
plt.close()
