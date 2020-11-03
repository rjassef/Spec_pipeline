#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.constants import c
from Spec_pipeline.Spec_Reader.iraf_spectrum1d import read_fits_spectrum1d

#B600
spec = read_fits_spectrum1d("gemini_sky_b600_w1049m2834_11may.f.fits")
fnu = spec[0].data * spec[0].unit
lam = spec[0].dispersion
flam = (fnu * c/lam**2).to(u.erg/u.s/u.cm**2/u.AA)
np.savetxt("template_sky_GMOS_B600_1.50arcsec.txt", np.array([lam.value, flam.value]).T)

#R400
spec = read_fits_spectrum1d("gemini_sky_r400_w0730m7218_11may.f.fits")
fnu = spec[0].data * spec[0].unit
lam = spec[0].dispersion
flam = (fnu * c/lam**2).to(u.erg/u.s/u.cm**2/u.AA)
np.savetxt("template_sky_GMOS_R400_1.50arcsec.txt", np.array([lam.value, flam.value]).T)
