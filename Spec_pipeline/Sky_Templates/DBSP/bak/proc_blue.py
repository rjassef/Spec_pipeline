#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.constants import c

from Spec_pipeline.Spec_Reader.iraf_spectrum1d import read_fits_spectrum1d as read_spec
from Spec_pipeline.Spec_Reader.rebin_spec import rebin_spec

#We will need to mix the two templates we have to cover the entire wavelength range. We will join them at lam_join with a buffer zone.
buffer = 100.*u.AA
lam_join = 4500.*u.AA

#Read the spectra
sky_b_d55_orig = read_spec("sky_palomar_b600_d55_12oct_b.f.fits")
sky_b_d68_orig = read_spec("sky_palomar_D68_b.w.fits")

#Resample all the spectra to be in 0.5 A binning.
lam_min = np.min(sky_b_d55_orig[0].dispersion)
lam_max = np.max(sky_b_d68_orig[0].dispersion)
wave_new = np.arange(lam_min.to(u.AA).value, lam_max.to(u.AA).value, 0.5) * u.AA

sky_b_d55 = rebin_spec(sky_b_d55_orig[0].dispersion, sky_b_d55_orig[0].data * sky_b_d55_orig[0].unit, wave_new)
sky_b_d68 = rebin_spec(sky_b_d68_orig[0].dispersion, sky_b_d68_orig[0].data * sky_b_d68_orig[0].unit, wave_new)

#Scale the D55 to the D68 in the buffer zone.
kbuff = (wave_new>=lam_join-0.5*buffer) & (wave_new<=lam_join+0.5*buffer)
norm1 = np.average(sky_b_d55[kbuff])
norm2 = np.average(sky_b_d68[kbuff])
norm = (norm2/norm1).to(1.)
sky_b_d55 *= norm.value

#Now, interpolate within the buffer zone to join the spectra.
lam1 = lam_join-0.5*buffer
lam2 = lam_join+0.5*buffer
f = np.zeros(wave_new.shape)
f[wave_new<lam1] = 1.0
f[wave_new>lam2] = 0.0
f[kbuff] = ((wave_new[kbuff]-lam2)/(lam1-lam2)).to(1).value

#Create the final spectrum. Since it is in units of fnu, convert it to flam and save it.
sky_b = sky_b_d55*f + sky_b_d68*(1.-f)
flam_sky_b = (sky_b * c/wave_new**2).to(u.erg/u.s/u.cm**2/u.AA)
np.savetxt("template_sky_DBSP_600-4000_1.50arcsec_b.txt", np.array([wave_new.to(u.AA).value, flam_sky_b.value]).T)

# import matplotlib.pyplot as plt
# plt.plot(wave_new, sky_b_d55, '-b')
# plt.plot(wave_new, sky_b_d68, '-r')
# plt.plot(wave_new, sky_b, '-k')
# plt.show()
