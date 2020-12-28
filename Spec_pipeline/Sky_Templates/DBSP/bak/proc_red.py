#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.constants import c, h

from Spec_pipeline.Spec_Reader.iraf_spectrum1d import read_fits_spectrum1d as read_spec

#Much easier than the blue side, this is just read and convert. Note, however, that the flux is in electrons.
sky_r = read_spec("sky_palomar_r.w.fits")
lam_sky = sky_r[0].dispersion
flam_sky = sky_r[0].data

#Find the bin size. Assume the lam_sky is the central wavelength
#of the bin. For most bin except the extremes, will take the bin
#size as half the distance between the center of the previous bin
#and the center of the next bin.
bin_size = np.zeros(lam_sky.shape)*u.AA
bin_size[ 0] = lam_sky[ 1]-lam_sky[ 0]
bin_size[-1] = lam_sky[-1]-lam_sky[-2]
bin_size[1:-1] = 0.5*(lam_sky[2: ]-lam_sky[0:-2])

#We'll assume a 15min exposure time (the standard for our observations).
D_t = 5.0*u.m #Telescope aperture.
A_t = np.pi*(D_t/2.)**2
texp = 15.*u.minute
flam_sky = flam_sky * (h*c/lam_sky)/(texp*A_t*bin_size)

#In the blue edge, the template has some negative bins. Turn them to 0.
flam_sky = np.where(flam_sky>0., flam_sky, 0.)

#Write file. Strip the units.
flam_sky_r = flam_sky.to(u.erg/u.s/u.cm**2/u.AA)
np.savetxt("template_sky_DBSP_316-7500_1.50arcsec_r.txt", np.array([lam_sky.to(u.AA).value, flam_sky_r.value]).T)

# import matplotlib.pyplot as plt
# plt.plot(wave_new, sky_b_d55, '-b')
# plt.plot(wave_new, sky_b_d68, '-r')
# plt.plot(wave_new, sky_b, '-k')
# plt.show()
