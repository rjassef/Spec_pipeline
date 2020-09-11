#!/usr/bin/env python

import numpy as np
import astropy.units as u
import re
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel, convolve

from Spec_pipeline.Spec_Reader.iraf_spectrum1d import read_fits_spectrum1d

#Array that holds all information about the gratings. Resolutions in angstroms/mm.
res_page = dict()
res_page['158-7560']  = [135.,201.]
res_page['300-3990']  = [140.,'-' ]
res_page['316-7500']  = ['-' ,102.]
res_page['600-4000']  =	[ 71., '-']
res_page['600-10000'] =	['-' , 54.]
#res_page['1200-4700'] =	[ 36., '-']
#res_page['1200-7100'] =	[ 36., 27.]
#res_page['1200-9400'] =	[ 35., 26.]

#Set the scale of arcseconds to microns. This should be independent of the detector, but we can easilty get it from the detector parameters combining the pixel size in microns with the plate_scale (arcsec/pixel).
scale = 0.389*u.arcsec / (15.*u.micron)

#Read the blue template.
spec = read_fits_spectrum1d("sky_palomar_D68_b.w.fits")

#Get the aperture of the template
slit_width_temp = float(spec[0].header['APERTURE']) * u.arcsec

#Get the grating of the template.
grating_temp = spec[0].header['GRATING']
grating_temp  = re.sub("7150","7500",grating_temp)
grating_temp  = re.sub("/","-",grating_temp)
if grating_temp == "3167500":
    grating_temp = "316-7500"

#Get the template grating resolution
res_temp = res_page[grating_temp][0] * u.angstrom/u.mm

#Set the template resolution
FWHM_res_temp  = (slit_width_temp/scale) * res_temp
sigma_res_temp = FWHM_res_temp/(2.*(2.*np.log(2.))**0.5)

#Now, go through every possible slit width and grating to construct the new templates. If the resolution is better than that of the template, display a warning.
slit_widths = np.array([1.5, 2.0])*u.arcsec
for slit_width in slit_widths:
    for grating in res_page.keys():

        if res_page[grating][0]=='-':
            continue

        #Get the template grating resolution
        res = res_page[grating][0] * u.angstrom/u.mm

        #Set the resolution to use
        FWHM_res  = (slit_width/scale) * res
        sigma_res = FWHM_res/(2.*(2.*np.log(2.))**0.5)

        #Compare to that of the template. If worse, convolve. If better, warn the user that nothing is being done.
        if sigma_res>sigma_res_temp:
            print("Creating template for {0:s} grating, {1:s} slit width".format(grating, slit_width))
            sigma_use = (sigma_res**2-sigma_res_temp**2)**0.5
            bin_size = np.mean(spec[0].dispersion[1:]-spec[0].dispersion[:-1])
            sigma_use_pix = (sigma_use/bin_size).to(1.).value
            gauss_kernel = Gaussian1DKernel(stddev=sigma_use_pix)
            sky_conv = convolve(spec[0].data, gauss_kernel)
        else:
            print("Warning. Template resolution is worse than requested resolution. Using unmodified template.")
            print("Template grating: {0:s}  Template aperture: {1:f}".format(grating_temp, slit_width_temp))
            print("Requested grating: {0:s} Requested aperture: {1:f}".format(grating, slit_width))
            sky_conv = spec[0].data

        #Write the resulting sky model.
        fname = "template_sky_DBSP_{0:s}_{1:.2f}arcsec_b.txt".format(grating, slit_width.to(u.arcsec).value)
        np.savetxt(fname,np.array([spec[0].dispersion, sky_conv]).T)
        plt.plot(spec[0].dispersion, sky_conv,'-b')
        plt.savefig(re.sub(".txt",".png",fname))
        plt.close()
