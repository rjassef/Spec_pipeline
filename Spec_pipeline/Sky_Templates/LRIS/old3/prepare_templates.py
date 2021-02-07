#!/usr/bin/env python

import numpy as np
import astropy.units as u
import re
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel, convolve
import subprocess
from copy import deepcopy

from Spec_pipeline import LRIS_Spec

#####

def decrease_resolution(spec, slit_width_targ):

    if spec.slit_width>slit_width_targ:
        print("Input spectrum has worse resolution.")
        return None

    sigma_use = spec.sigma_res * ((slit_width_targ/spec.slit_width)**2-1.)**0.5
    bin_size = np.mean(spec.lam_obs[1:]-spec.lam_obs[:-1])
    sigma_use_pix = (sigma_use/bin_size).to(1.).value
    gauss_kernel = Gaussian1DKernel(stddev=sigma_use_pix)
    flam_conv = convolve(spec.flam, gauss_kernel)

    spec_conv = deepcopy(spec)
    spec_conv.flam = flam_conv*spec.flam.unit
    spec_conv.slit_width = slit_width_targ

    return spec_conv


#####

#Create a fake local sky file that has a very long wavelength range, from 2000A to 1.5um.
fake_local_sky_lam = np.arange(2000., 15000., 100.)
fake_local_sky_flam = np.ones(fake_local_sky_lam.shape)
np.savetxt("fake_local_sky.txt", np.array([fake_local_sky_lam, fake_local_sky_flam]).T)

#Go through every spectrum in the data folder, and load them with the right properties. If more than one spectrum exists with the same grating+channel+slit_width, then load the one with the longest wavelength range.
aux = open("aux.dat","w")
subprocess.call("ls data/*.fits",shell=True,stdout=aux)
aux.close()
fnames = np.genfromtxt("aux.dat", dtype="U")
template = dict()
for fname in fnames:

    fname = re.sub("data/","",fname)
    print()
    print(fname)

    #Read the spectrum.
    blue = False
    red = False
    if re.search("_b\.", fname):
        blue = True
    else:
        red = True
    spec = LRIS_Spec("Sky_template", zspec=0., fits_file=fname, blue=blue, red=red, local_sky_file="fake_local_sky.txt")

    #If the spectrum does not have units, then we need to apply the sensitivity calibration.
    if spec.flam.unit == u.dimensionless_unscaled:
        spec.flam = spec.flam/spec.eps()

    temp_name = "{0:s}_{1:.02f}arcsec_{2:s}".format(spec.grating, spec.slit_width.to(u.arcsec).value, spec.channel)
    if temp_name in template.keys():
        #Select the one with the longest baseline.
        dlam_temp = np.max(template[temp_name].flam)-np.min(template[temp_name].flam)
        dlam = np.max(spec.flam)-np.min(spec.flam)
        if dlam > dlam_temp:
            template[temp_name] = spec
    else:
        template[temp_name] = spec

#Now, check that for each template that we have a 1.5" slit. If not, the we need to build it by convolving the 1.0" spectrum with a Gaussian.
temp_names = list(template.keys())
for temp_name in temp_names:
    if re.search("1.00arcsec", temp_name):
        temp_name2 = re.sub("1.00arcsec", "1.50arcsec", temp_name)
        if temp_name2 in template.keys():
            continue
        else:
            print("Decreasing resolution",temp_name,"to get",temp_name2)
            spec = decrease_resolution(template[temp_name], 1.5*u.arcsec)
            template[temp_name2] = spec

#Save the templates.
for temp_name in template.keys():
    lam = template[temp_name].lam_obs.to(u.AA).value
    flam = template[temp_name].flam.to(u.erg/u.cm**2/u.s/u.AA).value
    np.savetxt("template_sky_LRIS_"+temp_name+".txt", np.array([lam,flam]).T)

    plt.plot(lam, flam, 'b-')
    plt.savefig("template_sky_LRIS_"+temp_name+".png")
    plt.close()

#Erase the fake sky template.
subprocess.call("rm fake_local_sky.txt aux.dat", shell=True)
