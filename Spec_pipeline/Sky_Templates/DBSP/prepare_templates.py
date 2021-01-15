#!/usr/bin/env python

import numpy as np
import astropy.units as u
import re
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel, convolve
import subprocess
from copy import deepcopy

from Spec_pipeline import DBSP_Spec
from Spec_pipeline.Spec_Reader.rebin_spec import rebin_spec

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


def combine(spec1, spec2, lam_cross, dlam_cross):

    #Get the minimum and maximum wavelengths and the mean wavelength bin size.
    lam_min = np.min(np.concatenate([spec1.lam_obs, spec2.lam_obs]))
    lam_max = np.max(np.concatenate([spec1.lam_obs, spec2.lam_obs]))
    dlam = np.mean(np.concatenate([spec1.lam_obs[1:]-spec1.lam_obs[:-1], spec2.lam_obs[1:]-spec2.lam_obs[:-1]]))

    #Set the final wavelength grid and rebin the current spectra.
    lam_interp = np.arange(lam_min.to(u.AA).value, lam_max.to(u.AA).value, dlam.to(u.AA).value)*u.AA
    flam1 = rebin_spec(spec1.lam_obs, spec1.flam, lam_interp)
    flam2 = rebin_spec(spec2.lam_obs, spec2.flam, lam_interp)

    #Renormalize flam2 to match flam1 in the crossover region.
    kcross = (lam_interp>lam_cross-dlam_cross) & (lam_interp<lam_cross+dlam_cross)
    norm = np.sum(flam1[kcross]*flam2[kcross])/np.sum(flam2[kcross]**2)
    flam2 *= norm

    #Create the final spectrum.
    flam = np.zeros(len(flam1))*flam1.unit
    flam[lam_interp<=lam_cross-dlam_cross] = flam1[lam_interp<=lam_cross - dlam_cross]
    flam[lam_interp>=lam_cross+dlam_cross] = flam2[lam_interp>=lam_cross + dlam_cross]
    n = len(kcross[kcross==True])
    frac = np.arange(0,n,1)/float(n-1)
    flam[kcross] = (1.0-frac)*flam1[kcross] + frac*flam2[kcross]

    spec = deepcopy(spec1)
    spec.lam_obs = lam_interp
    spec.flam = flam

    return spec

#####

#Create a fake local sky file that has a very long wavelength range, from 2000A to 1.5um.
fake_local_sky_lam = np.arange(2000., 15000., 100.)
fake_local_sky_flam = np.ones(fake_local_sky_lam.shape)
np.savetxt("fake_local_sky.txt", np.array([fake_local_sky_lam, fake_local_sky_flam]).T)

#Go through every spectrum in the data folder, and load them with the right properties. If more than one spectrum exists with the same grating+channel+slit_width, then load the one with the longest wavelength range.
aux = open("aux.dat","w")
subprocess.call("ls data",shell=True,stdout=aux)
aux.close()
fnames = np.genfromtxt("aux.dat", dtype="U")
template = dict()
for fname in fnames:

    print()
    print(fname)

    #Read the spectrum.
    blue = False
    red = False
    if re.search("_b\.", fname):
        blue = True
    else:
        red = True
    spec = DBSP_Spec("Sky_template", zspec=0., fits_file=fname, blue=blue, red=red, local_sky_file="fake_local_sky.txt")

    #If the spectrum does not have units, then we need to apply the sensitivity calibration.
    if spec.flam.unit == u.dimensionless_unscaled:
        spec.flam = spec.flam/spec.eps()

    temp_name = "{0:s}_{1:.02f}arcsec_{2:s}".format(spec.grating, spec.slit_width.to(u.arcsec).value, spec.channel)

    if temp_name in template.keys():
        if spec.dichroic not in template[temp_name].keys():
            template[temp_name][spec.dichroic] = spec
        else:
            #Select the one with the longest baseline.
            temp = template[temp_name][spec.dichroic]
            dlam_temp = np.max(temp.flam)-np.min(temp.flam)
            dlam = np.max(spec.flam)-np.min(spec.flam)
            if dlam > dlam_temp:
                template[temp_name][spec.dichroic] = spec
    else:
        template[temp_name] = dict()
        template[temp_name][spec.dichroic] = spec

#Now, if for any template we have both dichroics, D55 and D68, then join them by interpolation in the 4900 +/- 100 AA range.
combined_template = dict()
for temp_name in template.keys():
    #If only one, just use that one.
    dichs = list(template[temp_name].keys())
    if len(dichs)==1:
        combined_template[temp_name] = template[temp_name][dichs[0]]

    #Otherwise, interpolate them.
    elif len(dichs)==2:
        spec1 = template[temp_name]['D55']
        spec2 = template[temp_name]['D68']
        combined_template[temp_name] = combine(spec1, spec2, 4900*u.AA, 100.*u.AA)

    else:
        print("Error. Too many dichroics for temp_name: ",temp_name)
        exit()

#Now, check that for each template that we have a 1.5" slit, we also have a 2.0" slit one. If not, or if the 2.0" one has a shorter wavelength range, then we need to build it by convolving the 1.5" spectrum with a Gaussian.
temp_names = list(combined_template.keys())
for temp_name in temp_names:
    if re.search("1.50arcsec", temp_name):
        temp_name2 = re.sub("1.50arcsec", "2.00arcsec", temp_name)
        if temp_name2 in combined_template.keys():
            lam1_range = np.max(combined_template[temp_name].lam_obs.to(u.AA).value) - np.min(combined_template[temp_name].lam_obs.to(u.AA).value)
            lam2_range = np.max(combined_template[temp_name2].lam_obs.to(u.AA).value) - np.min(combined_template[temp_name2].lam_obs.to(u.AA).value)
            if lam2_range>=lam1_range:
                continue

        print("Decreasing resolution",temp_name,"to get",temp_name2)
        spec = decrease_resolution(combined_template[temp_name], 2.*u.arcsec)
        combined_template[temp_name2] = spec

#Save the templates.
for temp_name in combined_template.keys():
    lam = combined_template[temp_name].lam_obs.to(u.AA).value
    flam = combined_template[temp_name].flam.to(u.erg/u.cm**2/u.s/u.AA).value
    np.savetxt("template_sky_DBSP_"+temp_name+".txt", np.array([lam,flam]).T)

    plt.plot(lam, flam, 'b-')
    plt.savefig("template_sky_DBSP_"+temp_name+".png")
    plt.close()

#Erase the fake sky template.
subprocess.call("rm fake_local_sky.txt aux.dat", shell=True)
