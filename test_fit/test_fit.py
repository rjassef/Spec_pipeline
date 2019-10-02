#!/usr/bin/env python

import warnings
warnings.simplefilter("ignore")

import numpy as np
import astropy.units as u
from astropy.constants import h,c
import matplotlib.pyplot as plt
import subprocess

from Spec_pipeline.Spec_Reader.Spec import Spec
from Spec_pipeline.Spec_Reader.rebin_spec import rebin_spec
from Spec_pipeline.Line_Fitter.default_lines import Default_Line_fit

#Start by creating a fake spectrum.
z = 3.0
fake = Spec("Fake",z)
fake.spec_err_name = "fake_err.txt"
subprocess.call("rm data/fake_err.txt",shell=True)
fake.sens = 1.0

#Set the wavelength range.
lam_rest_min = 1200.
lam_rest_max = 1800.
dlam_rest = 2.
fake.lam_obs = np.arange(lam_rest_min,lam_rest_max,dlam_rest)*(1+z)*u.AA

#First, add the CIV emission line.
civ_line = Default_Line_fit("CIV")
civ_line.lam_cen_fit = civ_line.line_center
civ_line.flam_line_cen_fit = 1e-17*u.erg/u.s/u.cm**2/u.AA
civ_line.sigma_v_fit = 1500.*u.km/u.s
flam_model = civ_line.flam_line_model(fake.lam_rest)

#The continuum will simply be a power-law continuum.
flam_norm = 8.e-18*u.erg/u.s/u.cm**2/u.AA
lam_norm = civ_line.line_center
flam_model += flam_norm*(fake.lam_rest/lam_norm)**-1.0

#Load the sky. We'll use the GMOS sky.
sky_temp = np.loadtxt(\
                      "../Spec_pipeline/Sky_Templates/template_sky_GMOS.dat")
lam_sky = sky_temp[:,0]*u.AA
flam_sky_orig = sky_temp[:,1]*u.erg/u.s/u.cm**2/u.AA

#Rebin the sky template to the object spectrum.
fake.flam_sky = rebin_spec(lam_sky, flam_sky_orig, fake.lam_obs)

#Now, simulate observations with a given telescope. We'll assume
#efficiency of 1 for now.
fake.RT = 1.*u.m
fake.texp = 1.*u.minute
fake.RON  = 3.0
fake.dlam = dlam_rest*(1.+z)*u.AA
fact = fake.texp*(np.pi*fake.RT**2)*fake.dlam/(h*c/fake.lam_obs)
N_Gamma_Source = flam_model*fact
N_Gamma_Sky = fake.flam_sky*fact
N_Gamma_err_model = (N_Gamma_Source + N_Gamma_Sky + fake.RON**2)**0.5
flam_err_model = (N_Gamma_err_model/fact).to(u.erg/u.s/u.cm**2/u.AA)

fake.flam = np.random.normal(flam_model,flam_err_model)*flam_model.unit

#Now, run the fit.
civ_fit = Default_Line_fit("CIV")
civ_fit.run_fit(fake)
print(civ_fit.FWHM_v,civ_line.FWHM_v)

plt.plot(fake.lam_rest,fake.flam,'-b')
plt.plot(fake.lam_rest,flam_model,'-k')
plt.plot(fake.lam_rest,civ_fit.flam_model(fake.lam_rest),'-r')
plt.show()

#civ_fit.run_MC(fake,1000)
#print(civ_fit.FWHM_v_low, civ_fit.FWHM_v_hig)

