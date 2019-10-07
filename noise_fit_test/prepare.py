#!/usr/bin/env python 

import numpy as np
import astropy.units as u
from astropy.constants import c,h
from astropy.table import Table
import matplotlib.pyplot as plt

from Spec_pipeline.Spec_Reader.rebin_spec import rebin_spec

###

flamunit = u.erg/u.s/u.cm**2/u.AA
waveunit = u.AA

RON_b = 3.82
RON_r = 4.64

###

#This scripts simulates the observations of the vandenberk composite
#template with Keck. 

def create_fake_LRIS_QSO_spec(lam,flam,flam_sky,sens,
                              Vmag = 24.0,
                              RT   =  5.*u.m,
                              texp = 15.*u.minute,
                              RON  = None):

    #Spectrum normalization
    lam_eff = 5500.*u.AA
    fnu_V   = 3631.*10.**(-0.4*Vmag)*u.Jy
    flam_V  = (fnu_V*c/lam_eff**2).to(flamunit)

    k = np.argmin(np.abs(lam-lam_eff))
    norm = flam_V/flam[k]
    flam *= norm
    flam = flam.to(flamunit)

    #Transform to counts.
    dlam = np.mean(lam[1:]-lam[:-1])
    flam_to_counts = sens * texp * (np.pi*RT)**2 * dlam / (h*c/lam)
    counts     = flam     * flam_to_counts
    counts_sky = flam_sky * flam_to_counts

    #Get the error spectrum, and then resample the original spectrum.
    flam_err = (counts + counts_sky + RON**2)**0.5 / flam_to_counts
    flam_err = flam_err.to(flamunit)
    flam_new = np.random.normal(flam,flam_err)

    return flam_new, flam_err

#####

#Read the spectrum. Note that, while the table lacks the units,
#Vanden Berk et al. (2001) shows that the flux density is flambda
#in arbitrary units.
qso = Table.read("vandenberk_composite.txt",format='ascii.cds')
lam_rest  = qso['Wave'].to(waveunit)
flam      = qso['FluxD'] * flamunit

z = 2.6
lam = lam_rest*(1.+z)

#Separate in blue and red components.
lam_b  = lam[ (lam>=3200.*u.AA) & (lam<= 5500.*u.AA)]
flam_b = flam[(lam>=3200.*u.AA) & (lam<= 5500.*u.AA)]

lam_r  = lam[ (lam>=5500.*u.AA) & (lam<=10000.*u.AA)]
flam_r = flam[(lam>=5500.*u.AA) & (lam<=10000.*u.AA)]

#Load the Sensitivities and rebin them.
sens_b_temp = np.loadtxt(
    "../Spec_pipeline/Sensitivity_Files/Sens_LRIS_B600.txt")
lam_sens_b_temp = sens_b_temp[:,0]*waveunit
sens_b_temp     = sens_b_temp[:,1]*u.dimensionless_unscaled
sens_b = rebin_spec(lam_sens_b_temp, sens_b_temp, lam_b)

sens_r_temp = np.loadtxt(
    "../Spec_pipeline/Sensitivity_Files/Sens_LRIS_R400.txt")
lam_sens_r_temp = sens_r_temp[:,0]*waveunit
sens_r_temp     = sens_r_temp[:,1]*u.dimensionless_unscaled
sens_r = rebin_spec(lam_sens_r_temp, sens_r_temp, lam_r)

#Load the sky templates and rebin them.
sky_b_temp  = np.loadtxt(
    "../Spec_pipeline/Sky_Templates/template_sky_LRIS_b.dat")
lam_sky_b_temp = sky_b_temp[:,0]*waveunit
sky_b_temp     = sky_b_temp[:,1]*flamunit
sky_b = rebin_spec(lam_sky_b_temp, sky_b_temp, lam_b)

sky_r_temp  = np.loadtxt(
    "../Spec_pipeline/Sky_Templates/template_sky_LRIS_r.dat")
lam_sky_r_temp = sky_r_temp[:,0]*waveunit
sky_r_temp     = sky_r_temp[:,1]*flamunit
sky_r = rebin_spec(lam_sky_r_temp, sky_r_temp, lam_r)


flam_new_b, flam_err_b = create_fake_LRIS_QSO_spec(lam_b,flam_b,sky_b,sens_b,
                                                   RON=RON_b)
flam_new_r, flam_err_r = create_fake_LRIS_QSO_spec(lam_r,flam_r,sky_r,sens_r,
                                                   RON=RON_r)

plt.plot(lam_b,flam_new_b,linestyle='solid',color='xkcd:blue')
plt.plot(lam_b,flam_b,linestyle='solid',color='xkcd:grey')

plt.plot(lam_r,flam_new_r,linestyle='solid',color='xkcd:red')
plt.plot(lam_r,flam_r,linestyle='solid',color='xkcd:grey')
plt.show()
