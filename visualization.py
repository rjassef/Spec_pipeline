#!/usr/bin/env python 

import numpy as np
import fit_chi2 as fit
import spec_reader
import astropy.units as u
from astropy.constants import c
import matplotlib.pyplot as plt
from spec_reader import read_spec
from scipy.special import betainc

def fit_civ(fnames,z,instrument):
    
    #Define the fitting regions for the emission line.
    continuum_regions = [[1425., 1470.],[1680.,1705.]]*u.AA
    line_velocity_region = 10000.*u.km/u.s

    #Read the spectrum.
    lam_cen_cat = 1549.*u.AA
    lam_obs, flam, flam_err, eps, RON = spec_reader.read_spec(fnames, z, 
                                                              instrument,
                                                              lam_cen_cat)
    lam_rest = lam_obs/(1.+z)

    #Initial guesses:
    lam_cen_0 = lam_cen_cat
    sigma_v_0  = 1000.*u.km/u.s

    #Fit the emission line.
    lam_cen, flam_line_cen, sigma_v, a, b, flam_mod = fit.fit(
        lam_rest, flam, flam_err, continuum_regions, 
        line_velocity_region, lam_cen_0, sigma_v_0)
    FWHM_v  = sigma_v*2.*(2.*np.log(2.))**0.5

    chi2 = fit.chi2_fit([lam_cen.value, flam_line_cen.value, sigma_v.value],
                        flam, flam_err, lam_rest, a, b, line_velocity_region,
                        lam_cen, flam_line_cen)
    chi2_no_line = fit.chi2_fit([lam_cen.value, 0., 
                                 sigma_v.value],
                                flam, flam_err, lam_rest, a, b, 
                                line_velocity_region, lam_cen, 0.)

    #Get the degrees of freedom.
    dv = (c*(lam_rest/lam_cen-1.)).to(u.km/u.s)
    dvabs = np.abs(dv)
    n_datapoints = len(flam[dvabs<line_velocity_region])
    nu = n_datapoints - 2 - 3
    nu_no_line = n_datapoints
    chi2_nu = chi2/float(nu)
    chi2_no_line_nu = chi2_no_line/float(nu_no_line)

    F = ((chi2_no_line-chi2)/float(nu_no_line-nu)) / (chi2/float(nu))
    F = F.value
    nnu1 = float(nu_no_line-nu)
    nnu2 = float(nu)
    w = nnu1*F/(nnu1*F+nnu2)
    print nnu1,nnu2,F,w
    p = 1.-betainc(nnu1/2.,nnu2/2.,w)

    print "CIV: ",lam_cen,FWHM_v,chi2_nu,chi2_no_line_nu,F,p
    lamx = lam_rest[np.abs(dv)<5.*sigma_v]
    flamx = flam_mod[np.abs(dv)<5.*sigma_v]
    return lamx, flamx

###

cat = open("data.txt")
for line in cat:
    
    if line[0]=="#":
        continue

    x = line.split()
    obj_id = x[0]
    z = float(x[1])
    instrument = x[2]
    fnames = x[3:]
    print obj_id

    #Plot the entire spectrum:
    for i in range(len(fnames)):
        lam, flam, flam_err, eps, RON = read_spec(
            fnames, z, instrument, spec_side=i+1)
        plt.plot(lam/(1.+z),flam,linestyle='solid',color='xkcd:grey')

    #Fit and plot the CIV line.
    lam_civ, flam_civ = fit_civ(fnames,z,instrument)
    plt.plot(lam_civ,flam_civ,'-b')

    plt.show()
