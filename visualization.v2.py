#!/usr/bin/env python 

import numpy as np
import fit_chi2 as fit
import spec_reader
import astropy.units as u
from astropy.constants import c
import matplotlib.pyplot as plt
from spec_reader import read_spec
from scipy.special import betainc

#####

show_plots = True
#show_plots = False

#####

def fit_lines(fnames,z,instrument,linecat):
    
    cat = open(linecat)
    output = []
    for line in cat:
        if line[0]=='#':
            continue
        x = line.split()
        x[1:] = [float(ix) for ix in x[1:]]
        
        lam_cen_cat = x[1]*u.AA
        line_velocity_region = x[2]*u.km/u.s
        n_cont = (len(x)-3)/2
        continuum_regions = []
        for i in range(n_cont):
            continuum_regions.append([x[2*i+3],x[2*i+4]])
        continuum_regions *= u.AA
            
        lam_cen, flam_line_cen, FWHM_v, F, p, lam_mod, flam_mod = line_fitter(
            fnames,z,instrument,lam_cen_cat,
            line_velocity_region,continuum_regions)
            
        if lam_cen is None:
            output_line = "{0:6s} ---".format(x[0])
        else:
            output_line = "{0:6s} {1:6.1f} {2:6.1f} {3:.3e}".format(
                x[0],lam_cen,FWHM_v,1.-p)
            if 1.-p>0.8 and show_plots:
                plt.plot(lam_mod,flam_mod,'-b')

        print output_line
        output.append(output_line+"\n")
        
    print
    return output

def line_fitter(fnames,z,instrument,lam_cen_cat,
                line_velocity_region,continuum_regions):
    
    #Read the spectrum.
    lam_obs, flam, flam_err, eps, RON = spec_reader.read_spec(fnames, z, 
                                                              instrument,
                                                              lam_cen_cat)
    if lam_obs is None:
        return None,None,None,None,None,None,None
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
    if F<0.:
        p = 1.0
    else:
        nnu1 = float(nu_no_line-nu)
        nnu2 = float(nu)
        w = nnu1*F/(nnu1*F+nnu2)
        p = 1.-betainc(nnu1/2.,nnu2/2.,w)

    lamx = lam_rest[np.abs(dv)<5.*sigma_v]
    flamx = flam_mod[np.abs(dv)<5.*sigma_v]
    return lam_cen, flam_line_cen, FWHM_v, F, p, lamx, flamx

###

cato = open("visualization.summary.txt","w")

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
    cato.write(obj_id+"\n")

    #Plot the entire spectrum:
    if show_plots:
        for i in range(len(fnames)):
            lam, flam, flam_err, eps, RON = read_spec(
                fnames, z, instrument, spec_side=i+1)
            plt.plot(lam/(1.+z),flam,linestyle='solid',color='xkcd:grey')

    #Fit the emission lines
    output = fit_lines(fnames,z,instrument,"lines.txt")

    if show_plots:
        plt.show()

    for line in output:
        cato.write(line)
    cato.write("\n")

