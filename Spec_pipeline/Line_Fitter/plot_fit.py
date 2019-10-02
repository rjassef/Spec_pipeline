#!/usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.constants import c
from scipy.signal import savgol_filter

from .MC_errors_general import get_error

###

def plot_line(lam,label):
    xmin, xmax = plt.xlim()
    if lam<xmin or lam>xmax:
        return
    ymin, ymax = plt.ylim()
    plt.plot([lam,lam],[ymin,0.9*ymax],linestyle='dashed',color='xkcd:sky blue')
    plt.text(lam,0.95*ymax,label,horizontalalignment='center',fontsize=8)
    return


###

#def plot_fit(lam_rest, flam, flam_mod, lam_cen,
#             continuum_regions, vmax, obj_id, FWHM_v,
#             plot_fname=None,chain=None):

def plot_fit(spec,line_fitter,plot_fname=None,chain=None):

    #Set the plot x-axis limits.
    xmin = np.min(line_fitter.continuum_regions).to(u.AA).value * (1.-0.005)
    xmax = np.max(line_fitter.continuum_regions).to(u.AA).value * (1.+0.005)
    plt.xlim([xmin,xmax])

    #Plot the spectrum.
    plt.plot(spec.lam_rest,spec.flam,
             linestyle='solid',linewidth=0.5,
             color='xkcd:grey')

    #Velocity region for plotting the model.
    lam_mod = spec.lam_rest[(spec.lam_rest>xmin*u.AA) & 
                            (spec.lam_rest<xmax*u.AA)]

    #Set the model
    #flam_mod = line_fitter.flam_cont_model(lam_mod)+\
    #           line_fitter.flam_line_model(lam_mod)
    flam_mod = line_fitter.flam_model(lam_mod)

    #If a chain is provided, plot the 1-sigma regions.
    if chain is not None:
        chain_output  = np.loadtxt(chain)
        #Unfortunately we'll have to go slowly about this to not
        #trigger a memory error.
        flam_mod_low1 = np.zeros(len(lam_mod))*u.erg/u.s/u.cm**2/u.AA
        flam_mod_hig1 = np.zeros(len(lam_mod))*u.erg/u.s/u.cm**2/u.AA
        flam_mod_low2 = np.zeros(len(lam_mod))*u.erg/u.s/u.cm**2/u.AA
        flam_mod_hig2 = np.zeros(len(lam_mod))*u.erg/u.s/u.cm**2/u.AA
        for k,lam_use in enumerate(lam_mod):
            lam_mod_chain = np.tile(lam_use,chain_output.shape[0])
            flam_mod_chain = line_fitter.flam_model(lam_mod_chain,
                                                    chain_output=chain_output)
            flam_mod_low1[k], flam_mod_hig1[k] = get_error(flam_mod_chain, 
                                                           flam_mod[k],
                                                           cf=68.3)
            flam_mod_low2[k], flam_mod_hig2[k] = get_error(flam_mod_chain, 
                                                           flam_mod[k],
                                                           cf=95.4)
        plt.fill_between(lam_mod,
                         flam_mod-flam_mod_low2,
                         flam_mod+flam_mod_hig2,
                         color='xkcd:cyan',
                         alpha=1.0)
        plt.fill_between(lam_mod,
                         flam_mod-flam_mod_low1,
                         flam_mod+flam_mod_hig1,
                         color='xkcd:orange',
                         alpha=1.0)

    #Plot the model.
    plt.plot(lam_mod,flam_mod,'-b')

    #Set the y-axis limits
    flam_min = np.min(flam_mod).value
    flam_max = np.max(flam_mod).value
    dflam = flam_max-flam_min
    flam_min -= dflam*1.5
    flam_max += dflam*1.5
    plt.ylim([flam_min,flam_max])

    #Plot the continuum fit regions.
    for i in range(len(line_fitter.continuum_regions[:][0])):
        lam1 = line_fitter.continuum_regions[i][0].value
        lam2 = line_fitter.continuum_regions[i][1].value
        plt.plot([lam1,lam1],[flam_min,flam_max],'--g')
        plt.plot([lam2,lam2],[flam_min,flam_max],'--g')

    #Plot the line-fitting regions.
    v = (c*(spec.lam_rest/line_fitter.lam_cen_fit-1.)).to(u.km/u.s)
    vabs = np.abs(v)
    vmax = line_fitter.line_velocity_region
    lam_line_fit = spec.lam_rest[vabs<vmax]
    lam_line_fit_min = np.min(lam_line_fit).value
    lam_line_fit_max = np.max(lam_line_fit).value
    plt.plot([lam_line_fit_min,lam_line_fit_min],[flam_min,flam_max],'--b')
    plt.plot([lam_line_fit_max,lam_line_fit_max],[flam_min,flam_max],'--b')

    #Labels
    plt.xlabel(r'$\lambda_{Rest}\ (\AA)$')
    plt.ylabel(r'$F_{\lambda}\ (erg/cm^2/s/\AA)$')
    plt.title("{0:s}  FWHM = {1:.1f}".format(spec.name,line_fitter.FWHM_v))

    if plot_fname is None:
        plt.show()
    else:
        plt.savefig(plot_fname,dpi=300)
    plt.close()

####

def plot_full_spec(lam_rest, flam, obj_id, plot_fname=None):

    #Set the plot x-axis limits.
    xmin = np.min(lam_rest).to(u.AA).value * (1.-0.005)
    xmax = np.max(lam_rest).to(u.AA).value * (1.+0.005)
    plt.xlim([xmin,xmax])

    #Set the y-axis limits
    flam_min = np.min(flam).value
    flam_max = np.max(flam).value
    #dflam = flam_max-flam_min
    #flam_min -= dflam*1.5
    #flam_max += dflam*1.5
    plt.ylim([flam_min,flam_max])

    #Plot the emission lines
    lines_cat = open("list_of_emission_lines.txt")
    for line in lines_cat:
        if line[0]=='#':
            continue
        x = line.split()
        plot_line(float(x[0]),r'${0:s}$'.format(x[1]))

    #Plot the spectrum.
    plt.plot(lam_rest,flam,linestyle='solid',color='xkcd:grey')

    #Plot the smoothed spectrum. Use a N pixel window with a mth order
    #polynomial.
    N = 11
    m = 3
    flam_smoothed = savgol_filter(flam,N,m)
    plt.plot(lam_rest,flam_smoothed,'-k')

    #Labels
    plt.xlabel(r'$\lambda_{Rest}\ (\AA)$')
    plt.ylabel(r'$F_{\lambda}\ (erg/cm^2/s/\AA)$')
    plt.title("{0:s}".format(obj_id))

    if plot_fname is None:
        plt.show()
    else:
        plt.savefig(plot_fname)
    plt.close()

###
