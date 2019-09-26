#!/usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.constants import c
from scipy.signal import savgol_filter

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

def plot_fit(lam_rest, flam, flam_mod, lam_cen,
             continuum_regions, vmax, obj_id, FWHM_v,
             plot_fname=None):

    #Set the plot x-axis limits.
    xmin = np.min(continuum_regions).to(u.AA).value * (1.-0.005)
    xmax = np.max(continuum_regions).to(u.AA).value * (1.+0.005)
    plt.xlim([xmin,xmax])

    #Plot the spectrum.
    plt.plot(lam_rest,flam,linestyle='solid',color='xkcd:grey')
    plt.plot(lam_rest,flam_mod,'-b')

    #Set the y-axis limits
    flam_min = np.min(flam_mod).value
    flam_max = np.max(flam_mod).value
    dflam = flam_max-flam_min
    flam_min -= dflam*1.5
    flam_max += dflam*1.5
    plt.ylim([flam_min,flam_max])

    #Plot the continuum fit regions.
    for i in range(len(continuum_regions[:][0])):
        lam1 = continuum_regions[i][0].value
        lam2 = continuum_regions[i][1].value
        plt.plot([lam1,lam1],[flam_min,flam_max],'--g')
        plt.plot([lam2,lam2],[flam_min,flam_max],'--g')

    #Plot the line-fitting regions.
    v = (c*(lam_rest/lam_cen-1.)).to(u.km/u.s)
    vabs = np.abs(v)
    lam_line_fit = lam_rest[vabs<vmax]
    lam_line_fit_min = np.min(lam_line_fit).value
    lam_line_fit_max = np.max(lam_line_fit).value
    plt.plot([lam_line_fit_min,lam_line_fit_min],[flam_min,flam_max],'--b')
    plt.plot([lam_line_fit_max,lam_line_fit_max],[flam_min,flam_max],'--b')

    #Labels
    plt.xlabel(r'$\lambda_{Rest}\ (\AA)$')
    plt.ylabel(r'$F_{\lambda}\ (erg/cm^2/s/\AA)$')
    plt.title("{0:s}  FWHM = {1:.1f}".format(obj_id,FWHM_v))

    if plot_fname is None:
        plt.show()
    else:
        plt.savefig(plot_fname)
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
