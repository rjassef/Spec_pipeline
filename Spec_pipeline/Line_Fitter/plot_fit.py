#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text',usetex=True)

import astropy.units as u
from astropy.constants import c
from scipy.signal import savgol_filter
import re

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

def plot_fit(spec,line_fitter,plot_fname=None,chain_file=None,chain=None):

    #Open the Figure. To not clash with figures already open, find the latest figure open, and open a new one.
    if len(plt.get_fignums())==0:
        nfig = 1
    else :
        nfig = plt.get_fignums()[-1]+1
    fig = plt.figure(nfig)
    ax = fig.add_subplot(111)
    #fig, ax = plt.subplots()

    #Set the plot x-axis limits.
    i_line = line_fitter.get_i_line(spec)
    i_cont = line_fitter.get_i_cont(spec)
    xminl = np.min(spec.lam_rest[i_line])
    xmaxl = np.max(spec.lam_rest[i_line])
    xminc = np.min(spec.lam_rest[i_cont])
    xmaxc = np.max(spec.lam_rest[i_cont])
    xmincr = np.min(line_fitter.continuum_regions)
    xmaxcr = np.max(line_fitter.continuum_regions)
    xmin = np.min([xminl.value,xminc.value,xmincr.value])
    xmax = np.max([xmaxl.value,xmaxc.value,xmaxcr.value])
    plt.xlim([xmin,xmax])

    #Plot the spectrum.
    plt.plot(spec.lam_rest,spec.flam,
             linestyle='solid',linewidth=0.5,
             color='xkcd:grey')
    plt.plot(spec.lam_rest,spec.flam_err,linestyle='dotted',linewidth=0.5,color='xkcd:grey')

    #Velocity region for plotting the model.
    lam_mod = spec.lam_rest[(spec.lam_rest>xmin*u.AA) &
                            (spec.lam_rest<xmax*u.AA)]

    #Set the model
    flam_mod = line_fitter.flam_model(lam_mod)
    flam_cont_mod = line_fitter.flam_cont_model(lam_mod)

    #If a chain is provided, plot the 1-sigma regions.
    if chain is not None or chain_file is not None:
        if chain_file is not None:
            chain_output  = np.loadtxt(chain_file)
        else:
            chain_output = chain
        #Unfortunately we'll have to go slowly about this to not
        #trigger a memory error.
        flam_mod_low1 = np.zeros(len(lam_mod))*line_fitter.flamunit
        flam_mod_hig1 = np.zeros(len(lam_mod))*line_fitter.flamunit
        flam_mod_low2 = np.zeros(len(lam_mod))*line_fitter.flamunit
        flam_mod_hig2 = np.zeros(len(lam_mod))*line_fitter.flamunit
        for k,lam_use in enumerate(lam_mod):
            lam_mod_chain = np.tile(lam_use,chain_output.shape[0])
            flam_mod_chain = line_fitter.flam_model(lam_mod_chain,
                                                    chain_output=chain_output)
            flam_mod_low1[k], flam_mod_hig1[k] = get_error(flam_mod_chain,
                                                           flam_mod[k],
                                                           cf=68.3)
            flam_mod_low2[k], flam_mod_hig2[k] = get_error(flam_mod_chain,
                                                           flam_mod[k],
                                                           cf=99.1)#cf=95.4)
        plt.fill_between(lam_mod,
                         flam_mod-flam_mod_low2,
                         flam_mod+flam_mod_hig2,
                         color='xkcd:cyan',
                         alpha=1.0, label="3$\sigma$")
        plt.fill_between(lam_mod,
                         flam_mod-flam_mod_low1,
                         flam_mod+flam_mod_hig1,
                         color='xkcd:orange',
                         alpha=1.0, label="1$\sigma$")

    #Plot the model.
    plt.plot(lam_mod,flam_cont_mod,'--r',label="Continuum")
    plt.plot(lam_mod,flam_mod,'-b',label="Best-fit")

    #Set the y-axis limits
    flam_min = np.min(flam_mod).value
    flam_max = np.max(flam_mod).value
    dflam = flam_max-flam_min
    flam_min -= dflam*1.5
    flam_max += dflam*1.5
    plt.ylim([flam_min,flam_max])

    #Plot the continuum fit regions.
    #for i in range(len(line_fitter.continuum_regions[:][0])):
    for i in range(line_fitter.ncont_reg):
        lam1 = line_fitter.continuum_regions[i][0].value
        lam2 = line_fitter.continuum_regions[i][1].value
        plt.plot([lam1,lam1],[flam_min,flam_max],'--g')
        plt.plot([lam2,lam2],[flam_min,flam_max],'--g')

    #Plot the line-fitting regions.
    lam_line_fit = spec.lam_rest[i_line]
    lam_line_fit_min = np.min(lam_line_fit).value
    lam_line_fit_max = np.max(lam_line_fit).value
    plt.plot([lam_line_fit_min,lam_line_fit_min],[flam_min,flam_max],'--b')
    plt.plot([lam_line_fit_max,lam_line_fit_max],[flam_min,flam_max],'--b')

    #Labels
    plt.xlabel(r'$\lambda_{Rest}\ (\AA)$')
    plt.ylabel(r'$F_{\lambda}\ (erg/cm^2/s/\AA)$')
    plt.title("{0:s} z={1:.3f} {2:s}".format(spec.name, spec.zspec, re.sub("_","\_",line_fitter.line_name)))
    plt.legend(loc='upper right')

    try:
        line_fitter.line_legends(ax)
    except AttributeError:
        pass

    # #Textbox with fit results. We will use one box per emission line, below the spectrum. We have up to three emission lines.
    # line_name = re.sub("red","",line_fitter.line_name)
    # lname = line_name.split("_")
    # if len(lname)<line_fitter.nlines:
    #     lname = [lname[0]]*line_fitter.nlines
    # props = dict(boxstyle='round', facecolor='white', alpha=0.75)
    # for i in range(line_fitter.nlines):
    #     textbox = "{0:s}".format(lname[i])
    #     textbox += "\nFWHM = {0:.0f}".format(line_fitter.FWHM_v[i])
    #     textbox += "\n$\Delta v$ = {0:.1f}".format(line_fitter.dv_fit[i])
    #     if line_fitter.MC_chain is not None:
    #         textbox += "\nSNR = {0:.1f}".format(line_fitter.line_SNR[i])
    #     if line_fitter.p is not None:
    #         textbox += "\np   = {0:.3f}".format(line_fitter.p[i])
    #     plt.text(0.025+i/3., 0.025, textbox, transform=ax.transAxes, fontsize=10, verticalalignment='bottom', horizontalalignment='left', bbox=props)


    if plot_fname is None:
        plt.show()
    else:
        plt.savefig(plot_fname,dpi=300)
        plt.close(nfig)

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
    lines_cat.close()

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
