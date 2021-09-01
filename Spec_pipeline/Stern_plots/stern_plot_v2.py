import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import astropy.units as u
import re
import os
import json

from ..Line_Fitter.multi_line import Multi_Line_fit

def stern_plot(specs, date, sp_lines_list=None, sv_wl=21, sv_polyorder=5, hardcopy=None, flam_units=(u.erg/u.s/u.cm**2/u.AA), lam_units=u.AA, legend_inside=False, xrange=None, yrange=None, sp_lines_conf=None):
    """
    This function receives a list of spec objects and makes a Stern-style spectrum plot with the possible emission/absorption lines marked. Note that all lines in em_lines are marked, regardless of whether they where found.

    Parameters
    ----------
    specs: list
        List of Spec objects from the Spec_pipeline.Spec_Reader.Spec class.

    date: str
        UT Date on which the observations were carried out.

    em_lines_list: list, optional
        List of emission line names to plot. Default is to load all the lines in Line_Fitter/multi_lines.txt.

    sv_wl: int, optional
        Window length for scipy.signal.savgol_filter smoothing. Default is 21.

    sv_polyorder: int, optional
        Polynomial order for scipy.signal.savgol_filter smoothing. Default is 5.

    hardcopy: str, optional
        File name of the hard copy of the plot. Can be any format supported by matplotlib.

    flam_units: astropy.units, optional
        Flux density per unit wavelength units on which to plot the spectrum. Default is (u.erg/u.s/u.cm**2/u.AA).

    lam_units: astropy.units, optional
        Units for plotting the observed wavelength. Default is u.AA .

    legend_inside: boolean, optional
        If True, puts the legend inside the plot. Default is False.

    xrange: numpy array with astropy.units, optional
        Observed wavelength range to plot. If None, full range covered by the spectrum is shown. Default is None.

    yrange: numpy array with astropy.units, optional
        Flam range range to plot. If None, range is autoscaled to the spectrum. Default is None.

    sp_lines_conf: dictionary, optional
        Dictionary with modifiers from default behavior. Default is None.

    """

    #If no wavelength range has been given, find the wavelength range of the spectra.
    if xrange is None:
        lam_rest = specs[0].lam_rest
        for spec in specs[1:]:
            lam_rest = np.concatenate([lam_rest, spec.lam_rest])
        xmin = np.min(lam_rest)
        xmax = np.max(lam_rest)
    else:
        xmin = xrange[0]
        xmax = xrange[1]

    #Smooth the spectra.
    for spec in specs:
        spec.flam_smooth = savgol_filter(spec.flam.to(flam_units).value, sv_wl, sv_polyorder)

    #Load the default list of spectral lines. 
    f = open("{}/Spec_pipeline/Stern_plots/default_lines.json".format(os.environ.get('SPEC_PIPE_LOC')))
    sp_lines = json.load(f)

    #Now, process the spectral lines list. For each spectral line we will need to add the new keywords, if provided, and determine whether they are to be skipped or not, and where the labels should be drawn. Keep track of the maximum line height.
    ymax = None
    for sp_line in sp_lines:
        #Load the user provided attributes.
        lid = sp_line['id']
        if sp_lines_conf is not None and lid in sp_lines_conf:
            for key in sp_lines_conf[lid].keys():
                sp_line[key] = sp_lines_conf[lid][key]

        #Add the units to the line central wavelength.
        if 'lam_unit' not in sp_line:
            sp_line['lam_unit'] = u.AA
        else:
            sp_line['lam_unit'] = u.Unit(sp_line['lam_unit'])
        sp_line['lam_rest'] = sp_line['lam_rest']*sp_line['lam_unit']

        #Add skip false is not loaded already.
        if 'skip' not in sp_line:
            sp_line['skip'] = False

        #Determine if the line is outside the range, and hence should be skipped. 
        if sp_line['lam_rest']<xmin or sp_line['lam_rest']>xmax:
            sp_line['skip'] = True

        #Do not process height for lines that should be skipped. 
        if sp_line['skip']:
            continue

        #Calculate the peak of the emission line.
        #lw = np.array([0.999, 1.001])*sp_line['lam_rest']
        lw = np.array([-10., 10.])*u.AA+sp_line['lam_rest']
        for spec in specs:
            cond = (spec.lam_rest>lw[0]) & (spec.lam_rest<lw[1])  
            if cond.sum()==0:
                continue
            sp_line['peak'] = np.max(spec.flam_smooth[cond])
            if ymax is None or sp_line['peak']>ymax:
                ymax = sp_line['peak']
                x_fmax = sp_line['lam_rest']
        if 'peak' not in sp_line:
            sp_line['skip'] = True

    if yrange is None:
        ymax *= 1.2
        ymin = -0.05*ymax
    else:
        ymax = yrange[1]
        ymin = yrange[0]

    #Create the figure.
    fig, ax = plt.subplots()

    #Set the axis ranges.
    xmin *= (1+specs[0].zspec)
    xmax *= (1+specs[0].zspec)
    ax.set_xlim([0.95*xmin.to(lam_units).value, 1.02*xmax.to(lam_units).value])
    ax.set_ylim([ymin, ymax])

    #Draw the spectral lines.
    for sp_line in sp_lines:

        if sp_line['skip']:
            continue

        if 'peak' not in sp_line:
            print(sp_line['id'])
            input()

        #Set the absolute y-position.
        if 'ypos' in sp_line:
            yline = sp_line['ypos']
        else:
            yline = 0.5*(ymax+ymin)
            if sp_line['peak']>yline:
                yline = ymax

        #Draw a line at the central wavelength.
        lam = sp_line['lam_rest']*(1+specs[0].zspec)
        linestyle='dashed'
        ax.plot(np.ones(2)*lam, [ymin, yline], linestyle=linestyle, color='xkcd:grey', linewidth=0.5)

        #Now, draw the label.
        if 'name' not in sp_line or sp_line['name'] is None:
            continue

        #Now, if a relative offset is set in y, then apply it.
        if 'ypos_rel' in sp_line:
            yline *= sp_line['ypos_rel']

        #Set the label left or right of the line.
        if 'xpos_rel' not in sp_line:
            if 'xpos' not in sp_line or sp_line['xpos']=='left':
                sp_line['xpos_rel'] = -0.015
            else:
                sp_line['xpos_rel'] = 0.005
        dlam_label = sp_line['xpos_rel']*(xmax-xmin)

        #Draw the label
        lam_label = (lam+dlam_label).to(lam_units).value
        ax.text(lam_label, 0.99*yline, sp_line['name'], rotation=90, fontsize=7, clip_on=True, va='top')

    #Draw a horizontal dashed line at the zero flux level.
    ax.plot(ax.get_xlim(), [0,0], linestyle='dashed', color='xkcd:grey', linewidth=0.5)

    #Plot the smoothed version of the spectrum.
    for k, spec in enumerate(specs):
        ax.plot(spec.lam_obs, spec.flam_smooth, 'k-', linewidth=0.5)

    #Make the ticks point inwards.
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)

    #Make the legend box.
    date_legend = "{} - UT {}".format(specs[0].instrument, date)
    if hasattr(specs[0],'dzspec'):
        z_legend = "$z = {0:.4f}\pm{1:.4f}$".format(specs[0].zspec, specs[0].dzspec)
    else:
        z_legend = "$z = {0:.3f}$".format(specs[0].zspec)
    if legend_inside:
        yloc   = 0.95
        dyloc  = 0.06
        va_leg = 'top'
        if x_fmax > 0.5*(xmax+xmin):
            xloc   = 0.05
            ha_leg = 'left'
        else:
            xloc   = 0.95
            ha_leg = 'right'
        xdate = xloc
        ydate = yloc-dyloc
        xz    = xloc
        yz    = yloc-1.8*dyloc
        ha_date = ha_leg
        va_date = va_leg
        ha_z    = ha_leg
        va_z    = va_leg
    else:
        xloc = 0.5
        yloc = 1.14
        xdate = 0.50 #0.30
        ydate = 1.085
        xz    = 0.50 #0.70
        yz    = 1.05
        ha_leg  = 'center'
        va_leg  = 'top'
        ha_date = 'center'
        va_date = 'top'
        ha_z    = 'center'
        va_z    = 'top'
    ax.annotate(specs[0].name, xy=(xloc , yloc ), xycoords='axes fraction', ha=ha_leg , va=va_leg , size=14)
    ax.annotate(date_legend  , xy=(xdate, ydate), xycoords='axes fraction', ha=ha_date, va=va_date, size=10)
    ax.annotate(z_legend     , xy=(xz   , yz   ), xycoords='axes fraction', ha=ha_z   , va=va_z   , size=11)


    #Draw the axis labels
    xlabel = r"Observed Wavelength ({})".format(lam_units)
    xlabel = re.sub("Angstrom",r"\\AA", xlabel)
    ax.set_xlabel(xlabel, fontsize=12)
    ylabel = r"$F_{{\lambda}}$ [{}]".format(flam_units)
    ylabel = re.sub("Angstrom",r"\\AA~", ylabel)
    ylabel = re.sub("cm2", "cm$^2$", ylabel)
    ax.set_ylabel(ylabel, fontsize=12)

    #Show or save the plot.
    if hardcopy is None:
        plt.show()
    else:
        fig.savefig(hardcopy, dpi=200)

    #Close the figure.
    plt.close(fig)

    return


