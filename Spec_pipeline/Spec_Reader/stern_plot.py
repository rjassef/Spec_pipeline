import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import astropy.units as u
import re
import os

from ..Line_Fitter.multi_line import Multi_Line_fit

#It would be great to have a more elegant way to convert from the code names of the lines to the latex pretty version, but we'll just do this for now.
linename_latex = {
    "Lyb"     : r"Ly$\beta$",
    "LyA"     : r"Ly$\alpha$",
    "Hbeta"   : r"H$\beta$",
    "Ha"      : r"H$\alpha$",
    "Hgamma"  : r"H$\gamma$",
    "Hdelta"  : r"H$\delta$",
    "Hepsilon": r"H$\epsilon$"
}

def stern_plot(specs, date, em_lines_list=None, sv_wl=21, sv_polyorder=5, hardcopy=None, flam_units=(u.erg/u.s/u.cm**2/u.AA), lam_units=u.AA, legend_inside=False):
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
    """

    #If no list of emission lines has been provided, load the full list. This is a bad idea though, much better to provide them.
    if em_lines_list is None:
        em_lines = np.genfromtxt("{}/Spec_pipeline/Line_Fitter/multi_lines.txt".format(os.environ.get('SPEC_PIPE_LOC')), usecols=[0], dtype='U')

    #Find the minimum line center and create all the line objects.
    em_lines = list()
    for k, em_line_name in enumerate(em_lines_list):
        #Start by loading the emission line in question.
        em_lines.append(Multi_Line_fit(em_line_name))
        #Check if we have a minimum.
        if k==0 or line_center_min > np.min(em_lines[-1].line_center):
            line_center_min = np.min(em_lines[-1].line_center)

    #Create the figure.
    fig, ax = plt.subplots()

    #Smooth the spectra.
    for k, spec in enumerate(specs):
        spec.flam_smooth = savgol_filter(spec.flam.to(flam_units).value, sv_wl, sv_polyorder)

    #Find the minimum and maximum value of F_lam and lambda for setting the plot x and y range. Only consider wavelengths longer than about the shortest wavelength in the list of line centers.
    for k, spec in enumerate(specs):
        cond = spec.lam_rest>0.999*line_center_min
        if k==0 or ymax<np.max(spec.flam_smooth[cond]):
            if cond.sum()>0:
                kw_max = np.argmax(spec.flam_smooth[cond])
                ymax = spec.flam_smooth[cond][kw_max]
                x_fmax = spec.lam_obs[cond][kw_max].to(lam_units).value
            else:
                ymax = 0
        if k==0 or xmin>np.min(spec.lam_obs):
            xmin = np.min(spec.lam_obs)
        if k==0 or xmax<np.min(spec.lam_obs):
            xmax = np.max(spec.lam_obs)
    xmin   = xmin.to(lam_units).value
    xmax   = xmax.to(lam_units).value

    #For the y range, set it to be 20% higher than the highest point in the plot. The minimum should be -5% of the maximum.
    ymax *= 1.2
    ymin  = -0.05*ymax
    ax.set_ylim([ymin, ymax])

    #Set the xlimits to match the spectrum range.
    ax.set_xlim([0.95*xmin, 1.02*xmax])

    #Mark the emission/absorption lines.
    for k, em_line in enumerate(em_lines):

        #Set the observed-frame line center.
        line_center = em_line.line_center.value*(1+specs[0].zspec)

        #The label with the line name should be either half way through the plot or at the top, depending on the flux values of the emission lines.
        yline = 0.5*(ymax+ymin)
        flam_max = None
        for spec in specs:
            box_width = 0.005
            cond = (spec.lam_obs.value>=np.min(line_center)*(1-box_width)) & (spec.lam_obs.value<=np.max(line_center)*(1+box_width))
            if len(spec.flam_smooth[cond])>0:
                flam_max = np.max(spec.flam_smooth[cond])
        if flam_max is not None and flam_max >= yline:
            yline = ymax

        #Go through the emission lines and separate the names of the multiple fit emission lines.
        if re.search("_", em_line.line_name):
            lnames = re.sub("red","",em_line.line_name)
            line_names = lnames.split("_")
        else:
            #The [NeV] lines are separate enough that we want to label both for clarity. Should do something more general based on the wavelength separation of the emission lines.
            if em_line.line_name == '[NeV]':
                line_names = [em_line.line_name]*em_line.nlines
            else:
                line_names = [em_line.line_name]#*em_line.nlines

        #Draw a vertical dashed gray line with a label at the location of every emission line within the x-axis range. For multiple line fits. draw the label of the first one on the left, and the rest on the right.
        #for k, lam in enumerate(line_center):
        for i,j in enumerate(np.argsort(line_center)):
            lam = line_center[j]
            if lam<xmin or lam>xmax:
                continue

            #If it is Hepsilon, push ymax down a bit to avoid clashing with other lines.
            yline_use = yline
            if len(line_names)>j and line_names[j]=='Hepsilon':
                yline_use = 0.6*yline

            ax.plot([lam, lam], [ymin, yline_use], linestyle='dashed', color='xkcd:grey', linewidth=0.5)

            #If the line is a complex with only one designation, like [OII] or [OIII], only plot the line name once.
            if len(line_names)==1 and j>0:
                continue

            #Otherwise, draw the label
            dlam_label = -0.015*(xmax-xmin) #120
            if i>=1 and line_names[j]!='[NeV]':
                dlam_label = 0.003*(xmax+xmin) #50
            line_name_use = line_names[j]
            if line_names[j] in linename_latex:
                line_name_use = linename_latex[line_names[j]]
            ax.text(lam+dlam_label, 0.9*yline_use, line_name_use, rotation=90, fontsize=7, clip_on=True)

    #Draw a horizontal dashed line at the zero flux level.
    ax.plot([xmin, xmax], [0,0], linestyle='dashed', color='xkcd:grey', linewidth=0.5)

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
    xlabel = re.sub("Angstrom","\AA", xlabel)
    ax.set_xlabel(xlabel, fontsize=12)
    ylabel = r"$F_{{\lambda}}$ [{}]".format(flam_units)
    ylabel = re.sub("Angstrom","\AA~", ylabel)
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
