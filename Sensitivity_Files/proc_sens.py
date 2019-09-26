#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
from specutils.io.read_fits import read_fits_spectrum1d
import astropy.units as u
from pysynphot import observation,spectrum
import matplotlib.pyplot as plt
import os

###

show_plots = True
#show_plots = False

save_plot = True
#save_plot = False.

###

##Taken from https://www.astrobetter.com/blog/2013/08/12/python-tip-re-sampling-spectra-with-pysynphot/
###

def rebin_sens(wave, specin, wavnew):
    spec = spectrum.ArraySourceSpectrum(wave=wave.to(u.AA).value, 
                                        flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave.to(u.AA).value, 
                                         f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, 
                                  binset=wavnew.to(u.AA).value, 
                                  force='taper')
    wmin = np.ceil(wave.min())
    wmax = np.floor(wave.max())
    sensx = ma.masked_where((wavnew<=wmin) | (wavnew>=wmax), obs.binflux)
    return sensx

###

#Read the master file with the Sensitivities. The first column is the
#folder name (separated by instrument and grism). The second is the
#degree of the polynomial fit to the median sensitivity if requested,
#and the final column says whether we are saving the fit to the median
#sensitivity or the median sensitivity directly.
cat = open("Sensitivities.txt")
folders = []
poly_fit_deg = []
type_save = []
for line in cat:
    if line[0]=="#":
        continue
    x = line.split()
    folders.append(x[0])
    poly_fit_deg.append(int(float(x[1])))
    type_save.append(x[2])
cat.close()

#Now, go through every folder. 
for i,folder in enumerate(folders):
    print folder
    for (dirpath,dirnames,filenames) in os.walk(folder):
        continue

    #Read each sensitivity file in the folder. 
    sens = []
    for fname in filenames:
        sens_aux = read_fits_spectrum1d(folder+"/"+fname,
                                        dispersion_unit=u.AA,
                                        flux_unit=None)
        sens.append(sens_aux)

    #Get the minimum value of all sensitivities and the maximum as
    #well.
    lmin_all = np.zeros(len(sens))
    lmax_all = np.zeros(len(sens))
    for k,s in enumerate(sens):
        lmin_all[k] = s.dispersion.min().value
        lmax_all[k] = s.dispersion.max().value
    lmin = np.ceil( lmin_all.min())*u.AA
    lmax = np.floor(lmax_all.max())*u.AA
    #Don't go too far on the UV or IR, not useful.
    if lmin<3000.*u.AA:
        lmin = 3000.*u.AA
    if lmax>11000.*u.AA:
        lmax = 11000.*u.AA
    
    #Now, we'll do rebinning of all the spectra in a single, finely
    #sampled, grid that spans the entire wavelength range of all
    #sensitivities. Mask the negative elements in the rebin function
    #we mask all of the elements that are not contained within the
    #wave grid.
    wave = np.arange(lmin.value,lmax.value,0.5)*u.AA
    sens_resampled = ma.zeros((len(sens),len(wave)))
    for k,s in enumerate(sens):
        sens_use = ma.masked_where(s.data<=0.,s.data)
        sens_resampled[k] = rebin_sens(s.dispersion,
                                       sens_use,
                                       wave)

    #Calculate the median and, if requested, the fit to it. 
    median_sens = ma.median(sens_resampled,axis=0)
    if type_save[i]=="fit":
        p = np.polyfit(wave,median_sens,poly_fit_deg[i])
        fit_sens = np.poly1d(p)

    #Plot the sensitivities.
    if show_plots is True:
        plt.title(folder)
        for k,s in enumerate(sens):
            l = s.dispersion
            f = s.data
            plt.plot(wave,sens_resampled[k],'--',linewidth=1.0)
        if type_save[i]=="fit":
            plt.plot(wave,fit_sens(wave),'-b',linewidth=3.0)
        else:
            plt.plot(wave,median_sens,'-k',linewidth=3.0)
        if save_plot is True:
            plot_name = "Sens_"+folder+".png"
            plt.savefig(plot_name)
            plt.close()
        else:
            plt.show()

    #Finally, save the sensitivities.
    fname = "Sens_"+folder+".txt"
    if type_save[i]=="fit":
        np.savetxt(fname,np.array([wave.value,fit_sens(wave)]).T)
    elif type_save[i]=="median":
        np.savetxt(fname,np.array([wave.value,median_sens]).T)
