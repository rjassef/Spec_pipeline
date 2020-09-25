#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
import astropy.units as u
import sys
import matplotlib.pyplot as plt
import re

import os
os.environ['PYSYN_CDBS'] = "."
import warnings
warnings.simplefilter("ignore")

from pysynphot import observation,spectrum

from Spec_pipeline.Spec_Reader.iraf_spectrum1d import read_fits_spectrum1d

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

##

if len(sys.argv)<=4:
    print("Correct use: python",sys.argv[0],"outputname.txt polyfit_degree lam_norm sens1 [sens2 [...sensN]]")
    sys.exit()

fname = sys.argv[1]
if fname[-4:] == 'fits':
    print("Output file must end in .txt, not in .fits.")
    sys.exit()

try:
    polyfit_degree = int(float(sys.argv[2]))
except ValueError:
    print("polyfit_degree must be a number.")
    sys.exit()

try:
    lam_norm = float(sys.argv[3])*u.AA
except ValueError:
    print("lam_norm must be a number in angstroms.")
    sys.exit()

#Read the sensitivity curves.
sens = list()
lmin = None
lmax = None
for k in range(4,len(sys.argv)):
    sens_aux = read_fits_spectrum1d(sys.argv[k], dispersion_unit=u.AA, flux_unit=None)
    sens.append(sens_aux[0])
    # plt.plot(sens_aux[0].dispersion,sens_aux[0].data,'-')
    # plt.title(re.sub("_"," ",sys.argv[k]))
    # plt.show()
    if lmin is None or np.min(sens_aux[0].dispersion)<lmin:
        lmin = np.min(sens_aux[0].dispersion)
    if lmax is None or np.max(sens_aux[0].dispersion)>lmax:
        lmax = np.max(sens_aux[0].dispersion)

#If only one curve, then just print it.
if len(sens)==1:
    plt.plot(sens[0].dispersion,sens[0].data,'b-')
    #plt.show()
    plt.savefig(fname+".png")

    np.savetxt(fname+".txt",np.array([sens[0].dispersion.value, sens[0].data*0.01]).T)
    sys.exit()

#Limit the range to 3000 - 11000AA, as we do not really cover more than that with spec data.
lmin = np.max([lmin.to(u.AA).value,3000.])*u.AA
lmax = np.min([lmax.to(u.AA).value,11000.])*u.AA

#Put them in the same wavelength grid.
wave = np.arange(np.min(lmin).value,np.max(lmax).value,0.5)*u.AA
sens_resampled = ma.zeros((len(sens),len(wave)))
for i,s in enumerate(sens):
    sens_use = ma.masked_where(s.data<=0.,s.data)
    sens_resampled[i] = rebin_sens(s.dispersion, sens_use, wave)

#Normalize all the curves at 4750A or 8000A, depending on which is observed.
# if lmin<4750.*u.AA:
#     lam_norm = 4750.*u.AA
# else:
#     lam_norm = 8000.*u.AA
ktarg = np.argmin(np.abs(wave.to(u.AA).value-lam_norm.to(u.AA).value))
for i in range(1,len(sens)):
    norm = sens_resampled[0][ktarg]/sens_resampled[i][ktarg]
    sens_resampled[i] *= norm

#Calculate the median.
median_sens = ma.median(sens_resampled,axis=0)

#Fit the median to obtain an smooth curve.
p = np.polyfit(wave,median_sens,polyfit_degree)
fit_sens = np.poly1d(p)

#Name of the combined sensitivity curve.
#fname = "Sens_LRIS_e2v_600-4000_560_b"

#Plot all the curves an the median.
for s in sens_resampled:
    plt.plot(wave,s,'--')
plt.plot(wave,fit_sens(wave.value),'b-')
#plt.show()
plt.savefig(fname+".png")

#Save the final curve.
np.savetxt(fname+".txt",np.array([wave.value, fit_sens(wave.value)*0.01]).T)
