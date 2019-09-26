import numpy as np
from pysynphot import observation
from pysynphot import spectrum
import astropy.units as u

####
#Taken from https://www.astrobetter.com/blog/2013/08/12/python-tip-re-sampling-spectra-with-pysynphot/
####

def rebin_spec(wave, specin, wavnew):
    spec = spectrum.ArraySourceSpectrum(wave=wave.to(u.AA).value, 
                                        flux=specin.value)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave.to(u.AA).value, 
                                         f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, 
                                  binset=wavnew.to(u.AA).value, 
                                  force='taper')

    return obs.binflux*specin.unit
