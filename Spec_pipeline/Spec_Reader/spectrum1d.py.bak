import numpy as np
import astropy.units as u
from astropy.constants import c
from astropy.io import fits
import re

def read_fits_spectrum1d(file_name,
                         dispersion_unit=u.dimensionless_unscaled,
                         flux_unit = None):

    try:
        s = fits.open(file_name)
    except IOError:
        print("Cannot open file {0:s}".format(file_name))
        return

    #Get the number of axis
    NAXIS = s[0].header['NAXIS']

    #Declare the list.
    spec = []

    #If only one axis, then we have a single spectrum, all is good.
    if NAXIS==1:
        spec.append(spectrum1d(None,s,dispersion_unit,flux_unit))
        return spec

    #If not, create a dummy array over which we'll
    #iterate.
    axis_size = np.zeros(NAXIS-1,dtype=np.int32)
    for i in range(NAXIS-1):
        naxis_j = s[0].header['NAXIS{0:d}'.format(NAXIS-i)]
        axis_size[i] = naxis_j
    dummy = np.ndarray(axis_size)

    #Now, go through the axes.
    it = np.nditer(dummy,["multi_index"])
    while not it.finished:
        spec1d = spectrum1d(it.multi_index,s,dispersion_unit,flux_unit)
        spec.append(spec1d)
        it.iternext()
    s.close()
    return spec

class spectrum1d(object):

    def __init__(self,multi_index,s,_dispersion_unit,_flux_unit):
        self.data = None
        self.dispersion = None
        self.dispersion_unit = _dispersion_unit
        self.unit = _flux_unit
        self.header = None
        self.load_spec(multi_index,s)

    def load_spec(self,multi_index,s):

        #Put here the spectrum.
        if multi_index:
            self.data = s[0].data[multi_index]
        else:
            self.data = s[0].data

        #Save the header
        self.header = s[0].header

        #Now, let's figure out the dispersion.  Not sure what is the
        #correct thing here, but this works for our spectra.
        if multi_index:
            i = multi_index[0]+1
        else:
            i = 1
        #ctype = s[0].header['CTYPE{0:d}'.format(i)]
        #try:
        #    s[0].header['CRVAL{0:d}'.format(i)]
        #except:
        #    i = 1
        #crval = s[0].header['CRVAL{0:d}'.format(i)]
        #cdii  = s[0].header['CD{0:d}_{0:d}'.format(i)]
        #crpix = s[0].header['CRPIX{0:d}'.format(i)]

        #Instead of doing the hack above, this is probably more correct. Try to load all of the header values needed, and if any of them fail, then do not load anything.
        try:
            ctype = s[0].header['CTYPE{0:d}'.format(i)]
            crval = s[0].header['CRVAL{0:d}'.format(i)]
            crval = s[0].header['CRVAL{0:d}'.format(i)]
            cdii  = s[0].header['CD{0:d}_{0:d}'.format(i)]
            crpix = s[0].header['CRPIX{0:d}'.format(i)]
        except KeyError:
            return
        if ctype=="LINEAR":
            l = np.array(range(1,len(self.data)+1))
            self.dispersion = crval + cdii*(l-crpix)
            self.dispersion = self.dispersion*self.dispersion_unit
        else:
            print("Unknown CTYPE {0:s}".format(ctype))

        #Finally, setup the units if there.
        self.unit = u.Unit(s[0].header['BUNIT'])

        return
