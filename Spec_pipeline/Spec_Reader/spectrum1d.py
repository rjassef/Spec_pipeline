import numpy as np
import astropy.units as u
from astropy.constants import c
from astropy.io import fits
import re
import sys
from warnings import warn

def read_fits_spectrum1d(file_name,
                         dispersion_unit=None, #u.dimensionless_unscaled,
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

    def __init__(self,multi_index,s,dispersion_unit,flux_unit=None):
        self.data = None
        self.dispersion = None
        self.dispersion_unit = dispersion_unit
        self.native_dispersion_unit = None
        self.unit = flux_unit
        self.native_unit = None
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

        #Try reading the wavelength units from the header. If no units are found on the headers, assume that they native units are the units in which the output is requested. If no output units are requested, assume that the native units are requested.
        try:
            wav_unit = re.search("units=([^\s]*)",s[0].header['WAT1_001'])[1]
            if wav_unit=="angstroms":
                wav_unit = "Angstrom"
            self.native_dispersion_unit = u.Unit(wav_unit)
            if self.dispersion_unit is None:
                self.dispersion_unit = self.native_dispersion_unit
        except (KeyError, ValueError):
            if self.dispersion_unit is not None:
                self.native_dispersion_unit = self.dispersion_unit
            else:
                print("No WAT1_001 in headers and no output dispersion units provided.")
                print("Using dimensionless_unscaled")
                self.dispersion_unit = u.dimensionless_unscaled
                self.native_dispersion_unit = u.dimensionless_unscaled

        #Instead of doing the hack above, this is probably more correct. Try to load all of the header values needed, and if any of them fail, then do not load anything.
        try:
            ctype = s[0].header['CTYPE{0:d}'.format(i)]
            crval = s[0].header['CRVAL{0:d}'.format(i)]
            cdii  = s[0].header['CD{0:d}_{0:d}'.format(i)]
            crpix = s[0].header['CRPIX{0:d}'.format(i)]
        except KeyError:
            warn("Could not read WCS headers for axis {0:d}".format(i))
            return
        if ctype=="LINEAR":
            l = np.array(range(1,len(self.data)+1))
            self.dispersion = crval + cdii*(l-crpix)
            self.dispersion = self.dispersion*self.native_dispersion_unit
            self.dispersion = self.dispersion.to(self.dispersion_unit)
        else:
            print("Unknown CTYPE {0:s}".format(ctype))

        #Finally, setup the flux units. If there are any in headers, do not assign any native units. If no requested units, then set the requested output units to the native flux units.
        try:
            self.native_unit = u.Unit(s[0].header['BUNIT'])
        except KeyError:
            print("No BUNIT in headers and no flux units provided.")
            print("Using dimensionless_unscaled")
            self.native_unit = u.dimensionless_unscaled

        if self.unit is None:
            self.unit = self.native_unit

        #Check if units are equivalent. If they are, apply conversion factor if any. If not, crash.
        try:
            self.native_unit.to(self.unit)
        except u.UnitConversionError:
            #Try to see if it is in flam and we are requesting fnu.
            try:
                self.native_unit.to(self.unit*u.Hz/u.AA)
                self.data *= self.dispersion**2/c
            except u.UnitConversionError:
                #Try the opposite. Maybe we are asking for flam and it is in fnu.
                try:
                    self.native_unit.to(self.unit*u.AA/u.Hz)
                    self.data *= c/self.dispersion**2
                except u.UnitConversionError:
                    print("Error: Flux units requested not equivalent to flux units in header.")
                    sys.exit()
        self.data *= self.native_unit.to(self.unit)

        return
