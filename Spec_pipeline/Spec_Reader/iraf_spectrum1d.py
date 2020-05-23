import numpy as np
import astropy.units as u
from astropy.constants import c
from astropy.io import fits
import re
from warnings import warn
import shlex
import sys

#This function is written following the wcs guidelines from the iraf specwcs routine: https://iraf.net/irafhelp.php?val=specwcs&help=Help+Page.


#This function returns a list, where each elemenf of the list is an spectrum1d object.
def read_fits_spectrum1d(file_name,
                         dispersion_unit=None, #u.dimensionless_unscaled,
                         flux_unit = None):

    #Try to open the file. Return none if the file is not found.
    try:
        s = fits.open(file_name)
    except IOError:
        print("Cannot open file {0:s}".format(file_name))
        return None


    #Get the number of axes
    NAXIS = s[0].header['NAXIS']

    #Declare the list.
    spec = []

    #Create an array over which we'll be iterating to load the spectra from all the axes. Note that, according to the IRAF specifications, NAXIS1 is always the dispersion axis. Specifically it says that in 3D images the aperture index wraps around the lowest non-dispersion axis.
    if NAXIS==1:
        axis_size = np.ones(1,dtype=np.int32)
    else:
        axis_size = np.zeros(NAXIS-1,dtype=np.int32)
        for i in range(NAXIS-1):
            axis_size[i] = s[0].header['NAXIS{0:d}'.format(NAXIS-i)]
    dummy = np.ndarray(axis_size)

    #Now, go through the axes and load the spectra.
    it = np.nditer(dummy,["multi_index"])
    while not it.finished:
        spec.append(spectrum1d(it.multi_index,s,dispersion_unit,flux_unit))
        it.iternext()
    s.close()

    return spec


#####

def parse_wat(i,header):

    #Put together the whole WAT header.
    aux = ""
    k = 1
    while k>0:
        try:
            aux += header['WAT{0:d}_{1:03d}'.format(i,k)]
            #I don't like this, but there seems to be a bug in astropy.io.fits that removes a trailing space in the first line of the WAT2_001 header that is needed to reconstruct the proper values.
            if len(header['WAT{0:d}_{1:03d}'.format(i,k)])<68:
                aux += " "
            k+=1
        except KeyError:
            break
    #If the WAT0 header is not there, we'll end up here.
    if k==1:
        return None

    #For the WAT2 or higher, we need to remove some certain unneeded spaces around equal signs.
    aux = re.sub(" = ","=",aux)

    #Assign all of the values found to a dictionary.
    wat = dict()
    for x in shlex.split(aux):
        m = re.match("(.*)=(.*)",x)
        wat[m.group(1)] = m.group(2)

    return wat


#####

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
        if s[0].header['NAXIS']==1:
            self.data = s[0].data
        else:
            self.data = s[0].data[multi_index]

        #Save the header
        self.header = s[0].header

        #Now, let's figure out the dispersion. Define i as the spectrum number.
        i = multi_index[0]+1

        #Load the wavelength dispersion
        lam = load_lam(i,s)

        #Get the units from the WAT1 header. If not there, then we'll assume the requested units, although this should never be the case.
        wat1 = parse_wat(1,s[0].header)
        try:
            wav_unit = wat1['units']
            if wav_unit=='angstroms':
                wav_unit = "Angstrom"
            self.native_dispersion_unit = u.Unit(wav_unit)
            if self.dispersion_unit is None:
                self.dispersion_unit = self.native_dispersion_unit
        except (KeyError, ValueError, TypeError):
            if self.dispersion_unit is not None:
                self.native_dispersion_unit = self.dispersion_unit
            else:
                print("No WAT1 headers and no output dispersion units provided.")
                print("Using dimensionless_unscaled")
                self.dispersion_unit = u.dimensionless_unscaled
                self.native_dispersion_unit = u.dimensionless_unscaled
        self.dispersion = lam*self.native_dispersion_unit
        self.dispersion = self.dispersion.to(self.dispersion_unit)

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
            conv = 1.0*u.dimensionless_unscaled
        except u.UnitConversionError:
            #Try to see if it is in flam and we are requesting fnu.
            try:
                conv = self.dispersion**2/c
                (self.native_unit*conv.unit).to(self.unit)
                self.data *= conv.value
            except u.UnitConversionError:
                #Try the opposite. Maybe we are asking for flam and it is in fnu.
                try:
                    conv = c/self.dispersion.to(u.AA)**2
                    (self.native_unit*conv.unit).to(self.unit)
                    self.data *= conv.value
                except u.UnitConversionError:
                    print("Error: Flux units requested not equivalent to flux units in header.")
                    sys.exit()
        self.data *= (self.native_unit*conv.unit).to(self.unit)

        return

####
def load_lam(i,s):

    #There should be a WAT0 header. First reconstruct the whole WAT0 header line, and then find whether we are dealing with an equispec or a multispec case.
    wat0 = parse_wat(0,s[0].header)

    #If the spectrum system is equispec, then all the spectra have the same dispersion solution. If not wat0 header, then we'll assume equispec.
    if wat0 is None or wat0['system']=='equispec':
        iuse = 1

    #If instead it is a multispec system, then will need the spectrum number.
    elif wat0['system']=='multispec':
        iuse = i

    #Now, we need to reconstruct the wavelength axis. We first try to find the entry for speci in WAT2. If WAT2 is not there, we'll check CTYPEi.
    wat2 = parse_wat(2,s[0].header)
    spec_i = 'spec{0:d}'.format(iuse)
    if wat2 is not None and spec_i in wat2:

        #Figure out the physical coordinates. We'll need them for all cases.
        LTM1_1 = 1.0
        if 'LTM1_1' in s[0].header:
            LTM1_1 = s[0].header['LTM1_1']
        LTV1 = 0.0
        if 'LTV1' in s[0].header:
            LTV1 = s[0].header['LTV1']
        l = np.arange(1,s[0].header['NAXIS1']+1)
        p = (l - LTV1) / LTM1_1

        #Get the parameters for the fitting functions.
        spN = parse_specN(wat2[spec_i])

        #Based on the params, we have a number of possible cases, all described in the specwcs help file. Let's implement them.

        #Linear
        if spN['dtype']==0:
            lam = (spN['w1'] + spN['dw'] * (p - 1)) / (1 + spN['z'])

        #Log-Linear
        elif spN['dtype']==1:
            lam = 10 ** ((spN['w1'] + spN['dw'] * (p - 1)) / (1 + spN['z']))

        #Non-LINEAR
        #w = sum from i=1 to nfunc {wt_i * (w0_i + W_i(p)) / (1 + z)}
        elif spN['dtype']==2:
            lam = np.zeros(s[0].header['NAXIS1'])
            for nfunc in range(1,spN['nfuncs']+1):

                #Chebyshev
                if spN['ftype_{0:d}'.format(nfunc)]==1:
                    pmin = spN['pmin_{0:d}'.format(nfunc)]
                    pmax = spN['pmax_{0:d}'.format(nfunc)]
                    n = (p - (pmax+pmin)/2.) / ((pmax-pmin) / 2.)
                    W = 0.
                    for i,coeff in enumerate(spN['coeffs_{0:d}'.format(nfunc)]):
                        W += coeff * cheby(i+1,n)
                lam += spN['wt_{0:d}'.format(nfunc)] * (spN['w0_{0:d}'.format(nfunc)] + W) / (1.0+spN['z'])

    #If not WAT2, then we'll need to get the dispersion values from the standard WCS headers. Only implemented for LINEAR right now.
    else:
        wat1 = parse_wat(1,s[0].header)
        if (wat1 is not None and wat1['wtype']=='linear') or (wat2 is not None and wat2['wtype']=='linear') or (s[0].header['CTYPE{0:d}'.format(iuse)=='LINEAR']):

            try:
                crval = s[0].header['CRVAL{0:d}'.format(iuse)]
                cdii  = s[0].header['CD{0:d}_{0:d}'.format(iuse)]
                crpix = s[0].header['CRPIX{0:d}'.format(iuse)]
            except KeyError:
                warn("Could not read WCS headers for axis {0:d}".format(iuse))
                return None
            l = np.array(range(1,s[0].header['NAXIS1']+1))
            lam = crval + cdii*(l-crpix)

        else:
            return None

    return lam

#specN = ap beam dtype w1 dw nw z aplow aphigh [functions_i]
#function_i =  wt_i w0_i ftype_i [parameters] [coefficients]
def parse_specN(text):
    x = [float(ix) for ix in text.split()]
    specN = dict()
    specN['ap']    = int(x[0])
    specN['beam']  = int(x[1])
    specN['dtype'] = int(x[2])
    specN['w1']    = x[3]
    specN['dw']    = x[4]
    specN['nw']    = int(x[5])
    specN['z']     = x[6]
    specN['aplow'] = x[7]
    specN['aphigh']= x[8]
    if specN['dtype']>=2:
        i = 1
        while i>0:
            j = i*9
            try:
                specN['wt_{0:d}'.format(i)]    = x[j]
                specN['w0_{0:d}'.format(i)]    = x[j+1]
                specN['ftype_{0:d}'.format(i)] = int(x[j+2])
                specN['order_{0:d}'.format(i)] = int(x[j+3])
                specN['pmin_{0:d}'.format(i)]  = x[j+4]
                specN['pmax_{0:d}'.format(i)]  = x[j+5]
                specN['coeffs_{0:d}'.format(i)]= x[j+6:j+6+int(x[j+3])]
                i+=1
            except (KeyError,IndexError):
                break
        specN['nfuncs'] = i-1
    return specN


def cheby(i,n):
    if i==1:
        return 1.0
    elif i==2:
        return n
    else:
        return 2 * n * cheby(i-1,n) - cheby(i-2,n)
