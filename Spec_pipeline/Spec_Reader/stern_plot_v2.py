import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import astropy.units as u
import re
import os

from ..Line_Fitter.multi_line import Multi_Line_fit

##
# Properties
#
# 1. Line list
# 2. Wavelength range of the spectra
# 3. Smooth properties
# 4. Wavelength range of the emission lines within spectra
# 5. Determine the peak of the region around the emission lines, both for the y-axis range as well as for the emission line label. 
# 6. 