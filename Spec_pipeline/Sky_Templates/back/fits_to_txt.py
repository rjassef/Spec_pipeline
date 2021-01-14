#!/usr/bin/env python

import numpy as np
import astropy.units as u
#import os
#import sys

from Spec_pipeline.Spec_Reader.iraf_spectrum1d import read_fits_spectrum1d

def wspectext(in_file, out_file):
    spec = read_fits_spectrum1d(in_file)
    np.savetxt(out_file, np.array([spec[0].dispersion.value,spec[0].data]).T)
    return

#wspectext("sky_keck_b.w.fits", "sky_keck_b.w.txt")
#wspectext("sky_keck_r.w.fits", "sky_keck_r.w.txt")
wspectext("sky_palomar_b.w.fits", "sky_palomar_D55_b.w.txt")
wspectext("sky_palomar_r.w.fits", "sky_palomar_D55_r.w.txt")
wspectext("sky_palomar_D68_b.w.fits", "sky_palomar_D68_b.w.txt")
wspectext("sky_palomar_D68_r.w.fits", "sky_palomar_D68_r.w.txt")

# if os.environ['CONDA_DEFAULT_ENV']!='iraf27':
#     print("Script need the IRAF environment to run.")
#     print("Please run 'conda activate iraf27' before running script")
#     sys.exit(0)
#
# from pyraf import iraf
#
# #Load the needed modules.
# iraf.noao(Stderr=1)
# iraf.onedspec(Stderr=1)
#
# #Run wspectext.
# iraf.wspectext.header = "no"
#
# iraf.wspectext("sky_keck_b.w.fits", "sky_keck_b.w.txt")
# iraf.wspectext("sky_keck_r.w.fits", "sky_keck_r.w.txt")
# iraf.wspectext("sky_palomar_b.w.fits", "sky_palomar_b.w.txt")
# iraf.wspectext("sky_palomar_r.w.fits", "sky_palomar_r.w.txt")
