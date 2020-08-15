#!/usr/bin/env python

import numpy as np
import os
import sys

if os.environ['CONDA_DEFAULT_ENV']!='iraf27':
    print("Script need the IRAF environment to run.")
    print("Please run 'conda activate iraf27' before running script")
    sys.exit(0)

from pyraf import iraf

#Load the needed modules.
iraf.noao(Stderr=1)
iraf.onedspec(Stderr=1)

#Run wspectext.
iraf.wspectext.header = "no"

iraf.wspectext("sky_keck_b.w.fits", "sky_keck_b.w.txt")
iraf.wspectext("sky_keck_r.w.fits", "sky_keck_r.w.txt")
iraf.wspectext("sky_palomar_b.w.fits", "sky_palomar_b.w.txt")
iraf.wspectext("sky_palomar_r.w.fits", "sky_palomar_r.w.txt")


