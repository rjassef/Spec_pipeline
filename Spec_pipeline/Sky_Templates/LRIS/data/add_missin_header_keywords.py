#!/usr/bin/env python

import numpy as np
from astropy.io import fits

#Add CCDGEOM to the blue side if missing.
fname = "lris_12apr_sky_b.f.fits"
h = fits.open(fname)
if 'CCDGEOM' not in h[0].header:
    h[0].header['CCDGEOM'] = 'e2v (Marconi) CCD44-82'
h.writeto(fname,overwrite=True)
