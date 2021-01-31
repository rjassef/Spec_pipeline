#!/usr/bin/env python

from astropy.io import fits

missing_kws = {
    'SLITNAME': 'long_1.0',
    'DICHNAME': '560     ',
    'GRANAME' : '400/8500',
    'GRISNAME': '600/4000',
    'EXPTIME' : 599.99920654,
}

missing_blue_kws = {
}

missing_red_kws = {
    'DETECTOR': 'LRIS-R Science mosaic Mark 2: CCD  #1: 19-3  CCD #2: 19-2'
}

for chan in ["b","r"]:
    h = fits.open("bak/sky_keck_{0:s}.w.fits".format(chan))

    for kw in missing_kws.keys():
        if kw not in h[0].header:
            h[0].header[kw] = missing_kws[kw]

    if chan=='b':
        for kw in missing_blue_kws.keys():
            if kw not in h[0].header:
                h[0].header[kw] = missing_blue_kws[kw]
    elif chan=='r':
        for kw in missing_red_kws.keys():
            if kw not in h[0].header:
                h[0].header[kw] = missing_red_kws[kw]

    h.writeto("sky_keck_{0:s}.w.fits".format(chan), overwrite=True)
    h.close()

#Now, add the CCDGEOM missing to w1016p5105_b.sky.f.fits
h = fits.open("bak/w1016p5105_b.sky.f.fits")
h[0].header['CCDGEOM'] = 'e2v (Marconi) CCD44-82'
h.writeto("w1016p5105_b.sky.f.fits".format(chan), overwrite=True)
h.close()
