#!/usr/bin/env python

import numpy as np
import sys
import re
import matplotlib.pyplot as plt

from Spec_pipeline.Spec_Reader.iraf_spectrum1d import read_fits_spectrum1d

for i in range(1,len(sys.argv)):
    spec = read_fits_spectrum1d(sys.argv[i])

    #Y axis is in percentages. Change into fraction.
    spec[0].data *= 0.01

    fname_out = re.sub(".fits",".txt",sys.argv[i])
    fname_out = re.sub("^s","S",fname_out)
    np.savetxt(fname_out,np.array([spec[0].dispersion.value,spec[0].data]).T)

    plt.plot(spec[0].dispersion,spec[0].data,'-b')

    title = re.sub(".fits","",sys.argv[i])
    title = re.sub("_"," ",title)
    title = re.sub("sens","",title)
    plt.title(title)

    fig_name = re.sub(".fits",".png",sys.argv[i])
    fig_name = re.sub("^s","S",fig_name)
    plt.savefig(fig_name)
    plt.close()
