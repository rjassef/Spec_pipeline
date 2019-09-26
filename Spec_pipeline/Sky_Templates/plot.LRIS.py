#!/usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt

sky_b = np.loadtxt("template_sky_LRIS_b.dat")
sky_r = np.loadtxt("template_sky_LRIS_r.dat")

plt.plot(sky_b[:,0],sky_b[:,1],'b-')
plt.plot(sky_r[:,0],sky_r[:,1],'r-')
plt.show()

