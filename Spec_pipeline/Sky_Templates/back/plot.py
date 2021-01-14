#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv)!=2:
    print("Correct use: python",sys.argv[0],"instrument")
    sys.exit()

#First assume it is dual spec.
try:
    sky_b = np.loadtxt("template_sky_{0:s}_b.dat".format(sys.argv[1]))
    sky_r = np.loadtxt("template_sky_{0:s}_r.dat".format(sys.argv[1]))
    plt.plot(sky_b[:,0],sky_b[:,1],'b-')
    plt.plot(sky_r[:,0],sky_r[:,1],'r-')
except IOError:
    #If not, try single spec.
    try:
        sky = np.loadtxt("template_sky_{0:s}.dat".format(sys.argv[1]))
        plt.plot(sky[:,0],sky[:,1],'g-')
    except IOError:
        print("Cannot open template for instrument {0:s}".format(sys.argv[1]))
        sys.exit()

plt.show()
