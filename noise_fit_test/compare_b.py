#!/usr/bin/env python 

import warnings
warnings.simplefilter("ignore")

import numpy as np
import matplotlib.pyplot as plt

from Fake_LRIS_Spec import Fake_LRIS_Spec
from Spec_pipeline.Spec_Reader.Spec import Spec

import sys

def get_error(x,xbf,cf=68.3):
    xlow = np.percentile(x[x<xbf], 100.-cf)
    xhig = np.percentile(x[x>xbf], cf)
    return xlow, xhig



#Create the fake object.
fake_spec_b = Fake_LRIS_Spec('blue',Vmag=24.0,z=2.5)

#Now, starting from the Spec object, recreate an LRIS spectrum to
#find the best-fit spectrum.
lris_b = Spec('fake_qso_blue',fake_spec_b.z)

#Load the basic information.
lris_b.RT    = np.copy(fake_spec_b.RT)*fake_spec_b.RT.unit
lris_b.dlam  = np.copy(fake_spec_b.dlam)*fake_spec_b.dlam.unit
lris_b.texp  = np.copy(fake_spec_b.texp)*fake_spec_b.texp.unit
lris_b.RON   = np.copy(fake_spec_b.RON)

#Copy the wavelength range, sky spectra and sensisitivity curve.
lris_b.lam_obs  = np.copy(fake_spec_b.lam_obs)*fake_spec_b.lam_obs.unit
lris_b.flam_sky = np.copy(fake_spec_b.flam_sky)*fake_spec_b.flam_sky.unit
lris_b.sens     = np.copy(fake_spec_b.sens)*fake_spec_b.sens.unit

#Set the switch to not save the error spectrum.
lris_b.save_err = False
lris_b.spec_err_name = "err.fake_qso_b.dat"

cato = open("err_b_comp.dat","w")
nrep = 1000
cato.write("{0:d} {1:d}\n".format(nrep,len(lris_b.lam_obs)))
np.savetxt(cato,lris_b.lam_obs.value)
for k in range(nrep):
    sys.stderr.write("\r{0:3d}".format(k))
    lris_b.flam = fake_spec_b.gen_obs_spec()
    np.savetxt(cato,lris_b.flam_err.value)
    del lris_b._flam_err
print()
cato.close()

cat = open("err_b_comp.dat")
nrep, nlam = [int(float(ix)) for ix in cat.readline().split()]
cat.close()

data = np.loadtxt("err_b_comp.dat",skiprows=1)
lam = data[:nlam]
flam_err_min = data[nlam:nlam*2]
flam_err_max = data[nlam:nlam*2]
for i in range(1,nrep):
    flam_err = data[nlam*(i+1):nlam*(i+2)]
    flam_err_min = np.where(flam_err<flam_err_min,flam_err,flam_err_min)
    flam_err_max = np.where(flam_err>flam_err_max,flam_err,flam_err_max)

flam_err_low = np.zeros(nlam)
flam_err_hig = np.zeros(nlam)
for k in range(nlam):
    flam_err_use = data[nlam+k::nlam]
    flam_err_low[k], flam_err_hig[k] = get_error(flam_err_use, 
                                                 fake_spec_b.flam_err[k].value,
                                                 cf=68.3)

plt.fill_between(lam,flam_err_min,flam_err_max,color='xkcd:cyan',alpha=1.0)
plt.fill_between(lam,flam_err_low,flam_err_hig,color='xkcd:red',alpha=1.0)
plt.plot(fake_spec_b.lam_obs, fake_spec_b.flam_err, '-k')
#plt.show()
plt.savefig("compare_b.png",dpi=300)
