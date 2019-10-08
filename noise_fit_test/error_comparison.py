#!/usr/bin/env python 

import warnings
warnings.simplefilter("ignore")

import numpy as np
import matplotlib.pyplot as plt

from Fake_LRIS_Spec import Fake_LRIS_Spec
from Spec_pipeline.Spec_Reader.Spec import Spec

#####

#show_plots = True
show_plots = False

#####


#Create the fake object.
fake_spec_b = Fake_LRIS_Spec('blue',Vmag=24.0,z=2.5)
#fake_spec_r = Fake_LRIS_Spec('red')

if show_plots:
    plt.plot(fake_spec_b.lam_obs, fake_spec_b.gen_obs_spec(),
             linestyle='solid',color='xkcd:blue')
    plt.plot(fake_spec_b.lam_obs, fake_spec_b.flam,
             linestyle='solid',color='xkcd:grey')

    plt.plot(fake_spec_r.lam_obs, fake_spec_r.gen_obs_spec(),
             linestyle='solid',color='xkcd:red')
    plt.plot(fake_spec_r.lam_obs, fake_spec_r.flam,
             linestyle='solid',color='xkcd:grey')

    plt.show()


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

#Now, assign it a flam and determine the error spectrum, and compare
#it to the original one.
lris_b.flam = fake_spec_b.gen_obs_spec()
lris_b.flam_err
print(lris_b.K1,lris_b.K2)

plt.plot(lris_b.lam_obs,lris_b.flam_err,'-b')
plt.plot(fake_spec_b.lam_obs, fake_spec_b.flam_err, '--r')
plt.show()
