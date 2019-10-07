#!/usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt

from Fake_LRIS_Spec import Fake_LRIS_Spec

#Create the fake object.
fake_spec_b = Fake_LRIS_Spec('blue')
fake_spec_r = Fake_LRIS_Spec('red')

plt.plot(fake_spec_b.lam_obs, fake_spec_b.gen_obs_spec(),
         linestyle='solid',color='xkcd:blue')
plt.plot(fake_spec_b.lam_obs, fake_spec_b.flam,
         linestyle='solid',color='xkcd:grey')

plt.plot(fake_spec_r.lam_obs, fake_spec_r.gen_obs_spec(),
         linestyle='solid',color='xkcd:red')
plt.plot(fake_spec_r.lam_obs, fake_spec_r.flam,
         linestyle='solid',color='xkcd:grey')

plt.show()
