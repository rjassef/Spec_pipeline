#!/usr/bin/env python 

import warnings
warnings.simplefilter("ignore")

import numpy as np

from Fake_LRIS_Spec import Fake_LRIS_Spec
from Spec_pipeline.Spec_Reader.Spec import Spec

import general_compare as gc

z_values = np.arange(0.5,3.0+0.1,0.5)
rmags    = np.arange(20.0,25.0+0.1,1.0)

for z in z_values:
    for rmag in rmags:

        print("z={0:.1f} r={1:.1f}".format(z,rmag))

        #Create the fake object.
        fake_spec_r = Fake_LRIS_Spec('red',rmag=rmag,z=z)

        #Now, starting from the Spec object, recreate an LRIS spectrum to
        #find the best-fit spectrum.
        lris_r = Spec('fake_qso_red',fake_spec_r.z)
        gc.load_new_object(lris_r,fake_spec_r)

        #Run the MC.
        nrep = 100
        MC_filename = "MCfiles/err_r.z{0:.1f}.r{1:.1f}.dat".format(z,rmag)
        gc.run_MC(lris_r,nrep,fake_spec_r,MC_filename)

        #Make the plot
        plot_name = "plots/err_r.z{0:.1f}.r{1:.1f}.png".format(z,rmag)
        gc.plot_MC(MC_filename,fake_spec_r,plot_name=plot_name)
        print()


