#!/usr/bin/env python 

import warnings
warnings.simplefilter("ignore")

import numpy as np

from Fake_LRIS_Spec import Fake_LRIS_Spec
from Spec_pipeline.Spec_Reader.Spec import Spec

import general_compare as gc

z_values = np.arange(0.5,3.0+0.1,0.5)
gmags    = np.arange(20.0,25.0+0.1,1.0)

for z in z_values:
    for gmag in gmags:

        print("z={0:.1f} g={1:.1f}".format(z,gmag))

        #Create the fake object.
        fake_spec_b = Fake_LRIS_Spec('blue',gmag=gmag,z=z)

        #Now, starting from the Spec object, recreate an LRIS spectrum to
        #find the best-fit spectrum.
        lris_b = Spec('fake_qso_blue',fake_spec_b.z)
        gc.load_new_object(lris_b,fake_spec_b)

        #Run the MC.
        nrep = 100
        MC_filename = "MCfiles/err_LRIS_b.z{0:.1f}.g{1:.1f}.dat".format(z,gmag)
        gc.run_MC(lris_b,nrep,fake_spec_b,MC_filename)

        #Make the plot
        plot_name = "plots/err_LRIS_b.z{0:.1f}.g{1:.1f}.png".format(z,gmag)
        gc.plot_MC(MC_filename,fake_spec_b,plot_name=plot_name)
        print()


