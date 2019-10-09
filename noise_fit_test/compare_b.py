#!/usr/bin/env python 

import warnings
warnings.simplefilter("ignore")

import numpy as np

from Fake_LRIS_Spec import Fake_LRIS_Spec
from Spec_pipeline.Spec_Reader.Spec import Spec

import general_compare as gc

z_values = np.arange(0.5,3.0+0.1,0.5)
Vmags    = np.arange(20.0,25.0+0.1,1.0)

for z in z_values:
    for Vmag in Vmags:

        print("z={0:.1f} V={1:.1f}".format(z,Vmag))

        #Create the fake object.
        fake_spec_b = Fake_LRIS_Spec('blue',Vmag=Vmag,z=z)

        #Now, starting from the Spec object, recreate an LRIS spectrum to
        #find the best-fit spectrum.
        lris_b = Spec('fake_qso_blue',fake_spec_b.z)
        gc.load_new_object(lris_b,fake_spec_b)

        #Run the MC.
        nrep = 100
        MC_filename = "MCfiles/err_b.z{0:.1f}.V{1:.1f}.dat".format(z,Vmag)
        gc.run_MC(lris_b,nrep,fake_spec_b,MC_filename)

        #Make the plot
        plot_name = "plots/err_b.z{0:.1f}.V{1:.1f}.png".format(z,Vmag)
        gc.plot_MC(MC_filename,fake_spec_b,plot_name=plot_name)
        print()


