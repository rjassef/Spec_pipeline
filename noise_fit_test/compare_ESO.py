#!/usr/bin/env python 

#import warnings
#warnings.simplefilter("ignore")

import numpy as np
import sys

from Fake_ESO_Spec import Fake_ESO_Spec
from Spec_pipeline.Spec_Reader.Spec import Spec
import general_compare as gc

if len(sys.argv)!=2:
    print("Correct use: python",sys.argv[0],"moon_degs")
    sys.exit()
moon_template = int(float(sys.argv[1]))

z_values = np.arange(0.5,3.0+0.1,0.5)
gmags    = np.arange(20.0,25.0+0.1,1.0)
moons = [0,90,180]

for z in z_values:
    for gmag in gmags:

        #Create the object for the sky template.
        template_spec = Fake_ESO_Spec(gmag=gmag,z=z,moon=moon_template)

        for moon in moons:

            print("z={0:.1f} g={1:.1f} moon={2:d}".format(z,gmag,moon))

            #Create the fake object.
            fake_spec = Fake_ESO_Spec(gmag=gmag,z=z,moon=moon)

            #Now, starting from the Spec object, recreate an ESO spectrum to
            #find the best-fit spectrum.
            eso = Spec('fake_qso',fake_spec.z)
            gc.load_new_object(eso,fake_spec,
                               flam_sky_use=template_spec.flam_sky)

            #Run the MC.
            nrep = 100
            MC_filename = "MCfiles/err_ESO.z{0:.1f}.g{1:.1f}.m{2:d}.mt{3:d}.dat".format(z,gmag,moon,moon_template)
            gc.run_MC(eso,nrep,fake_spec,MC_filename)

            #Make the plot
            plot_name = "plots/err_ESO.z{0:.1f}.g{1:.1f}.m{2:d}.mt{3:d}.png".format(z,gmag,moon,moon_template)
            gc.plot_MC(MC_filename,fake_spec,plot_name=plot_name)
            print()


