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

def load_new_object(obj_new, obj_old):

    #Load the basic information.
    obj_new.RT    = np.copy(obj_old.RT)*obj_old.RT.unit
    obj_new.dlam  = np.copy(obj_old.dlam)*obj_old.dlam.unit
    obj_new.texp  = np.copy(obj_old.texp)*obj_old.texp.unit
    obj_new.RON   = np.copy(obj_old.RON)

    #Copy the wavelength range, sky spectra and sensisitivity curve.
    obj_new.lam_obs  = np.copy(obj_old.lam_obs)*obj_old.lam_obs.unit
    obj_new.flam_sky = np.copy(obj_old.flam_sky)*obj_old.flam_sky.unit
    obj_new.sens     = np.copy(obj_old.sens)*obj_old.sens.unit

    #Set the switch to not save the error spectrum.
    obj_new.save_err = False
    obj_new.spec_err_name = "err.fake_qso_b.dat"

    return

######

def run_MC(lris_obj, nrep, MC_filename):

    cato = open(MC_filename,"w")
    cato.write("{0:d} {1:d}\n".format(nrep,len(lris_obj.lam_obs)))
    np.savetxt(cato,lris_obj.lam_obs.value)
    for k in range(nrep):
        sys.stderr.write("\r{0:3d}".format(k))
        lris_obj.flam = fake_spec_b.gen_obs_spec()
        np.savetxt(cato,lris_obj.flam_err.value)
        del lris_obj._flam_err
    print()
    cato.close()
    return

#######

def plot_MC(MC_filename,fake_obj,plot_name=None):

    cat = open(MC_filename)
    nrep, nlam = [int(float(ix)) for ix in cat.readline().split()]
    cat.close()

    data = np.loadtxt(MC_filename,skiprows=1)
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
        flam_err_low[k], \
            flam_err_hig[k] = get_error(flam_err_use, 
                                        np.median(flam_err_use),
                                        #fake_spec_b.flam_err[k].value,
                                        cf=68.3)

    plt.fill_between(lam,flam_err_min,flam_err_max,color='xkcd:cyan',alpha=1.0)
    plt.fill_between(lam,flam_err_low,flam_err_hig,color='xkcd:red',alpha=1.0)
    plt.plot(fake_spec_b.lam_obs, fake_spec_b.flam_err, '-k')
    if plot_name is None:
        plt.show()
    else:
        plt.savefig(plot_name,dpi=300)
    plt.close()

    return


#######

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
        load_new_object(lris_b,fake_spec_b)

        #Run the MC.
        nrep = 100
        MC_filename = "err_b.z{0:.1f}.V{1:.1f}.dat".format(z,Vmag)
        run_MC(lris_b,nrep,MC_filename)

        #Make the plot
        plot_name = "err_b.z{0:.1f}.V{1:.1f}.png".format(z,Vmag)
        plot_MC(MC_filename,fake_spec_b,plot_name=plot_name)
        print()


