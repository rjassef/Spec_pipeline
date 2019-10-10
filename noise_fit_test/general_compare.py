import numpy as np
import matplotlib.pyplot as plt

import sys

def get_error(x,xbf,cf=68.3):
    xlow = np.percentile(x[x<xbf], 100.-cf)
    xhig = np.percentile(x[x>xbf], cf)
    return xlow, xhig

def load_new_object(obj_new, obj_old, flam_sky_use=None):

    #Load the basic information.
    obj_new.RT    = np.copy(obj_old.RT)*obj_old.RT.unit
    obj_new.dlam  = np.copy(obj_old.dlam)*obj_old.dlam.unit
    obj_new.texp  = np.copy(obj_old.texp)*obj_old.texp.unit
    obj_new.RON   = np.copy(obj_old.RON)

    #Copy the wavelength range, sky spectra and sensisitivity curve.
    obj_new.lam_obs  = np.copy(obj_old.lam_obs)*obj_old.lam_obs.unit
    obj_new.sens     = np.copy(obj_old.sens)*obj_old.sens.unit
    
    if flam_sky_use is None:
        obj_new.flam_sky = np.copy(obj_old.flam_sky)*obj_old.flam_sky.unit
    else:
        obj_new.flam_sky = flam_sky_use

    #Set the switch to not save the error spectrum.
    obj_new.save_err = False
    obj_new.spec_err_name = "err.fake_qso.dat"

    return

######

def run_MC(lris_obj, nrep, fake_obj, MC_filename):

    cato = open(MC_filename,"w")
    cato.write("{0:d} {1:d}\n".format(nrep,len(lris_obj.lam_obs)))
    np.savetxt(cato,lris_obj.lam_obs.value)
    for k in range(nrep):
        sys.stderr.write("\r{0:3d}".format(k))
        lris_obj.flam = fake_obj.gen_obs_spec()
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
                                        cf=68.3)

    plt.fill_between(lam,flam_err_min,flam_err_max,color='xkcd:cyan',alpha=1.0)
    plt.fill_between(lam,flam_err_low,flam_err_hig,color='xkcd:red',alpha=1.0)
    plt.plot(fake_obj.lam_obs, fake_obj.flam_err, '-k')
    if plot_name is None:
        plt.show()
    else:
        plt.savefig(plot_name,dpi=300)
    plt.close()

    return


#######
