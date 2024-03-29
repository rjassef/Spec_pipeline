import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import multiprocessing as mp
from functools import partial
import copy
import psutil

from . import new_fit_general as fit

####

def my_hist(x,xbf,xlow,xhig,ax,nbins=None,cum=False):


    if cum is False:
        ax.hist(x,nbins)
    else:
        ax.hist(x,nbins,cumulative=True,histtype='step',density=True)
        ax.hist(x,nbins,cumulative=-1  ,histtype='step',density=True)
    ymin,ymax = ax.ylim()
    x1  = xbf-xlow
    x2  = xbf+xhig
    ax.plot([xbf,xbf],[ymin,ymax],'k-')
    ax.plot([x1 ,x1 ],[ymin,ymax],'b--')
    ax.plot([x2 ,x2 ],[ymin,ymax],'b--')

def plot_MC_chain(chain,bf_par, err_low, err_hig):

    print("{0:.2e} -{1:.2f} +{2:.2f}".format(
        bf_par,err_low,err_hig))

    #First non-cumulative
    f, ax = plt.subplots(2)
    my_hist(chain,bf,err_low,err_hig,ax[0],nbins=50,cum=False)

    #The the cumulative below.
    my_hist(chain,bf,err_low,err_hig,ax[1],nbins=50,cum=True)

    plt.show()

####


def get_error(x,xbf,cf=68.3):
    try:
        xlow = xbf - np.percentile(x[x<xbf], 100.-cf)
    except IndexError:
        try:
            xlow = 0.*xbf.unit
        except AttributeError:
            xlow = 0.
    try:
        xhig = np.percentile(x[x>xbf], cf) - xbf
    except IndexError:
        try:
            xhig = 0.*xbf.unit
        except AttributeError:
            xhig = 0.
    return xlow, xhig


def fit_func(spec, line_fitter, flam_resamp):

    nrep_x = len(flam_resamp)
    output = np.zeros((nrep_x,line_fitter.npar_fit))
    for i in range(nrep_x):
        new_spec = copy.deepcopy(spec)
        new_spec.flam = flam_resamp[i]
        x0 = line_fitter.xopt
        xopt = fit.fit(new_spec, line_fitter, x0=x0)
        output[i,:] = xopt
        del new_spec

    return output


def MC_errors(nrep, spec, line_fitter,
              save_chain=None,Ncpu=None):


    #Create the resampled spectra.
    flam_wide     = np.tile(spec.flam,(nrep,1))
    flam_err_wide = np.tile(spec.flam_err,(nrep,1))
    flam_resamp = np.random.normal(flam_wide,flam_err_wide)*spec.flam.unit

    #Star the multiprocessing.
    if Ncpu is None:
        Ncpu = mp.cpu_count()
    Pool = mp.Pool(Ncpu)
    func = partial(fit_func, spec, line_fitter)

    #Produce the data chunks.
    flam_resamp_split = np.array_split(flam_resamp,Ncpu)
    Output = Pool.map(func,flam_resamp_split)
    Output = np.vstack(Output)
    Pool.close()

    line_fitter.MC_chain = Output

    line_fitter.parse_chain_output(Output)

    if save_chain is not None:
        np.savetxt(save_chain,Output)

    return
