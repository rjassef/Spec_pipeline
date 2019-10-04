import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import multiprocessing as mp
from functools import partial
import copy

from . import fit_general as fit

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
    xlow = xbf - np.percentile(x[x<xbf], 100.-cf)*xbf.unit
    xhig = np.percentile(x[x>xbf], cf)*xbf.unit - xbf
    return xlow, xhig


def fit_func(spec, line_fitter, flam_resamp):
    
    nrep_x = len(flam_resamp)
    output = np.zeros((nrep_x,line_fitter.npar_fit))
    for i in range(nrep_x):
        new_spec = copy.deepcopy(spec)
        new_spec.flam = flam_resamp[i]
        x0_cont = line_fitter.xopt_cont
        x0_line = line_fitter.xopt_line
        xopt_line, xopt_cont = fit.fit(new_spec, line_fitter,
                                       x0_cont=x0_cont,
                                       x0_line=x0_line)
        output[i,:] = np.concatenate((xopt_line, xopt_cont))

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

    line_fitter.parse_chain_output(Output)
    
    if save_chain is not None:
        np.savetxt(save_chain,Output)

    line_fitter.parse_chain_output(Output)

    return

