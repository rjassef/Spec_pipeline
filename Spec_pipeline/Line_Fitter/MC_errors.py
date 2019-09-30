import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import multiprocessing as mp
from functools import partial
import copy

from . import fit_chi2 as fit

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


def get_error(x,xbf):
    cf = 68.3 #%
    xlow = xbf - np.percentile(x[x<xbf], 100.-cf)*xbf.unit
    xhig = np.percentile(x[x>xbf], cf)*xbf.unit - xbf
    return xlow, xhig


def fit_func(spec, line_fitter, flam_resamp):
    
    nrep_x = len(flam_resamp)
    output = np.zeros((nrep_x,5))
    for i in range(nrep_x):
        new_spec = copy.deepcopy(spec)
        new_spec.flam = flam_resamp[i]
        lam_cenx, flam_line_cenx, sigma_vx, ax, bx, flam_mod = fit.fit(
            new_spec, line_fitter, 
            line_fitter.lam_cen_fit, line_fitter.sigma_v_fit)
        output[i,:] = [lam_cenx.to(u.AA).value, 
                       flam_line_cenx.to(u.erg/u.s/u.cm**2/u.AA).value, 
                       sigma_vx.to(u.km/u.s).value,
                       ax.to(u.erg/u.s/u.cm**2/u.AA**2).value,
                       bx.to(u.erg/u.s/u.cm**2/u.AA).value]
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
    lam_cen       = Output[:,0] * u.AA
    flam_line_cen = Output[:,1] * u.erg/u.s/u.cm**2/u.AA
    sigma_v       = Output[:,2] * u.km/u.s
    
    if save_chain is not None:
        #np.savetxt(save_chain,np.array([lam_cen.value,flam_line_cen.value,
        #                                sigma_v.value]).T)
        np.savetxt(save_chain,Output)

    lam_cen_low, lam_cen_hig = get_error(lam_cen, line_fitter.lam_cen_fit)
    flam_line_cen_low, \
        flam_line_cen_hig = get_error(flam_line_cen, 
                                      line_fitter.flam_line_cen_fit)
    sigma_v_low, sigma_v_hig = get_error(sigma_v, 
                                         line_fitter.sigma_v_fit)

    return [lam_cen_low, lam_cen_hig, flam_line_cen_low, flam_line_cen_hig, sigma_v_low, sigma_v_hig]

