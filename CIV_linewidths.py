#!/usr/bin/env python 

import warnings
warnings.simplefilter("ignore")

import numpy as np
import astropy.units as u
import multiprocessing as mp
import copy

from Spec_pipeline.Spec_Reader.read_spec import read_spec
from Spec_pipeline.Line_Fitter.default_lines import Default_Line_fit

######

verbose = True
#verbose = False

Ncpu = mp.cpu_count()-2
#nrep = 1000
#nrep = 100
nrep = 10

#estimate_errors = True
estimate_errors = False

#######

cato = open("results.txt","w")
cato.write("{0:s} {1:s} {2:s} {3:s} ".format(
        "obj_id","lam_cen", "flam_line_cen", "FWHM_v"))
if estimate_errors is True:
    cato.write("{0:s} {1:s} ".format("lam_cen_low", "lam_cen_hig"))
    cato.write("{0:s} {1:s} ".format("flam_line_cen_low",
                                     "flam_line_cen_hig"))
    cato.write("{0:s} {1:s} ".format("FWHM_v_low", "FWHM_v_hig"))
cato.write("\n")

cat = open("data.txt")
for line in cat:
    
    if line[0]=="#":
        continue

    #Set the fitter object.
    civ_fit = Default_Line_fit("CIV")

    #Read the line and set up the correct object.
    x = line.split()
    spec = read_spec(x[0],float(x[1]),x[2],x[3:],
                     line_center=civ_fit.line_center)
    if verbose is True:
        print
        print(spec.name)

    #Run the fit
    civ_fit.run_fit(spec)
    cato.write("{0:s} {1:.3f} {2:.3e} {3:.3f} ".format(
        spec.name, civ_fit.lam_cen_fit.value, 
        civ_fit.flam_line_cen_fit.value, 
        civ_fit.FWHM_v.value))

    #Get the fit errors.
    if estimate_errors is True:
        if verbose is True:
            print("Running MC with ",nrep,"steps in ",Ncpu,"cores...")
        civ_fit.run_MC(spec,nrep,Ncpu=Ncpu)
        cato.write("{0:.3f} {1:.3f} ".format(civ_fit.lam_cen_low.value,
                                             civ_fit.lam_cen_hig.value))
        cato.write("{0:.3e} {1:.3e} ".format(civ_fit.flam_line_cen_low.value,
                                             civ_fit.flam_line_cen_hig.value))
        cato.write("{0:.3f} {1:.3f} ".format(civ_fit.FWHM_v_low.value,
                                             civ_fit.FWHM_v_hig.value))

    #Plot the fit.
    if verbose is True:
        print("Making the plot...")
    #plot_fname = "{0:s}.CIV.eps".format(obj_id)
    civ_fit.plot(spec)

    cato.write("\n")

    #exit()

cato.close()
cat.close()
