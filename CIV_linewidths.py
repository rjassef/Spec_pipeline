#!/usr/bin/env python 

import warnings
warnings.simplefilter("ignore")

import numpy as np
import astropy.units as u
import multiprocessing as mp
import copy

import sys
sys.path.append("Spec_Reader")
from read_spec import read_spec
sys.path.append("Line_Fitter")
from default_lines import Default_Line_fit

######

verbose = True
#verbose = False

Ncpu = mp.cpu_count()-2
#nrep = 1000
#nrep = 100
nrep = 10

estimate_errors = True
#estimate_errors = False

######

#Define the continuum regions.
continuum_regions = [[1425., 1470.],[1680.,1705.]]*u.AA

#Define the velocity region to consider for the line fit. 
line_velocity_region = 10000.*u.km/u.s

#Line center
line_center = 1549.*u.AA

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

    #Read the line and set up the correct object.
    x = line.split()
    spec = read_spec(x[0],float(x[1]),x[2],x[3:],line_center=line_center)
    if verbose is True:
        print
        print(spec.name)

    #Set the fitter object.
    civ_fit = Default_Line_fit("CIV")

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

    exit()

cato.close()
cat.close()
