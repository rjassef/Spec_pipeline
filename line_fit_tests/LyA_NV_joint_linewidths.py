#!/usr/bin/env python 

import warnings
warnings.simplefilter("ignore")

import numpy as np
import astropy.units as u
import multiprocessing as mp
import copy
import matplotlib.pyplot as plt

from Spec_pipeline.Spec_Reader.read_spec import read_spec
from Spec_pipeline.Line_Fitter.joint_lines import Joint_Line_fit

######

verbose = True
#verbose = False

Ncpu = mp.cpu_count()-2
nrep = 1000
#nrep = 100
#nrep = 10

estimate_errors = True
#estimate_errors = False

#######

cato = open("LyA_NV_results.txt","w")
cato.write("{0:s}".format("obj_id"))
cato.write(" {0:s} {1:s} {2:s} ".format("dv1", "flam_line1", "FWHM_v1"))
cato.write(" {0:s} {1:s} {2:s} ".format("dv2", "flam_line2", "FWHM_v2"))
if estimate_errors is True:
    cato.write("{0:s} {1:s} ".format("dv1_low", "dv1_hig"))
    cato.write("{0:s} {1:s} ".format("flam_line1_low", "flam_line1_hig"))
    cato.write("{0:s} {1:s} ".format("FWHM_v1_low", "FWHM_v1_hig"))
    cato.write("{0:s} {1:s} ".format("dv2_low", "dv2_hig"))
    cato.write("{0:s} {1:s} ".format("flam_line2_low", "flam_line2_hig"))
    cato.write("{0:s} {1:s} ".format("FWHM_v2_low", "FWHM_v2_hig"))
cato.write("\n")

cat = open("data.txt")
for line in cat:
    
    if line[0]=="#":
        continue

    #Set the fitter object.
    line_fit = Joint_Line_fit("LyA_NV")

    #Read the line and set up the correct object.
    x = line.split()
    spec = read_spec(x[0],float(x[1]),x[2],x[3:],
                     line_center=line_fit.line1_center)
    if verbose is True:
        print
        print(spec.name)
    
    #Run the fit
    line_fit.run_fit(spec)
    cato.write("{0:s}".format(spec.name))
    cato.write(" {0:.3f} {1:.3e} {2:.3f} ".format(
        line_fit.dv1_fit.value, 
        line_fit.flam_line1_fit.value, 
        line_fit.FWHM_v1.value))
    cato.write(" {0:.3f} {1:.3e} {2:.3f} ".format(
        line_fit.dv2_fit.value, 
        line_fit.flam_line2_fit.value, 
        line_fit.FWHM_v2.value))

    #Get the fit errors.
    chain_name = None
    if estimate_errors is True:
        if verbose is True:
            print("Running MC with ",nrep,"steps in ",Ncpu,"cores...")
        chain_name = spec.name+".LyA_NV_chain.txt"
        line_fit.run_MC(spec,nrep,Ncpu=Ncpu,save_chain=chain_name)
        cato.write("{0:.3f} {1:.3f} ".format(line_fit.dv1_low.value,
                                             line_fit.dv1_hig.value))
        cato.write("{0:.3f} {1:.3f} ".format(line_fit.flam_line1_low.value,
                                             line_fit.flam_line1_hig.value))
        cato.write("{0:.3f} {1:.3f} ".format(line_fit.FWHM_v1_low.value,
                                             line_fit.FWHM_v1_hig.value))
        cato.write("{0:.3f} {1:.3f} ".format(line_fit.dv2_low.value,
                                             line_fit.dv2_hig.value))
        cato.write("{0:.3f} {1:.3f} ".format(line_fit.flam_line2_low.value,
                                             line_fit.flam_line2_hig.value))
        cato.write("{0:.3f} {1:.3f} ".format(line_fit.FWHM_v2_low.value,
                                             line_fit.FWHM_v2_hig.value))

    #Plot the fit.
    if verbose is True:
        print("Making the plot...")
    plot_fname = "{0:s}.LyA_NV_fit.png".format(spec.name)
    line_fit.plot(spec,chain=chain_name,plot_fname=plot_fname)

    cato.write("\n")

cato.close()
cat.close()
