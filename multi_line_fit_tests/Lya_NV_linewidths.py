#!/usr/bin/env python 

#import warnings
#warnings.simplefilter("ignore")

import numpy as np
import multiprocessing as mp
import copy
import matplotlib.pyplot as plt

from Spec_pipeline.Spec_Reader.read_spec import read_spec
from Spec_pipeline.Line_Fitter.multi_line import Multi_Line_fit

######

verbose = True
#verbose = False

Ncpu = mp.cpu_count()-2
#nrep = 1000
#nrep = 100
nrep = 10

estimate_errors = False
#estimate_errors = True

#make_plot = False
make_plot = True

#######

cato = open("Lya_NV_results.txt","w")
cato.write("{0:s} ".format("obj_id"))

cat = open("data.txt")
for i,line in enumerate(cat):
    
    if line[0]=="#":
        continue

    #Set the fitter object.
    #lya_nv_fit = Multi_Line_fit("LyA_NV")
    lya_nv_fit = Multi_Line_fit("LyA_NVred")

    if i==0:
        cato.write(lya_nv_fit.print_fit_header())
        if estimate_errors:
            cato.write(lya_nv_fit.print_MC_header())
        cato.write("\n")
    
    #Read the line and set up the correct object.
    x = line.split()
    spec = read_spec(x[0],float(x[1]),x[2],x[4:],grname=x[3],
                     line_center=lya_nv_fit.line_center[0])
    if verbose is True:
        print
        print(spec.name)
    cato.write("{0:s} ".format(spec.name))
        
    #Run the fit
    lya_nv_fit.run_fit(spec)
    cato.write(lya_nv_fit.print_fit())

    #Get the fit errors.
    chain_name = None
    if estimate_errors is True:
        if verbose is True:
            print("Running MC with ",nrep,"steps in ",Ncpu,"cores...")
        chain_name = spec.name+".CIV_chain.txt"
        lya_nv_fit.run_MC(spec,nrep,Ncpu=Ncpu,save_chain=chain_name)
        cato.write(lya_nv_fit.print_MC())

    #Plot the fit.
    if make_plot:
        if verbose is True:
            print("Making the plot...")
        plot_fname = "{0:s}.Lya_NV_fit.png".format(spec.name)
        lya_nv_fit.plot(spec,chain=chain_name,plot_fname=plot_fname)

    cato.write("\n")

cato.close()
cat.close()
