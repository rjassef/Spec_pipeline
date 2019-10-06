#!/usr/bin/env python 

import warnings
warnings.simplefilter("ignore")

import numpy as np
import astropy.units as u
import multiprocessing as mp
import copy
import matplotlib.pyplot as plt
from astropy.constants import c
from scipy.signal import savgol_filter

from Spec_pipeline.Spec_Reader.read_spec import read_spec
from Spec_pipeline.Line_Fitter.multi_line import Multi_Line_fit


em_lines = ["LyA_NVred", "CIV", "CIII", "MgII", "OII", "Hb", "OIII"]

cato = open("visualization.summary.txt","w")

cat = open("data.txt")
for line in cat:
    
    if line[0]=="#":
        continue

    x = line.split()

    #First, plot the spectrum.
    spec_b = read_spec(x[0],float(x[1]),x[2],x[3:],blue=True)
    spec_r = read_spec(x[0],float(x[1]),x[2],x[3:],red=True)
    spec_b.flam[spec_b.flam<0] = 0.
    spec_r.flam[spec_r.flam<0] = 0.
    plt.plot(spec_b.lam_rest,
             savgol_filter(spec_b.flam,5,2),
             linestyle='solid',color='xkcd:grey',linewidth=0.5)    
    plt.plot(spec_r.lam_rest,
             savgol_filter(spec_r.flam,5,2),
             linestyle='solid',color='xkcd:grey',linewidth=0.5)

    #Now, fit the line in turn.
    for k,em_line in enumerate(em_lines):
        print(x[0],em_line)
        line_fit = Multi_Line_fit(em_line)
        spec = read_spec(x[0],float(x[1]),x[2],x[3:],
                         line_center=line_fit.line_center[0])
        line_fit.run_fit(spec)

        #Check if fit was done. If not, move on to next emission line.
        if line_fit.xopt_line is None:
            continue
        
        #Run the Ftest
        line_fit.run_Ftest(spec)

        #Write the results.
        cato.write("{0:15s} {1:10s} {2:8.3f}".format(x[0],em_line,spec.zspec))
        cato.write("{0:7.3f}".format(1.-line_fit.p))
        if 1.-line_fit.p > 0.5:
            for i in range(line_fit.nlines):
                cato.write("{0:10.1f}".format(line_fit.dv_fit[i].value))
                cato.write("{0:10.1f}".format(line_fit.FWHM_v[i].value))
                zline = spec.zspec+line_fit.dv_fit[i]/c
                cato.write("{0:8.3f}".format(zline))
        else:
            for i in range(line_fit.nlines):
                cato.write("{0:^10s} {0:^10s} {0:^8s} {0:^8s}".format("--"))
        cato.write("\n")

        #If line is well detected, plot. 
        if line_fit.xopt_line is not None and 1.-line_fit.p>0.5:
            i_line = line_fit.get_i_line(spec)
            lam_use = spec.lam_rest[i_line]
            flam_model = line_fit.flam_model(lam_use)
            plt.plot(lam_use,flam_model)
            
    #plt.show()
    fig_name = x[0]+".png"
    plt.savefig(fig_name,dpi=300)
    plt.close()
    
