#!/usr/bin/env python 

import warnings
warnings.simplefilter("ignore")

import numpy as np
import astropy.units as u
import multiprocessing as mp
import copy

import sys
import os
sys.path.append(os.environ['SPEC_PIPE_LOC'])
from Spec_pipeline.Spec_Reader.read_spec import read_spec
from Spec_pipeline.Line_Fitter.default_lines import Default_Line_fit


em_lines = ["LyA", "NV", "CIV", "CIII", "Hb"]

cato = open("visualization.summary.txt","w")

cat = open("data.txt")
for line in cat:
    
    if line[0]=="#":
        continue

    x = line.split()
    for em_line in em_lines:
        print(em_line)
        line_fit = Default_Line_fit(em_line)
        spec = read_spec(x[0],float(x[1]),x[2],x[3:],
                         line_center=line_fit.line_center)
        line_fit.run_fit(spec,sigma_v_0=1000.*u.km/u.s)
        if line_fit.sigma_v_fit is not None:
            line_fit.run_Ftest(spec)
            print("{0:6.2f} {1:.3e}".format(line_fit.F, 1.-line_fit.p))
            print("{0:6.1f} {1:6.1f}".format(
                line_fit.lam_cen_fit,line_fit.FWHM_v))
        print()
    exit()

