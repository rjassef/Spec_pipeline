#!/usr/bin/env python

import warnings
warnings.simplefilter("ignore")

import numpy as np
import astropy.units as u
import multiprocessing as mp
import copy
import matplotlib.pyplot as plt

from Spec_pipeline.Spec_Reader.read_spec import read_spec
from Spec_pipeline.Line_Fitter.multi_line import Multi_Line_fit

cat = open("data.txt")
x = cat.readline().split()

spec = read_spec(x[0],float(x[1]),x[2],x[4:],grname=x[3])

line_fit = Multi_Line_fit("CIV",spec=spec)

line_fit.run_fit()
line_fit.run_MC(100,Ncpu=7)
