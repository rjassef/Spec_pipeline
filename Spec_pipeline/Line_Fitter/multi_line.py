#Class to jointly fit N emission lines. 

import numpy as np
import astropy.units as u
from astropy.constants import c
import os
import re

from .line_class import Line_fit
from .MC_errors_general import get_error

def comb(n,k):
    c = np.math.factorial(n)/(np.math.factorial(n-k)*np.math.factorial(k))
    return np.int32(c)

class Multi_Line_fit(Line_fit):

    def __init__(self,_line_name):

        #Search the list for the line in question.
        cat = open(os.environ['SPEC_PIPE_LOC']+\
                   "/Spec_pipeline/Line_Fitter/multi_lines.txt","r")
        for line in cat:
            if not line.strip():
                continue
            x = line.split()
            if x[0]==_line_name:
                break
        cat.close()
        x[1:] = [float(ix) for ix in x[1:]]
        
        #Parse the data
        self.nlines = int(x[1])
        self.line_center = np.array(x[2:2+self.nlines*1])*u.AA
        k = 2+self.nlines
        if self.nlines>1:
            self.joint_dv = np.zeros((self.nlines,self.nlines),dtype=np.int32)
            for i in range(self.nlines):
                j1 = np.sum(np.arange(self.nlines-i,self.nlines))+k
                j2 = np.sum(np.arange(self.nlines-(i+1),self.nlines))+k
                self.joint_dv[i][i+1:] = x[j1:j2]
            k += comb(self.nlines,2)
            self.joint_sigma = np.zeros((self.nlines,self.nlines),
                                        dtype=np.int32)
            for i in range(self.nlines):
                j1 = np.sum(np.arange(self.nlines-i,self.nlines))+k
                j2 = np.sum(np.arange(self.nlines-(i+1),self.nlines))+k
                self.joint_sigma[i][i+1:] = x[j1:j2]
            k += comb(self.nlines,2)
            self.fixed_ratio= np.zeros((self.nlines,self.nlines))
            for i in range(self.nlines):
                j1 = np.sum(np.arange(self.nlines-i,self.nlines))+k
                j2 = np.sum(np.arange(self.nlines-(i+1),self.nlines))+k
                self.fixed_ratio[i][i+1:] = x[j1:j2]
        self.line_velocity_region = x[k+1]*u.km/u.s
        n_creg = int((len(x)-(k+2))/2)
        self.continuum_regions = np.zeros((n_creg,2))
        for i in range(n_creg):
            self.continuum_regions[i][0] = x[i*2+k+2]
            self.continuum_regions[i][1] = x[i*2+k+2+1]
        self.continuum_regions = self.continuum_regions*u.AA
        
