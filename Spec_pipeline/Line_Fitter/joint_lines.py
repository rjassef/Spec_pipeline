#This is a class that will fit two lines simultaneously, allowing for
#different types of constraints. 


import numpy as np
import astropy.units as u
from astropy.constants import c
import os

from .line_class import Line_fit
from .MC_errors_general import get_error

class Default_Line_fit(Line_fit):

    def __init__(self,_line_name)

        #Values that will come from the fitting.
        self.dv1_fit        = None
        self.flam_line1_fit = None
        self.sigma_v1_fit   = None
        self.dv2_fit        = None
        self.flam_line2_fit = None
        self.sigma_v2_fit   = None
        self.a              = None
        self.b              = None

        #Values that will come from the MC

        
        #Set the default constraints for the fitting.


        #Search the list for the line in question.
        cat = open(os.environ['SPEC_PIPE_LOC']+\
                   "/Spec_pipeline/Line_Fitter/jointlines.txt","r")
        for line in cat:
            x = line.split()
            if x[0]==_line_name:
                break
        cat.close()
        
        #Parse the data
        x[1:] = [float(ix) for ix in x[1:]]
        self.line1_center = x[1]*u.AA
        self.line2_center = x[2]*u.AA
        self.joint_dv     = int(x[3])
        self.joint_sigma  = int(x[4])
        self.fixed_ratio  = int(x[5])
        self.line_velocity_region = x[2]*u.km/u.s
        self.continuum_regions = np.zeros((int((len(x)-3)/2),2))
        for i in range((int((len(x)-3)/2))):
            self.continuum_regions[i][0] = x[i*2+3]
            self.continuum_regions[i][1] = x[i*2+4]
        self.continuum_regions = _continuum_regions*u.AA

        #Initialize the class
        super(Default_Line_fit,self).__init__(_line_name)
