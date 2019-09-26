import numpy as np
import astropy.units as u
import os

from .line_class import Line_fit

class Default_Line_fit(Line_fit):
    
    def __init__(self,_line_name):
        
        #Search the list for the line in question.
        cat = open(os.environ['SPEC_PIPE_LOC']+"/Line_Fitter/lines.txt","r")
        for line in cat:
            x = line.split()
            if x[0]==_line_name:
                break
        cat.close()

        #Parse the data
        x[1:] = [float(ix) for ix in x[1:]]
        _line_name = x[0]
        _line_center = x[1]*u.AA
        _line_velocity_region = x[2]*u.km/u.s
        _continuum_regions = np.zeros((int((len(x)-3)/2),2))
        for i in range((int((len(x)-3)/2))):
            _continuum_regions[i][0] = x[i*2+3]
            _continuum_regions[i][1] = x[i*2+4]
        _continuum_regions = _continuum_regions*u.AA

        #Initialize the class
        super(Default_Line_fit,self).__init__(_line_name,_line_center,
                                              _line_velocity_region,
                                              _continuum_regions)


