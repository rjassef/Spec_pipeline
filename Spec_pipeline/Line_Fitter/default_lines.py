import numpy as np
import astropy.units as u
from astropy.constants import c
import os

from .line_class import Line_fit

class Default_Line_fit(Line_fit):
    
    def __init__(self,_line_name):
        
        #Search the list for the line in question.
        cat = open(os.environ['SPEC_PIPE_LOC']+\
                   "/Spec_pipeline/Line_Fitter/lines.txt","r")
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

    #Our default line fitting class will have a Gaussian Profile.
    def flam_line_model(self, lam, lam_cen=None, flam_line_cen=None,
                        sigma_v=None):
        if lam_cen is None: lam_cen = self.lam_cen_fit
        if flam_line_cen is None: flam_line_cen = self.flam_line_cen_fit
        if sigma_v is None: sigma_v = self.sigma_v_fit
        v = c*(lam/lam_cen-1.)
        return flam_line_cen * np.exp(-0.5*(v/sigma_v)**2)

    #Continuum will be a line.
    def flam_cont_model(self,lam,a=None,b=None):
        if a is None: a = self.a
        if b is None: b = self.b
        return a*lam+b

