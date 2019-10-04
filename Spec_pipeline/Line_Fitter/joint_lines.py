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
        self.fixed_ratio  = x[5]
        self.line_velocity_region = x[2]*u.km/u.s
        self.continuum_regions = np.zeros((int((len(x)-3)/2),2))
        for i in range((int((len(x)-3)/2))):
            self.continuum_regions[i][0] = x[i*2+3]
            self.continuum_regions[i][1] = x[i*2+4]
        self.continuum_regions = _continuum_regions*u.AA

        #Initialize the class
        super(Default_Line_fit,self).__init__(_line_name)

    ###################
    # Fit Model
    ###################
    
    #Order of the parameters is always:
    #dv1_fit
    #flam_line1_fit
    #sigma_v1_fit
    #dv2_fit (if join_dv is 0)
    #flam_line2_fir (if fixed_ratio is 0)
    #sigma_v2_fit (if joint_sigma is 0)
    #a
    #b

    #Note that line 1 is always bluer than line 2.
    
    #Complete model.
    def flam_model(self,lam,x_line=None,x_cont=None,chain_output=None):

        npar_line = 6
        if joint_dv==1   : npar_line -= 1
        if fixed_ratio>0 : npar_line -= 1
        if joint_sigma==1: npar_line -= 1
        
        if chain_output is not None:
            x_line = chain_output[:,:npar_line].T
            x_cont = chain_output[:,npar_line:].T
            
        return self.flam_cont_model(lam,x_cont)+\
            self.flam_line_model(lam,x_line)

    #A Gaussian per line.
    def flam_line_model(self,lam,x_line):

        if x_line is not None:
            dv1, flam_line1, sigma_v1,\
                dv2, flam_line2, sigma_v2 = self.line_par_parser(x_line)
        else:
            dv1 = self.dv1_fit
            dv2 = self.dv2_fit
            flam_line1 = self.flam_line1_fit
            flam_line2 = self.flam_line2_fit
            sigma_v1 = self.sigma_v1
            sigma_v2 = self.sigma_v2
            
        v1 = c*(lam/self.line1_center-1.)
        v2 = c*(lam/self.line2_center-1.)
        
        line1 = flam_line1 * np.exp(-0.5*((v1-dv1)/sigma_v1)**2)
        line2 = flam_line2 * np.exp(-0.5*((v2-dv2)/sigma_v2)**2)
        return line1 + line2

    #Continuum will be a straight line.
    def flam_cont_model(self,lam,x_cont=None):
        if x_cont is not None:
            a = x_cont[0]*u.erg/u.cm**2/u.s/u.AA**2
            b = x_cont[1]*u.erg/u.cm**2/u.s/u.AA
        else:
            a = self.a
            b = self.b
        return a*lam+b


    def line_par_paser(self,x_line):

        dv1        = x_line[0]*u.km/u.s
        flam_line1 = x_line[1]*u.erg/u.s/u.cm**2/u.AA
        sigma_v1   = x_line[2]*u.km/u.s
            
        k=3
        if joint_dv==0:
            dv2 = x_line[k]*u.km/u.s
            k+=1
        else:
            dv2 = dv1
            if fixed_ratio>0:
                flam_line2 = x_line[k]*u.erg/u.s/u.cm**2/u.AA
                k+=1
            else:
                flam_line2 = flam_line1 * fixed_ratio
            if joint_sigma==0:
                sigma_v2 = x_line[k]*u.km/u.s
                k+=1
            else:
                sigma_v2 = sigma_v1
        return dv1, flam_line1, sigma_v1, dv2, flam_line2, sigma_v2


    ###############
    # Constraints
    ###############

    #Constraints on the continuum fit parameters.
    def meet_cont_constraints(self,x_cont):
        return True

    #Constraints on the emission line fit parameters.
    def meet_line_constraints(self,x_line):
        return True

    #This function is called to determine that the fit can indeed be run.
    def can_fit_be_run(self,spec,verbose=True):
        
        #Check if centroid of the emission line is within the
        #spectrum.
        if spec.lam_rest is None or \
           self.line1_center<np.min(spec.lam_rest) or \
           self.line2_center>np.max(spec.lam_rest):
            if verbose:
                print("Lines not within spectrum")
            return False
        return True

    #This function is called to determine the indices of the spectrum
    #to be used for fitting the continuum.
    def get_i_cont(self,spec):
        i_cont = np.argwhere(
            ((spec.lam_rest>=self.continuum_regions[0][0]) &
             (spec.lam_rest<=self.continuum_regions[0][1])) |
            ((spec.lam_rest>=self.continuum_regions[1][0]) &
             (spec.lam_rest<=self.continuum_regions[1][1])))
    return i_cont

    #This is one is called to determine the indices of the spectrum
    #used for fitting the emission line.
    def get_i_line(self,spec):
        #Only consider regions within a certain velocity range of the
        #canonical emission line center.
        v1 = (c*(spec.lam_rest/self.line1_center-1.)).to(u.km/u.s)
        v1abs = np.abs(v1)
        v2 = (c*(spec.lam_rest/self.line2_center-1.)).to(u.km/u.s)
        v2abs = np.abs(v2)
        i_line = np.argwhere((v1abs<self.line_velocity_region) |
                             (v2abs<self.line_velocity_region))
    return i_line


