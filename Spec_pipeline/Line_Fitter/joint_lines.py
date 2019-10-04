#This is a class that will fit two lines simultaneously, allowing for
#different types of constraints. 


import numpy as np
import astropy.units as u
from astropy.constants import c
import os

from .line_class import Line_fit
from .MC_errors_general import get_error

class Joint_Line_fit(Line_fit):

    def __init__(self,_line_name):

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
        self.dv1_low        = None
        self.dv1_hig        = None
        self.dv2_low        = None
        self.dv2_hig        = None
        self.flam_line1_low = None
        self.flam_line1_hig = None
        self.flam_line2_low = None
        self.flam_line2_hig = None
        self.sigma_v1_low   = None
        self.sigma_v1_hig   = None
        self.sigma_v2_low   = None
        self.sigma_v2_hig   = None
        
        #Set the default constraints for the fitting.
        self.dv1_max      = 500.*u.km/u.s
        self.dv2_max      = 300.*u.km/u.s
        self.sigma_v1_min =  100.*u.km/u.s
        self.sigma_v1_max = 5000.*u.km/u.s
        self.sigma_v2_min =  100.*u.km/u.s
        self.sigma_v2_max = 5000.*u.km/u.s

        #Search the list for the line in question.
        cat = open(os.environ['SPEC_PIPE_LOC']+\
                   "/Spec_pipeline/Line_Fitter/joint_lines.txt","r")
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
        self.line_velocity_region = x[6]*u.km/u.s
        self.continuum_regions = np.zeros((int((len(x)-7)/2),2))
        for i in range((int((len(x)-7)/2))):
            self.continuum_regions[i][0] = x[i*2+7]
            self.continuum_regions[i][1] = x[i*2+8]
        self.continuum_regions = self.continuum_regions*u.AA

        #Initialize the class
        super(Joint_Line_fit,self).__init__(_line_name)

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

        if chain_output is not None:
            x_line = chain_output[:,:self.npar_line].T
            x_cont = chain_output[:,self.npar_line:].T
            
        return self.flam_cont_model(lam,x_cont)+\
            self.flam_line_model(lam,x_line)

    #A Gaussian per line.
    def flam_line_model(self,lam,x_line):

        dv1, flam_line1, sigma_v1,\
            dv2, flam_line2, sigma_v2 = self.line_par_parser(x_line)
            
        v1 = c*(lam/self.line1_center-1.)
        v2 = c*(lam/self.line2_center-1.)
        
        line1 = flam_line1 * np.exp(-0.5*((v1-dv1)/sigma_v1)**2)
        line2 = flam_line2 * np.exp(-0.5*((v2-dv2)/sigma_v2)**2)
        return line1 + line2

    #Continuum will be a straight line.
    def flam_cont_model(self,lam,x_cont=None):
        a, b = self.cont_par_parser(x_cont)
        return a*lam+b

    #
    def cont_par_parser(self,x_cont):
        if x_cont is None:
            a = self.a
            b = self.b
        else:
            a = x_cont[0]*u.erg/u.cm**2/u.s/u.AA**2
            b = x_cont[1]*u.erg/u.cm**2/u.s/u.AA
        return a, b
    
    def line_par_parser(self,x_line):
        if x_line is None:
            dv1 = self.dv1_fit
            dv2 = self.dv2_fit
            flam_line1 = self.flam_line1_fit
            flam_line2 = self.flam_line2_fit
            sigma_v1 = self.sigma_v1_fit
            sigma_v2 = self.sigma_v2_fit
        else:
            dv1        = x_line[0]*u.km/u.s
            flam_line1 = x_line[1]*u.erg/u.s/u.cm**2/u.AA
            sigma_v1   = x_line[2]*u.km/u.s
            
            k=3
            if self.joint_dv==0:
                dv2 = x_line[k]*u.km/u.s
                k+=1
            else:
                dv2 = dv1
            if self.fixed_ratio<0:
                flam_line2 = x_line[k]*u.erg/u.s/u.cm**2/u.AA
                k+=1
            else:
                flam_line2 = flam_line1 * self.fixed_ratio
            if self.joint_sigma==0:
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
        dv1, flam_line1, sigma_v1,\
            dv2, flam_line2, sigma_v2 = self.line_par_parser(x_line)

        #Do no allow huge shifts on the line centers
        if np.abs(dv1)>self.dv1_max or np.abs(dv2)>self.dv2_max:
            return False

        #Do not allow too narrow or too broad emission lines.
        if sigma_v1 < self.sigma_v1_min or  sigma_v1 > self.sigma_v1_max: 
            return False
        if sigma_v2 < self.sigma_v2_min or  sigma_v2 > self.sigma_v2_max: 
            return False

        #Do no allow negative emission lines.
        if flam_line1<0. or flam_line2<0.:
            return False

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


    #####################
    # Useful properties.
    #####################


    def sigma_to_FWHM(self,sigma_v):
        if sigma_v is None:
            return None
        return sigma_v*2.*(2.*np.log(2.))**0.5

    @property
    def FWHM_v1(self):
        return self.sigma_to_FWHM(self.sigma_v1_fit)
    @property
    def FWHM_v1_low(self):
        return self.sigma_to_FWHM(self.sigma_v1_low)
    @property
    def FWHM_v1_hig(self):
        return self.sigma_to_FWHM(self.sigma_v1_hig)

    @property
    def FWHM_v2(self):
        return self.sigma_to_FWHM(self.sigma_v2_fit)
    @property
    def FWHM_v2_low(self):
        return self.sigma_to_FWHM(self.sigma_v2_low)
    @property
    def FWHM_v2_hig(self):
        return self.sigma_to_FWHM(self.sigma_v2_hig)

    @property
    def npar_fit(self):
        _npar_fit = self.npar_line + self.npar_cont
        return _npar_fit

    @property
    def npar_line(self):
        _npar_line = 3
        if self.joint_dv==0    : _npar_line += 1
        if self.fixed_ratio<0  : _npar_line += 1
        if self.joint_sigma==0 : _npar_line += 1
        return _npar_line

    @property
    def npar_cont(self):
        return 2

    ######################
    # Parameter setters. #
    ######################

    def set_line_pars(self,x_line):
        self.dv1_fit, self.flam_line1_fit, self.sigma_v1_fit,\
            self.dv2_fit, self.flam_line2_fit, \
            self.sigma_v2_fit = self.line_par_parser(x_line)
        return

    def set_cont_pars(self,x_cont):
        self.a, self.b = self.cont_par_parser(x_cont)


    def set_initial_fit_values(self, spec):

        #Check that the fit can be run.
        if not self.can_fit_be_run(spec,verbose=False):
            return

        #Set up the initial values.
        dv1_0 = 0. #km/s
        dv2_0 = 0. #km/s
        sigma_v1_0 = 1000. #km/s
        sigma_v2_0 = 1000. #km/s

        flam_line1_0 = \
            (np.max(
                spec.flam[np.abs(spec.lam_rest-self.line1_center)<3*u.AA]
            ))*0.5
        flam_line2_0 = \
            (np.max(
                spec.flam[np.abs(spec.lam_rest-self.line1_center)<3*u.AA]
            ))*0.5
        

        self.x0_line = np.zeros(self.npar_line)
        self.x0_line[:3] = [dv1_0, flam_line1_0.value,sigma_v1_0]
        k = 3
        if self.joint_dv==0:
            self.x0_line[k] = dv2_0
            k+=1
        if self.fixed_ratio<0:
            self.x0_line[k] = flam_line2_0.value
            k+=1
        if self.joint_sigma==0:
            self.x0_line[k] = sigma_v2_0
            k+=1

        self.x0_cont = [1.,0.]
        return

    def parse_chain_output(self,Output):
        
        x_line = Output[:,:self.npar_line].T
        x_cont = Output[:,self.npar_line:].T

        dv1, flam_line1, sigma_v1,\
            dv2, flam_line2, \
            sigma_v2 = self.line_par_parser(x_line)
        a, b = self.cont_par_parser(x_cont)

        self.dv1_low, self.dv1_hig = get_error(dv1, self.dv1_fit)
        self.dv2_low, self.dv2_hig = get_error(dv2, self.dv2_fit)
        self.flam_line1_low, self.flam_line1_hig = \
            get_error(flam_line1, self.flam_line1_fit)
        self.flam_line2_low, self.flam_line2_hig = \
            get_error(flam_line2, self.flam_line2_fit)
        self.sigma_v1_low, self.sigma_v1_hig = \
            get_error(sigma_v1, self.sigma_v1_fit)
        self.sigma_v2_low, self.sigma_v2_hig = \
            get_error(sigma_v2, self.sigma_v2_fit)

        self.a_low, self.a_hig = get_error(a, self.a)
        self.b_low, self.b_hig = get_error(b, self.b)

        return
        
