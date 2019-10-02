import numpy as np
import astropy.units as u
from astropy.constants import c
import os

from .line_class import Line_fit
from .MC_errors_general import get_error

class Default_Line_fit(Line_fit):
    
    def __init__(self,_line_name):

        #Values that will come from the fitting.
        self.lam_cen_fit       = None
        self.flam_line_cen_fit = None
        self.sigma_v_fit       = None
        self.a                 = None
        self.b                 = None

        #Values that will come from the MC
        self.lam_cen_low = None
        self.lam_cen_hig = None
        self.flam_line_cen_low = None
        self.flam_line_cen_hig = None
        self.sigma_v_low = None
        self.sigma_v_hig = None
        self.a_low = None
        self.a_hig = None
        self.b_low = None
        self.b_hig = None

        #Set the default constraints for the fitting.
        self.sigma_v_min =  100.*u.km/u.s
        self.sigma_v_max = 5000.*u.km/u.s
        self.delta_lam_cen_max = 2.*u.AA
        
        #Minimum peak height fraction compared to initial guess.
        self.frac_flam_line_cen_min = 1e-1

        
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

    ###################
    # Fit Model
    ###################

    #Complete model.
    def flam_model(self,lam,x_line=None,x_cont=None,chain_output=None):
        
        if chain_output is not None:
            x_opt = chain_output[:,:3].T
            x_cont = chain_output[:,3:].T
            
        return self.flam_cont_model(lam,x_cont)+\
            self.flam_line_model(lam,x_line)


    #Our default line fitting class will have a Gaussian Profile.
    def flam_line_model(self, lam, x_line=None):
        
        if x_line is not None:
            lam_cen       = x_line[0]*u.AA
            flam_line_cen = x_line[1]*u.erg/u.cm**2/u.s/u.AA
            sigma_v       = x_line[2]*u.km/u.s
        else:
            lam_cen       = self.lam_cen_fit
            flam_line_cen = self.flam_line_cen_fit
            sigma_v       = self.sigma_v_fit
            
        v = c*(lam/lam_cen-1.)
        return flam_line_cen * np.exp(-0.5*(v/sigma_v)**2)

    #Continuum will be a straight line.
    def flam_cont_model(self,lam,x_cont=None):
        if x_cont is not None:
            a = x_cont[0]*u.erg/u.cm**2/u.s/u.AA**2
            b = x_cont[1]*u.erg/u.cm**2/u.s/u.AA
        else:
            a = self.a
            b = self.b
        return a*lam+b

    ###############
    # Constraints
    ###############

    #Constraints on the continuum fit parameters.
    def meet_cont_constraints(self,x_cont):
        return True

    #Constraints on the emission line fit parameters.
    def meet_line_constraints(self,x_line):
            
        lam_cen       = x_line[0] * u.AA
        flam_line_cen = x_line[1] * u.erg/u.s/u.cm**2/u.AA
        sigma_v       = x_line[2] * u.km/u.s
   
        lam_cen_0       = self.x0_line[0] * u.AA
        flam_line_cen_0 = self.x0_line[1] * u.erg/u.s/u.cm**2/u.AA

        #Velocity width is within bounds.
        if sigma_v<self.sigma_v_min or sigma_v>self.sigma_v_max:
            return False

        #Centroid displacement is within bounds.
        if np.abs(lam_cen-lam_cen_0)>self.delta_lam_cen_max:
            return False

        #Do not allow absorption lines.
        if flam_line_cen.value<0.:
            return False

        #Do not allow the peak of the emission line to drift much below
        #the guess value. There is a failure mode on which the lines
        #are fit at any central wavelength with any width but with a flux
        #of effectively 0.
        if flam_line_cen_0>0 and \
           (flam_line_cen/flam_line_cen_0).to(1.)<\
           self.frac_flam_line_cen_min:
            return False

        return True


    #This function is called to determine that the fit can indeed be run.
    def can_fit_be_run(self,spec,verbose=True):
        
        #Check if centroid of the emission line is within the
        #spectrum.
        if spec.lam_rest is None or \
           self.line_center<np.min(spec.lam_rest) or \
           self.line_center>np.max(spec.lam_rest):
            if verbose:
                print("Line not within spectrum")
            return False
        return True
        

    #####################
    # Useful properties.
    #####################

    @property
    def FWHM_v_low(self):
        if self.sigma_v_low is None:
            print("First run the MC")
            return
        return self.sigma_v_low*2.*(2.*np.log(2.))**0.5

    @property
    def FWHM_v_hig(self):
        if self.sigma_v_hig is None:
            print("First run the MC")
            return
        return self.sigma_v_hig*2.*(2.*np.log(2.))**0.5

    @property
    def FWHM_v(self):
        if self.sigma_v_fit is None:
            print("First fit the line")
            return
        return self.sigma_v_fit*2.*(2.*np.log(2.))**0.5

        
     
    ######################
    # Parameter setters. #
    ######################
    
    def set_cont_pars(self,x_cont):
        self.a = x_cont[0]*u.erg/u.cm**2/u.s/u.AA**2
        self.b = x_cont[1]*u.erg/u.cm**2/u.s/u.AA
        return

    def set_line_pars(self,x_line):
        self.lam_cen_fit       = x_line[0]*u.AA
        self.flam_line_cen_fit = x_line[1]*u.erg/u.cm**2/u.s/u.AA
        self.sigma_v_fit       = x_line[2]*u.km/u.s
        return

    def set_initial_fit_values(self, spec, lam_cen_0=None,
                              sigma_v_0=3000.*u.km/u.s):

        #Check that the fit can be run.
        if not self.can_fit_be_run(spec,verbose=False):
            return

        #Set up the initial values.
        if lam_cen_0 is None:
            lam_cen_0 = self.line_center
        flam_line_cen_0 = \
            (np.max(spec.flam[np.abs(spec.lam_rest-lam_cen_0)<3*u.AA]))*0.5

        
        self.x0_line = [lam_cen_0.to(u.AA).value,
                        flam_line_cen_0.to(u.erg/u.cm**2/u.s/u.AA).value,
                        sigma_v_0.to(u.km/u.s).value]
        self.x0_cont = [1.,0.]
        return

    
    def parse_chain_output(self,Output):

        lam_cen       = Output[:,0] * u.AA
        flam_line_cen = Output[:,1] * u.erg/u.s/u.cm**2/u.AA
        sigma_v       = Output[:,2] * u.km/u.s
        a             = Output[:,3] * u.erg/u.s/u.cm**2/u.AA**2
        b             = Output[:,4] * u.erg/u.s/u.cm**2/u.AA
        
        self.lam_cen_low, self.lam_cen_hig = \
            get_error(lam_cen, self.lam_cen_fit)
        
        self.flam_line_cen_low, \
            self.flam_line_cen_hig = get_error(flam_line_cen, 
                                               self.flam_line_cen_fit)
        self.sigma_v_low, self.sigma_v_hig = get_error(sigma_v, 
                                                       self.sigma_v_fit)
        self.a_low, self.a_hig = get_error(a, self.a)
        self.b_low, self.b_hig = get_error(b, self.b)

        return
