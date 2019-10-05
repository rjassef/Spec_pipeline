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

        #Define all the fit control variables that might not get
        #defined later.
        self.joint_dv = None
        self.joint_sigma = None
        self.fixed_ratio = None

        #Parse the data
        self.nlines = int(x[1])

        #Assign parameters for the fit
        self.line_center = np.array(x[2:2+self.nlines*1])*u.AA
        k = 2+self.nlines
        if self.nlines>1:
            self.joint_dv = np.zeros((self.nlines,self.nlines),dtype=np.int32)
            for i in range(self.nlines-1):
                j1 = np.sum(np.arange(self.nlines-i,self.nlines))+k
                j2 = np.sum(np.arange(self.nlines-(i+1),self.nlines))+k
                self.joint_dv[i][i+1:] = x[j1:j2]
            k+= np.sum(np.arange(self.nlines))
            self.joint_sigma = np.zeros((self.nlines,self.nlines),
                                        dtype=np.int32)
            for i in range(self.nlines-1):
                j1 = np.sum(np.arange(self.nlines-i,self.nlines))+k
                j2 = np.sum(np.arange(self.nlines-(i+1),self.nlines))+k
                self.joint_sigma[i][i+1:] = x[j1:j2]
            k+= np.sum(np.arange(self.nlines))
            self.fixed_ratio= np.zeros((self.nlines,self.nlines))
            for i in range(self.nlines-1):
                j1 = np.sum(np.arange(self.nlines-i,self.nlines))+k
                j2 = np.sum(np.arange(self.nlines-(i+1),self.nlines))+k
                self.fixed_ratio[i][i+1:] = x[j1:j2]
            k+= np.sum(np.arange(self.nlines))
            
        self.line_velocity_region = x[k]*u.km/u.s
        n_creg = int((len(x)-(k+1))/2)
        self.continuum_regions = np.zeros((n_creg,2))
        for i in range(n_creg):
            self.continuum_regions[i][0] = x[i*2+k+1]
            self.continuum_regions[i][1] = x[i*2+k+2]
        self.continuum_regions = self.continuum_regions*u.AA

        #Set the default constraints for the fitting.
        self.dv_max      = np.ones(self.nlines)*500.*u.km/u.s
        self.sigma_v_min = np.ones(self.nlines)*100.*u.km/u.s
        self.sigma_v_max = np.ones(self.nlines)*5000.*u.km/u.s

        
        #Finally, load the main class.
        super(Multi_Line_fit,self).__init__(_line_name)


    #Parameter translator.
    #x_line_fit is what is used in the fits,
    #x_line_use is what is used to construct the model.

    #Order is always [dv1, flam_line1, sigma_v1, dv2, flam_line2, ...]
    
    def line_par_translator(self,x_line_fit):

        if len(x_line_fit.shape)==1:
            x_line_use = np.zeros(self.nlines*3)
        else:
            x_line_use = np.zeros((self.nlines*3,x_line_fit.shape[1]))
                                  
        #All parameters for the first line are always fit.
        x_line_use[:3] = x_line_fit[:3]

        k = 3
        for i in range(1,self.nlines):

            #For line i, see if any of the dv values are greater than
            #0. Stop at the first one different than 0 if any. If none
            #are greater than 0, then the dv parameter is meant to be
            #fit. Repeat then for the line ratio and sigma.
            if self.joint_dv[self.joint_dv[:,i]>0,i].size==0:
                x_line_use[i*3] = x_line_fit[k]
                k+=1
            else:
                j = np.nonzero(self.joint_dv[:,i])[0]
                x_line_use[i*3] = x_line_use[j*3]
                
            if self.joint_sigma[self.joint_sigma[:,i]>0,i].size==0:
                x_line_use[i*3+1] = x_line_fit[k]
                k+=1
            else:
                j = np.nonzero(self.joint_sigma[:,i])[0]
                x_line_use[i*3+1] = x_line_use[j*3]

            if self.fixed_ratio[self.fixed_ratio[:,i]>=0,i].size==0:
                x_line_use[i*3+2] = x_line_fit[k]
                k+=1
            else:
                j = np.nonzero(self.fixed_ratio[:,i])[0]
                x_line_use[i*3+2] = x_line_use[j*3]*self.fixed_ratio[j,i]

        return x_line_use
        
        
    ###################
    # Fit Model
    ###################
        
    
    #Complete model.
    def flam_model(self,lam,x_line=None,x_cont=None,chain_output=None):

        if chain_output is not None:
            x_line = chain_output[:,:self.npar_line].T
            x_cont = chain_output[:,self.npar_line:].T

        if x_line is not None:
            x_line_use = self.line_par_translator(x_line)
        else:
            x_line_use = None
            
        flam_model = self.flam_cont_model(lam,x_cont)
        for i in range(self.nlines):
            flam_model += self.flam_line_model(lam,i,x_line_use)
            
        return flam_model


    #A Gaussian per line.
    def flam_line_model(self,lam,i,x_line_use=None):
        dv, flam_line, sigma_v = self.line_par_parser(i,x_line_use)
        v = c*(lam/self.line_center[i]-1.)
        return flam_line * np.exp(-0.5*((v-dv)/sigma_v)**2)

    #Continuum will be a straight line.
    def flam_cont_model(self,lam,x_cont=None):
        a, b = self.cont_par_parser(x_cont)
        return a*lam+b

    
    def cont_par_parser(self,x_cont):
        if x_cont is None:
            a = self.a
            b = self.b
        else:
            a = x_cont[0]*u.erg/u.cm**2/u.s/u.AA**2
            b = x_cont[1]*u.erg/u.cm**2/u.s/u.AA
        return a, b

    def line_par_parser(self,i,x_line_use):
        if x_line_use is None:
            dv        = self.dv_fit[i]
            flam_line = self.flam_line_fit[i]
            sigma_v   = self.sigma_v_fit[i]
        else:
            dv        = x_line_use[i*3]*u.km/u.s
            flam_line = x_line_use[i*3+1]*u.erg/u.cm**2/u.s/u.AA
            sigma_v   = x_line_use[i*3+2]*u.km/u.s 
        return dv, flam_line, sigma_v


    #############
    # Constraints
    #############

    #Constraints on the continuum fit parameters.
    def meet_cont_constraints(self,x_cont):
        return True

    #Constraints on the emission line fit parameters.
    def meet_line_constraints(self,x_line):

        x_line_use = self.line_par_translator(x_line)

        for i in range(self.nlines):
            dv, flam_line, sigma_v = self.line_par_parser(i,x_line)

            #Do no allow huge shifts on the line centers
            if np.abs(dv)>self.dv_max[i]:
                return False

            #Do not allow too narrow or too broad emission lines.
            if sigma_v < self.sigma_v_min[i] or \
               sigma_v > self.sigma_v_max[i]: 
                return False

            #Do no allow negative emission lines.
            if flam_line<0.:
                return False

        return True

    #This function is called to determine that the fit can indeed be run.
    def can_fit_be_run(self,spec,verbose=True):
        
        #Check if centroid of the emission line is within the
        #spectrum.
        min_lam_targ = np.min(self.line_center)
        max_lam_targ = np.max(self.line_center)
        if spec.lam_rest is None or \
           min_lam_targ<np.min(spec.lam_rest) or \
           max_lam_targ>np.max(spec.lam_rest):
            if verbose:
                print("Lines not within spectrum")
            return False
        return True

    ##################
    # Index functions
    ##################
    
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
        lam_min = np.min(self.line_center)
        lam_max = np.max(self.line_center)
        
        v1 = (c*(spec.lam_rest/lam_min-1.)).to(u.km/u.s)
        v2 = (c*(spec.lam_rest/lam_max-1.)).to(u.km/u.s)
        v1abs = np.abs(v1)
        v2abs = np.abs(v2)
        
        i_line1 = np.argwhere(v1abs<self.line_velocity_region)
        i_line2 = np.argwhere(v2abs<self.line_velocity_region)

        k_min = np.min([np.min(i_line1),np.min(i_line2)])
        k_max = np.max([np.max(i_line1),np.max(i_line2)])

        i_line = np.arange(k_min,k_max+1,dtype=np.int32)

        return i_line


    #####################
    # Useful properties.
    #####################


    def sigma_to_FWHM(self,sigma_v):
        if sigma_v is None:
            return None
        return sigma_v*2.*(2.*np.log(2.))**0.5

    @property
    def FWHM_v(self):
        return self.sigma_to_FWHM(self.sigma_v_fit)
    @property
    def FWHM_v_low(self):
        return self.sigma_to_FWHM(self.sigma_v_low)
    @property
    def FWHM_v_hig(self):
        return self.sigma_to_FWHM(self.sigma_v_hig)

    @property
    def npar_fit(self):
        _npar_fit = self.npar_line + self.npar_cont
        return _npar_fit

    @property
    def npar_line(self):
        _npar_line = 3
        for i in range(1,self.nlines):
            if self.joint_dv[self.joint_dv[:,i]>0,i].size==0: _npar_line+=1
            if self.joint_sigma[self.joint_sigma[:,i]>0,i].size==0: \
               _npar_line+=1
            if self.fixed_ratio[self.fixed_ratio[:,i]>=0,i].size==0: \
               _npar_line+=1
        return _npar_line

    @property
    def npar_cont(self):
        return 2

    ######################
    # Parameter setters. #
    ######################

    def set_line_pars(self,x_line):

        self.dv_fit        = np.zeros(self.nlines)*u.km/u.s
        self.flam_line_fit = np.zeros(self.nlines)*u.erg/u.s/u.cm**2/u.AA
        self.sigma_v_fit   = np.zeros(self.nlines)*u.km/u.s
        
        x_line_use = self.line_par_translator(x_line)
        
        for i in range(self.nlines):
            self.dv_fit[i], self.flam_line_fit[i],\
                self.sigma_v_fit[i] = self.line_par_parser(i,x_line)

        return

    def set_cont_pars(self,x_cont):
        self.a, self.b = self.cont_par_parser(x_cont)


    def set_initial_fit_values(self, spec):

        #Check that the fit can be run.
        if not self.can_fit_be_run(spec,verbose=False):
            return

        #Set up the initial values.
        dv_0      = np.zeros(self.nlines) #km/s
        sigma_v_0 = 1000.*np.ones(self.nlines) #km/s
        
        flam_line_0 = np.zeros(self.nlines)
        for i in range(self.nlines):
            aux = np.max(
                spec.flam[np.abs(spec.lam_rest-self.line_center[i])<3*u.AA]
            )*0.5
            flam_line_0[i] = aux.value
        
        self.x0_line = np.zeros(self.npar_line)
        self.x0_line[:3] = [dv_0[0], flam_line_0[0], sigma_v_0[0]]
        k = 3
        for i in range(1,self.nlines):
            if self.joint_dv[self.joint_dv[:,i]>0,i].size==0:
                self.x0_line[k] = dv_0[i]
                k+=1
            if self.joint_sigma[self.joint_sigma[:,i]>0,i].size==0:
                self.x0_line[k] = flam_line_0[i]
                k+=1
            if self.fixed_ratio[self.fixed_ratio[:,i]>=0,i].size==0:
                self.x0_line[k] = sigma_v_0[i]
                k+=1

        self.x0_cont = [1.,0.]
        return

    def parse_chain_output(self,Output):
        
        x_line = Output[:,:self.npar_line].T
        x_cont = Output[:,self.npar_line:].T

        self.dv_low = np.zeros(self.nlines)*u.km/u.s
        self.dv_hig = np.zeros(self.nlines)*u.km/u.s
        self.flam_line_low = np.zeros(self.nlines)*u.erg/u.s/u.cm**2/u.AA
        self.flam_line_hig = np.zeros(self.nlines)*u.erg/u.s/u.cm**2/u.AA
        self.sigma_v_low = np.zeros(self.nlines)*u.km/u.s
        self.sigma_v_hig = np.zeros(self.nlines)*u.km/u.s
        
        x_line_use = self.line_par_translator(x_line)
        for i in range(self.nlines):
            dv, flam_line, sigma_v = self.line_par_parser(i,x_line_use)
            self.dv_low[i], self.dv_hig[i] = get_error(dv, self.dv_fit[i])
            self.flam_line_low[i], \
                self.flam_line_hig[i] = get_error(flam_line,
                                                  self.flam_line_fit[i])
            self.sigma_v_low[i], \
                self.sigma_v_hig[i] = get_error(sigma_v,
                                                self.sigma_v_fit[i])

        a, b = self.cont_par_parser(x_cont)
        self.a_low, self.a_hig = get_error(a, self.a)
        self.b_low, self.b_hig = get_error(b, self.b)

        return
     
