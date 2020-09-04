#Class to jointly fit N emission lines.

import numpy as np
import astropy.units as u
from astropy.constants import c
from astropy.stats import sigma_clipped_stats
import os
import re
import sys

from .line_class import Line_fit
from .MC_errors_general import get_error

class Multi_Line_fit(Line_fit):

    def __init__(self,_line_name,lines_file=None,lines_center_file=None,spec=None):

        #Load the main class.
        super(Multi_Line_fit,self).__init__(_line_name,spec)

        #Basic units to be used.
        self.flamunit = u.erg/(u.s*u.cm**2*u.AA)
        self.waveunit = u.AA
        self.vunit    = u.km/u.s
        self.flux_unit = u.erg/(u.s*u.cm**2)

        #Initial default guess value for line widths
        self.sigma_v_0 = 1000. * self.vunit

        #First, read the line centers.
        if lines_center_file is None:
            lines_center_file = os.environ['SPEC_PIPE_LOC'] + "/Spec_pipeline/Line_Fitter/line_centers.txt"
        line_centers = dict()
        cat = open(lines_center_file,"r")
        for line in cat:
            if line[0]=="#":
                continue
            x = line.split()
            line_centers[x[0]] = float(x[1])
        cat.close()

        #Search the list for the line in question.
        if lines_file is None:
            lines_file = os.environ['SPEC_PIPE_LOC']+\
                   "/Spec_pipeline/Line_Fitter/multi_lines.txt"
        cat = open(lines_file,"r")
        line_found = False
        for line in cat:
            if not line.strip():
                continue
            x = line.split()
            if x[0]==_line_name:
                line_found = True
                break
        cat.close()
        if not line_found:
            print("Error: Line",_line_name,"not found in file",lines_file)
            sys.exit()
        #x[1:] = [float(ix) for ix in x[1:]]

        #Define all the fit control variables that might not get
        #defined later.
        self.joint_dv = None
        self.joint_sigma = None
        self.fixed_ratio = None

        #Parse the data
        self.nlines = int(float(x[1]))

        #Assign parameters for the fit
        #Line centers - the lines file could provide the name of the line or its central wavelength.
        #self.line_center = np.array(x[2:2+self.nlines*1])*self.waveunit
        self.line_center = np.zeros(self.nlines)*self.waveunit
        for j,jj in enumerate(range(2,2+self.nlines*1)):
            try:
                self.line_center[j] = float(x[jj])*self.waveunit
            except:
                if x[jj] in line_centers.keys():
                    self.line_center[j] = line_centers[x[jj]]*self.waveunit
                else:
                    print("Warning: Line",x[jj],"not found in line centers file")
                    return

        #Joint constraints
        k = 2+self.nlines
        x[k:] = [float(ix) for ix in x[k:]]
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

        #Fit regions.
        self.line_velocity_region = x[k]*u.km/u.s
        self.ncont_reg = int((len(x)-(k+1))/2)
        self.continuum_regions = np.zeros((self.ncont_reg,2))
        for i in range(self.ncont_reg):
            self.continuum_regions[i][0] = x[i*2+k+1]
            self.continuum_regions[i][1] = x[i*2+k+2]
        self.continuum_regions = self.continuum_regions*u.AA

        #Set the default constraints for the fitting.
        self.dv_max      = np.ones(self.nlines)*500.*u.km/u.s
        self.sigma_v_max = np.ones(self.nlines)*5000.*u.km/u.s

        #If default spectrum has been defined, then we should use its resolution to set sigma_v_min.
        if spec is not None and spec.sigma_res is not None:
            self.sigma_v_min = (spec.sigma_res/(self.line_center*(1+spec.zspec))* c).to(self.vunit)
        else:
            self.sigma_v_min = np.ones(self.nlines)*100.*u.km/u.s

        return


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
                j = np.nonzero(self.joint_dv[:,i])[0][0]
                x_line_use[i*3] = x_line_use[j*3]

            if self.fixed_ratio[self.fixed_ratio[:,i]>0,i].size==0:
                x_line_use[i*3+1] = x_line_fit[k]
                k+=1
            else:
                #j = np.nonzero(self.fixed_ratio[:,i])[0][0]
                j = np.argwhere(self.fixed_ratio[:,i]>0)[0][0]
                x_line_use[i*3+1] = x_line_use[j*3+1]/self.fixed_ratio[j,i]

            if self.joint_sigma[self.joint_sigma[:,i]>0,i].size==0:
                x_line_use[i*3+2] = x_line_fit[k]
                k+=1
            else:
                j = np.nonzero(self.joint_sigma[:,i])[0][0]
                x_line_use[i*3+2] = x_line_use[j*3+2]


        return x_line_use


    ###################
    # Fit Model
    ###################

    #Complete model.
    def flam_model(self,lam,x=None,chain_output=None):

        if chain_output is not None:
            x_line = chain_output[:,:self.npar_line].T
            x_cont = chain_output[:,self.npar_line:].T
        elif x is not None:
            x_line = x[:self.npar_line]
            x_cont = x[self.npar_line:]
        else:
            x_line = None
            x_cont = None

        flam_model = self.flam_cont_model(lam,x_cont)
        flam_model += self.flam_line_model(lam,x_line)

        return flam_model


    #A Gaussian per line.
    def flam_line_model(self,lam,x_line=None):
        x_line_use = None
        if x_line is not None:
            x_line_use = self.line_par_translator(x_line)

        flam_line_model = None
        for i in range(self.nlines):
            dv, flam_line, sigma_v = self.line_par_parser(i,x_line_use)
            v = c*(lam/self.line_center[i]-1.)
            if flam_line_model is None:
                flam_line_model  = flam_line * np.exp(-0.5*((v-dv)/sigma_v)**2)
            else:
                flam_line_model += flam_line * np.exp(-0.5*((v-dv)/sigma_v)**2)

        return flam_line_model

    #Continuum will be a straight line.
    def flam_cont_model(self,lam,x_cont=None):
        a, b = self.cont_par_parser(x_cont)
        return a*lam+b


    def cont_par_parser(self,x_cont):
        if x_cont is None:
            a = self.a
            b = self.b
        else:
            a = x_cont[0]*self.flamunit/self.waveunit
            b = x_cont[1]*self.flamunit
        return a, b

    def line_par_parser(self,i,x_line_use):
        if x_line_use is None:
            dv        = self.dv_fit[i]
            flam_line = self.flam_line_fit[i]
            sigma_v   = self.sigma_v_fit[i]
        else:
            dv        = x_line_use[i*3]  *self.vunit
            flam_line = x_line_use[i*3+1]*self.flamunit
            sigma_v   = x_line_use[i*3+2]*self.vunit
        return dv, flam_line, sigma_v

    @property
    def line_SNR(self):
        #Can only run if an MC chain exists.
        if self.MC_chain is None:
            print("Need to run an MC first")
            return None

        #Get the line fluxes and their dispersion (used as noise).
        S = self.line_flux()
        #For the dispersion do 3 sigma clipping. Should not do less than about 500 realizations.
        #N = np.std(self.line_flux(MC=True),axis=1)
        N = sigma_clipped_stats(self.line_flux(MC=True),axis=1)[2]
        return (S/N).to(1.)

    #This implementation, instead of doing sigma clipping, estimates the noise using the steepness of the wings. Specifically, the noise is calculated as the 95.4% range - 68.3% range. For a Gaussian distribution, this should be equal to sigma.
    @property
    def line_SNR_wing(self):
        #Can only run if an MC chain exists.
        if self.MC_chain is None:
            print("Need to run an MC first")
            return None

        #Get the line fluxes and their dispersion (used as noise).
        S = self.line_flux()
        #For the dispersion do 3 sigma clipping. Should not do less than about 500 realizations.
        N = np.zeros(self.nlines)*S.unit
        for k in range(self.nlines):
            N_low, N_hig = get_error(self.line_flux(MC=True)[k],self.line_flux()[k])
            N2_low, N2_hig = get_error(self.line_flux(MC=True)[k],self.line_flux()[k],cf=95.4)
            N[k] = N2_low-N_low
        return (S/N).to(1.)

    #Line flux or fluxes, depending on the case.
    def line_flux(self,x_line=None,chain_output=None,MC=False):

        if MC:
            chain_output = self.MC_chain

        if chain_output is not None:
            x_line = chain_output[:,:self.npar_line].T

        if x_line is not None:
            x_line_use = self.line_par_translator(x_line)
        else:
            x_line_use = None

        if chain_output is not None:
            integrated_flux = np.zeros((self.nlines,len(chain_output)))*self.flux_unit
        else:
            integrated_flux = np.zeros(self.nlines)*self.flux_unit
        for i in range(self.nlines):
            integrated_flux[i] = self._get_line_flux(i,x_line_use)

        return integrated_flux

    def _get_line_flux(self,i,x_line_use=None):
        dv, flam_line, sigma_v = self.line_par_parser(i,x_line_use)
        flux = flam_line * self.line_center[i]*(2.*np.pi)**0.5 * \
               sigma_v/c
        flux = flux.to(self.flux_unit)
        return flux

    #Note that this expression for the EW is only valid for cases in which the continuum is either constant with wavelength or antisymmetric around the peak wavelength (such as the case for the linear continuum used here).
    def EW(self,x=None,chain_output=None,MC=False):

        if MC:
            chain_output = self.MC_chain

        if chain_output is not None:
            x_line = chain_output[:,:self.npar_line].T
            x_cont = chain_output[:,self.npar_line:].T
        elif x is not None:
            x_line = x[:self.npar_line]
            x_cont = x[self.npar_line:]
        else:
            x_line = None
            x_cont = None

        Fline = self.line_flux(x_line,chain_output=chain_output)
        if chain_output is not None:
            flam_cont = np.zeros((self.nlines,len(chain_output)))*self.flamunit
        else:
            flam_cont = np.zeros(self.nlines)*self.flamunit
        lam_peak = self.line_center*(1+self.dv_fit/c)
        for i in range(self.nlines):
            flam_cont[i] = self.flam_cont_model(lam_peak[i],x_cont)
        return (Fline/flam_cont).to(self.waveunit)

    #############
    # Constraints
    #############

    #Check if constraints are met.
    def meet_constraints(self,x):
        x_line = x[:self.npar_line]
        x_cont = x[self.npar_line:]
        return (self.meet_line_constraints(x_line) and self.meet_cont_constraints(x_cont))

    #Constraints on the continuum fit parameters.
    def meet_cont_constraints(self,x_cont):
        #Require continuum to be non-negative over the entire range. Since it is just a straight line, it is enough to check the edges of the regions.
        #a, b = self.cont_par_parser(x_cont)
        #edge_fluxes = a*self.continuum_regions.flatten()+b
        a = x_cont[0]
        b = x_cont[1]
        lam = self.continuum_regions.flatten().to(self.waveunit).value
        edge_fluxes = a*lam+b
        if len(edge_fluxes[edge_fluxes<0])>0:
            return False
        return True

    #Constraints on the emission line fit parameters.
    def meet_line_constraints(self,x_line):

        x_line_use = self.line_par_translator(x_line)

        for i in range(self.nlines):
            dv, flam_line, sigma_v = self.line_par_parser(i,x_line_use)

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

    def get_i_fit(self,spec):
        i_cont = self.get_i_cont(spec)
        i_line = self.get_i_line(spec)
        i_ban  = self.get_i_ban(spec)
        i_all  = np.unique(np.concatenate((i_cont,i_line)))
        if len(i_ban)>0:
            i_all = i_all[~np.in1d(i_all,i_ban)]
        return i_all

    #This function is called to determine the indices of the spectrum
    #to be used for fitting the continuum.
    def get_i_cont(self,spec):

        i_cont_use = []
        for k in range(self.ncont_reg):
            i_cont_use.append(
                np.argwhere((spec.lam_rest>=self.continuum_regions[k][0]) &
             (spec.lam_rest<=self.continuum_regions[k][1]))
            )
        i_cont_flat = [item for sublist in i_cont_use for item in sublist]
        i_cont = np.unique(i_cont_flat)
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

    #This function is called to get the indices of wavelengths ranges that should be excised from the any fit. The first obvious regions to be omitted are the atmospheric A-band and B-band absorption lines. Taken from Smette et al. (2015), section 2.
    #Updated on Sep 1, 2020 to numbers suggested by Dan Stern over email.
    def get_i_ban(self,spec):
        #A-band
        i_banA = np.argwhere((spec.lam_obs>=7592.*u.AA) & (spec.lam_obs<=7677.*u.AA))
        #B-band
        i_banB = np.argwhere((spec.lam_obs>=6864.*u.AA) & (spec.lam_obs<=6920.*u.AA))
        i_ban = np.concatenate([i_banA.flatten(),i_banB.flatten()])
        return i_ban

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
            if self.fixed_ratio[self.fixed_ratio[:,i]>0,i].size==0: \
               _npar_line+=1
        return _npar_line

    @property
    def npar_cont(self):
        return 2

    ######################
    # Parameter setters. #
    ######################

    def set_pars(self,x):
        x_line = x[:self.npar_line]
        x_cont = x[self.npar_line:]
        self.set_line_pars(x_line)
        self.set_cont_pars(x_cont)
        return


    def set_line_pars(self,x_line):

        self.dv_fit        = np.zeros(self.nlines)*self.vunit
        self.flam_line_fit = np.zeros(self.nlines)*self.flamunit
        self.sigma_v_fit   = np.zeros(self.nlines)*self.vunit

        x_line_use = self.line_par_translator(x_line)

        for i in range(self.nlines):
            self.dv_fit[i], self.flam_line_fit[i],\
                self.sigma_v_fit[i] = self.line_par_parser(i,x_line_use)

        return

    def set_cont_pars(self,x_cont):
        self.a, self.b = self.cont_par_parser(x_cont)


    def set_initial_fit_values(self, spec):

        #Check that the fit can be run.
        if not self.can_fit_be_run(spec,verbose=False):
            return

        self.x0_line = self.lines_initial_fit_values(spec)
        self.x0_cont = self.cont_initial_fit_values(spec)

        self.x0 = np.concatenate((self.x0_line,self.x0_cont))

        return

    def lines_initial_fit_values(self,spec):

        #Set up the initial values.
        dv_0      = np.zeros(self.nlines) #km/s
        sigma_v_0 = self.sigma_v_0.to(u.km/u.s).value*np.ones(self.nlines) #km/s

        flam_line_0 = np.zeros(self.nlines)
        for i in range(self.nlines):
            aux = np.max(spec.flam[np.abs(spec.lam_rest-self.line_center[i])<3*u.AA])*0.5
            flam_line_0[i] = aux.value

        x0_line = np.zeros(self.npar_line)
        x0_line[:3] = [dv_0[0], flam_line_0[0], sigma_v_0[0]]
        k = 3
        for i in range(1,self.nlines):
            if self.joint_dv[self.joint_dv[:,i]>0,i].size==0:
                x0_line[k] = dv_0[i]
                k+=1
            if self.fixed_ratio[self.fixed_ratio[:,i]>0,i].size==0:
                x0_line[k] = flam_line_0[i]
                k+=1
            if self.joint_sigma[self.joint_sigma[:,i]>0,i].size==0:
                x0_line[k] = sigma_v_0[i]
                k+=1
        return x0_line

    def cont_initial_fit_values(self,spec):
        i_cont = self.get_i_cont(spec)
        mean_cont = np.mean(spec.flam[i_cont]).to(self.flamunit)
        mean_lam  = np.mean(spec.lam_rest[i_cont]).to(self.waveunit)
        a0 = mean_cont.value/mean_lam.value
        b0 = 0.
        x0_cont = [a0,b0]
        return x0_cont

    def parse_chain_output(self,Output):

        x_line = Output[:,:self.npar_line].T
        x_cont = Output[:,self.npar_line:].T

        self.dv_low = np.zeros(self.nlines)*self.vunit
        self.dv_hig = np.zeros(self.nlines)*self.vunit
        self.flam_line_low = np.zeros(self.nlines)*self.flamunit
        self.flam_line_hig = np.zeros(self.nlines)*self.flamunit
        self.sigma_v_low = np.zeros(self.nlines)*self.vunit
        self.sigma_v_hig = np.zeros(self.nlines)*self.vunit

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

        self.line_flux_low = np.zeros(self.nlines)*self.flux_unit
        self.line_flux_hig = np.zeros(self.nlines)*self.flux_unit
        line_flux_MC = self.line_flux(MC=True)
        line_flux_fit = self.line_flux()
        self.EW_hig = np.zeros(self.nlines)*self.waveunit
        self.EW_low = np.zeros(self.nlines)*self.waveunit
        EW_MC = self.EW(MC=True)
        EW_fit = self.EW()
        for i in range(self.nlines):
            self.line_flux_low[i], self.line_flux_hig[i] = get_error(line_flux_MC[i],line_flux_fit[i])
            self.EW_low[i], self.EW_hig[i] = get_error(EW_MC[i],EW_fit[i])

        a, b = self.cont_par_parser(x_cont)
        self.a_low, self.a_hig = get_error(a, self.a)
        self.b_low, self.b_hig = get_error(b, self.b)

        return

    ##########
    # Plotting
    ##########

    #Textbox with fit results. We will use one box per emission line, below the spectrum. We have up to three emission lines.
    def line_legends(self, ax):

        import matplotlib.pyplot as plt

        line_name = re.sub("red","",self.line_name)
        lname = line_name.split("_")
        if len(lname)<self.nlines:
            lname = [lname[0]]*self.nlines
        props = dict(boxstyle='round', facecolor='white', alpha=0.75)
        for i in range(self.nlines):
            textbox = "{0:s}".format(lname[i])
            textbox += "\nFWHM = {0:.0f}".format(self.FWHM_v[i])
            textbox += "\n$\Delta v$ = {0:.1f}".format(self.dv_fit[i])
            if self.MC_chain is not None:
                textbox += "\nSNR = {0:.1f}".format(self.line_SNR[i])
            if self.p is not None:
                textbox += "\np   = {0:.3f}".format(self.p[i])
            plt.text(0.025+i/3., 0.025, textbox, transform=ax.transAxes, fontsize=10, verticalalignment='bottom', horizontalalignment='left', bbox=props)

        return

    ##########
    # Printing
    ##########

    @property
    def print_fit_header(self):
        print_output = ""
        for i in range(self.nlines):
            print_output += "{1:s}{0:d} {2:s}{0:d} {3:s}{0:d} ".format(
                i+1,"dv", "flam_line", "FWHM_v", "line_flux", "EW")
        return print_output

    @property
    def print_fit(self):
        print_output = ""
        for i in range(self.nlines):
            print_output += "{0:.3f} {1:.3e} {2:.3f} ".format(
                self.dv_fit[i].value,
                self.flam_line_fit[i].value,
                self.FWHM_v[i].value, self.line_flux()[i].value,
                self.EW()[i].value)
        return print_output

    @property
    def print_MC_header(self):
        print_output = ""
        for i in range(self.nlines):
            print_output += "{1:s}{0:d} ".format(i+1,"SNR")
            print_output += "{1:s}{0:d}_low {1:s}{0:d}_hig ".format(
                i+1,"dv")
            print_output += "{1:s}{0:d}_low {1:s}{0:d}_hig ".format(
                i+1,"flam_line")
            print_output += "{1:s}{0:d}_low {1:s}{0:d}_hig ".format(
                i+1,"FWHM_v")
            print_output += "{1:s}{0:d}_low {1:s}{0:d}_hig ".format(
                i+1,"line_flux")
            print_output += "{1:s}{0:d}_low {1:s}{0:d}_hig ".format(
                i+1,"EW")
        return print_output

    @property
    def print_MC(self):
        print_output = ""
        for i in range(self.nlines):
            print_output += "{0:.3e} ".format(self.line_SNR[i])
            print_output += "{0:.3f} {1:.3f} ".format(
                self.dv_low[i].value,
                self.dv_hig[i].value)
            print_output += "{0:.3e} {1:.3e} ".format(
                self.flam_line_low[i].value,
                self.flam_line_hig[i].value)
            print_output += "{0:.3f} {1:.3f} ".format(
                self.FWHM_v_low[i].value,
                self.FWHM_v_hig[i].value)
            print_output += "{0:.3f} {1:.3f} ".format(
                self.line_flux_low[i].value,
                self.line_flux_hig[i].value)
            print_output += "{0:.3f} {1:.3f} ".format(
                self.EW_low[i].value,
                self.EW_hig[i].value)
        return print_output

    def print_header(self,MC=False):
        print_output = "{0:15s} {1:7s} ".format("Obj_ID","z_spec")
        print_output += "{0:10s} ".format("line_id")
        print_output += "{0:10s} {1:10s} {2:10s} {3:10s} {4:10s} {5:10s} {6:10s} ".format("dv", "flam_line", "FWHM_v", "lam_cen", "line_flux", "EW","p")
        if MC:
            print_output += "{0:10s} ".format("SNR")
            print_output += "{0:10s} ".format("SNR_wing")
            print_output += "{0:14s} {1:14s} ".format("dv_low","dv_hig")
            print_output += "{0:14s} {1:14s} ".format( "flam_line_low", "flam_line_hig")
            print_output += "{0:14s} {1:14s} ".format( "FWHM_v_low", "FWHM_v_hig")
            print_output += "{0:14s} {1:14s} ".format( "lam_cen_low", "lam_cen_hig")
            print_output += "{0:14s} {1:14s} ".format( "line_flux_low", "line_flux_hig")
            print_output += "{0:14s} {1:14s} ".format( "EW_low", "EW_hig")
        print_output += "\n"
        return print_output

    def print(self,spec_use=None,MC=False):
        spec = self.get_spec_use(spec_use)
        self.run_Ftest()
        print_output = ""
        lname = self.line_name.split("_")
        if len(lname)<self.nlines:
            lname = [lname[0]]*self.nlines
        for i in range(self.nlines):
            print_output += "{0:15s} {1:7.3f} ".format(spec.name, spec.zspec)
            print_output += "{0:10s} ".format(lname[i])
            print_output += "{0:10.3f} {1:10.3e} {2:10.3f} {3:10.3f} {4:10.3e} {5:10.3f} {6:10.3f}  ".format( self.dv_fit[i].value, self.flam_line_fit[i].value, self.FWHM_v[i].value, self.line_center[i].value*(1+self.zline()[i]), self.line_flux()[i].value, self.EW()[i].value, self.p[i])
            if MC:
                print_output += "{0:10.3e} ".format(self.line_SNR[i])
                print_output += "{0:10.3e} ".format(self.line_SNR_wing[i])
                print_output += "{0:14.3f} {1:14.3f} ".format( self.dv_low[i].value, self.dv_hig[i].value)
                print_output += "{0:14.3e} {1:14.3e} ".format( self.flam_line_low[i].value,self.flam_line_hig[i].value)
                print_output += "{0:14.3f} {1:14.3f} ".format( self.FWHM_v_low[i].value, self.FWHM_v_hig[i].value)
                print_output += "{0:14.3f} {1:14.3f} ".format( self.line_center[i].value*(self.dv_low[i]/c).to(1.).value, self.line_center[i].value*(self.dv_hig[i]/c).to(1.).value)
                print_output += "{0:14.3e} {1:14.3e} ".format( self.line_flux_low[i].value, self.line_flux_hig[i].value)
                print_output += "{0:14.3f} {1:14.3f} ".format( self.EW_low[i].value, self.EW_hig[i].value)
            print_output += "\n"
        return print_output
