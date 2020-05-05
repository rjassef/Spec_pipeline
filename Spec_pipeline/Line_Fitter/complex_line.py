import numpy  as np
import astropy.units as u
from astropy.constants import c

from .multi_line import Multi_Line_fit
from .line_class import Line_fit
from .MC_errors_general import get_error


#This will be its own subclass of Line_fit, actually, but will be made from combinations of Multi_Line_fit objects.

#We will defferiate between broad and narrow emission lines. In this model, all narrow emission lines have the same width and systemic velocity shift, while this can be different for each of the broad emission lines.

class Complex_Line(Line_fit):

    def __init__(self, line_name, spec=None):
        super(Complex_Line).__init__(line_name, spec)
        self.multi_line = []
        self.ntot_lines = 0
        self.n_multilines = 0

        self.x0 = None
        self.xopt = None

        self.flamunit = u.erg/(u.s*u.cm**2*u.AA)
        self.waveunit = u.AA
        self.vunit    = u.km/u.s

        self.dv_fit = None
        self.flam_line_fit = None
        self.sigma_v_fit = None
        self.a = None
        self.b = None

        return

    ###

    #The continuum regions will be the combination of all continuum regions.
    @property
    def continuum_regions(self):
        continuum_regions = self.multi_line[0].continuum_regions
        for i in range(1,self.n_multilines):
            continuum_regions = np.concatenate((continuum_regions,self.multi_line[i].continuum_regions))
        return continuum_regions

    @propety
    def ncont_reg(self):
        return len(continuum_regions)

    def add_line(self,line_name,width_type=None):

        #Setup the emission line.
        self.multi_line.append(Multi_Line_fit(line_name))
        self.multi_line[-1].width_type = width_type

        #If broad, FWHM>=1000 km/s. If narrow, FWHM<=1000 km/s
        if width_type == 'broad':
            self.multi_line[-1].sigma_v_min = 1000.*self.vunit /(2.*(2.*np.log(2.))**0.5)
        elif width_type == 'narrow':
            self.multi_line[-1].sigma_v_max = 1000.*self.vunit /(2.*(2.*np.log(2.))**0.5)
        elif width_type is not None:
            print("Unknown emission line type: ",width_type)
            print("Assuming full bounds for line ",line_name)

        self.ntot_lines += multi_line[-1].nlines
        self.n_multilines += 1

        return

    def set_pars(self,x):
        x_line = x[:self.npar_line]
        x_cont = x[self.npar_line:]
        self.set_line_pars(x_line)
        self.set_cont_pars(x_cont)
        return

    def get_i_fit(self,spec):
        i_cont = self.get_i_cont(spec)
        i_line = self.get_i_line(spec)
        i_all  = np.unique(np.concatenate((i_cont,i_line)))
        return i_all

    def meet_constraints(self,x):
        x_line = x[:self.npar_line]
        x_cont = x[self.npar_line:]
        return (self.meet_line_constraints(x_line) and self.meet_cont_constraints(x_cont))

    ###

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

        flam_model  = self.flam_cont_model(lam,x_cont)
        flam_model += self.flam_line_model(lam,x_line)

        return flam_model

    def set_initial_fit_values(self, spec):
        x0_line = self.lines_initial_fit_values(spec)
        x0_cont = self.cont_initial_fit_values(spec)
        self.x0 = np.concatenate((x0_line,x0_cont))
        return

    @property
    def npar_fit(self):
        return self.npar_line + self.npar_cont

    #####
    # Continuum Model

    #For now, we define a linear continuum.
    def flam_cont_model(lam,x_cont=None):
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

    def cont_initial_fit_values(self,spec):
        i_cont = self.get_i_cont(spec)
        mean_cont = np.mean(spec.flam[i_cont]).to(self.flamunit)
        mean_lam  = np.mean(spec.lam_rest[i_cont]).to(self.waveunit)
        a0 = mean_cont.value/mean_lam.value
        b0 = 0.
        self.x0_cont = [a0,b0]

    @property
    def npar_cont(self):
        return 2

    def set_cont_pars(self,x_cont):
        self.a, self.b = self.cont_par_parser(x_cont)

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

    def meet_cont_constraints(self,x_cont):
        return True

    ####
    # Lines Model

    def flam_line_model(self,lam,x_line=None):

        #Note that n_multilines counts the number of Multi_Line_fit objects, but each of those can be composed of more than 1 line.
        j1=0
        flam_line_model = np.zeros(self.ntot_lines)*self.flamunit
        dv_narrow = None #We'll set it to the width of the first narrow line.
        j1 = 0
        for k,line in enumerate(self.multi_line):
            j0=j1
            j1+=line.npar_line
            x_line_use = line.line_par_translator(x_line[j0:j1])
            for i in range(line.nlines):
                dv, flam_line, sigma_v = line.line_par_parser(i,x_line_use)
                if line.width_type == "narrow":
                    if dv_narrow is not None:
                        dv = dv_narrow
                    else:
                        dv_narrow = dv
                v = c*(lam/line.line_center[i]-1.)
                flam_line_model[k+i] = flam_line * np.exp(-0.5*((v-dv)/sigma_v)**2)
        return np.sum(flam_line_model,axis=0)

    def lines_initial_fit_values(self,spec):
        #For each line declared, get the starting values.
        x0 = self.multi_line[0].lines_initial_fit_values(spec)
        for i in range(1,self.n_multilines):
            x0 = np.concatenate((x0,self.multi_line[i].lines_initial_fit_values(spec)))
        return x0

    @property
    def npar_line(self):
        npar_line = 0
        for line in self.multi_line:
            npar_line += line.npar_line
        return npar_line

    def set_line_pars(self,x_line):

        self.dv_fit        = np.zeros(self.ntot_lines)*self.vunit
        self.flam_line_fit = np.zeros(self.ntot_lines)*self.flamunit
        self.sigma_v_fit   = np.zeros(self.ntot_lines)*self.vunit

        j1 = 0
        for k,line in enumerate(self.multi_line):
            j0=j1
            j1+=line.npar_line
            x_line_use = line.line_par_translator(x_line[j0:j1])
            for i in range(line.nlines):
                dv, flam_line, sigma_v = line.line_par_parser(i,x_line_use)
                self.dv_fit[k+i] = dv
                self.flam_line_fit[k+i] = flam_line
                self.sigma_v_fit[k+i] = sigma_v

        return

    def get_i_line(self,spec):
        i_all = self.multi_line[0].get_i_line(spec)
        for i in range(1,self.n_multilines):
            i_all = np.concat((i_all,self.multi_line[i].get_i_line(spec)))
        return np.unique(i_all)

    def meet_line_constraints(self,x_line):
        j1=0
        dv_narrow = None
        for k,line in enumerate(self.multi_line):
            j0=j1
            j1+=line.npar_line
            x_line_use = line.line_par_translator(x_line[j0:j1])
            for i in range(line.nlines):
                dv, flam_line, sigma_v = line.line_par_parser(i,x_line_use)
                if line.width_type == "narrow":
                    if dv_narrow is not None:
                        dv = dv_narrow
                    else:
                        dv_narrow = dv

                #Do no allow huge shifts on the line centers
                if np.abs(dv)>line.dv_max[i]:
                    return False

                #Do not allow too narrow or too broad emission lines.
                if sigma_v<line.sigma_v_min[i] or sigma_v>line.sigma_v_max[i]:
                    return False

                #Do no allow negative emission lines.
                if flam_line<0.:
                    return False

        return True

    ####

    def parse_chain_output(self,Output):

        x_line = Output[:,:self.npar_line].T
        x_cont = Output[:,self.npar_line:].T

        self.dv_low        = np.zeros(self.ntot_lines)*self.vunit
        self.dv_hig        = np.zeros(self.ntot_lines)*self.vunit
        self.flam_line_low = np.zeros(self.ntot_lines)*self.flamunit
        self.flam_line_hig = np.zeros(self.ntot_lines)*self.flamunit
        self.sigma_v_low   = np.zeros(self.ntot_lines)*self.vunit
        self.sigma_v_hig   = np.zeros(self.ntot_lines)*self.vunit

        j1 = 0
        for k,line in enumerate(self.multi_line):
            j0 = j1
            j1 += line.npar_line
            x_line_use = line.line_par_translator(x_line[j0:j1])
            for i in range(line.nlines):
                dv, flam_line, sigma_v = line.line_par_parser(i,x_line_use)
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


    #####

    #Check if fit can be run. If it can be run for each of the lines, it can be run.
    def can_fit_be_run(self,spec):
        if self.n_multilines==0:
            print("Please add an emission line before running.")
            return False
        can_all = True
        for line in self.multi_line:
            can_all = can_all and line.can_fit_be_run(spec)
        return can_all
