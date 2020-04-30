import numpy  as np

from multi_line import Multi_Line_fit
from line_class import Line_fit

#This will be its own subclass of Line_fit, actually, but will be made from combinations of Multi_Line_fit objects.

#We will defferiate between broad and narrow emission lines. In this model, all narrow emission lines have the same width and systemic velocity shift, while this can be different for each of the broad emission lines.

class Complex_Line(Line_fit):

    def __init__(self, line_name, spec=None):
        super(Complex_Line).__init__(line_name, spec)
        self.multi_line = []
        return

    ###
    #To be defined

    def set_initial_fit_values(self,spec=None):
        return

    def self.set_cont_pars():
        return

    def self.set_line_pars():
        return

    def self.get_i_cont():
        return

    def self.get_i_line():
        return

    @property
    def npar_line(self):
        return

    def meet_cont_constraints():
        return

    def meet_line_constraints():
        return

    def flam_model():
        return

    def cont_par_parser():
        return

    ###

    @property
    def n_multilines(self):
        return len(self.multi_line)

    def add_line(self,line_name,width_type=None):

        #Setup the emission line.
        self.multi_line.append(Multi_Line_fit(line_name))
        self.multi_line[-1].width_type = width_type

        #If broad, FWHM>=1000 km/s. If narrow, FWHM<=1000 km/s
        if width_type == 'broad':
            self.multi_line[-1].sigma_v_min = 1000.*u.km/u.s/(2.*(2.*np.log(2.))**0.5)
        elif width_type == 'narrow':
            self.multi_line[-1].sigma_v_max = 1000.*u.km/u.s/(2.*(2.*np.log(2.))**0.5)
        elif width_type is not None:
            print("Unknown emission line type: ",width_type)
            print("Assuming full bounds for line ",line_name)

        return


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



        return

    #For now, we define a linear continuum.
    def flam_cont_model(lam,x_cont=None):
        a, b = self.cont_par_parser(x_cont)
        return a*lam+b



    #Check if fit can be run. If it can be run for each of the lines, it can be run.
    def can_fit_be_run(self,spec):
        if self.n_multilines==0:
            print("Please add an emission line before running.")
            return False
        can_all = True
        for line in self.multi_line:
            can_all = can_all and line.can_fit_be_run(spec)
        return can_all
