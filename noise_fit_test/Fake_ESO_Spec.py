import numpy as np
import astropy.units as u

from Fake_Spec import Fake_Spec

class Fake_ESO_Spec(Fake_Spec):

    def __init__(self,z=2.5,gmag=24.0,moon=0):

        self.RT    = 4.*u.m       #Telescope radius
        self.texp  = 15.*u.minute #Exposure time.

        self.set_side(gmag,moon)
        super(Fake_ESO_Spec,self).__init__(z=z)

        return

    
    def set_side(self,gmag,moon):

        self.lam_spec_min = 3500.*u.AA
        self.lam_spec_max = 9500.*u.AA
        self.RON          = 3.82
        self.mag_norm     = gmag
        self.lam_eff_norm = 4770.*u.AA
        self.sky_fname    = "ESO_Sky_Models/template_sky_ESO_{0:d}.dat".format(moon)
        self.sens_fname   = "ESO_Sky_Models/Sens_ESO.dat"

        return
