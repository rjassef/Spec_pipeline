import numpy as np
import astropy.units as u

from Fake_Spec import Fake_Spec

class Fake_LRIS_Spec(Fake_Spec):

    def __init__(self,_side,z=2.5,gmag=24.0,rmag=24.0):

        self.RT    = 5.*u.m       #Telescope radius
        self.texp  = 15.*u.minute #Exposure time.

        #Sky is super bright in the LRIS template made by
        #Dan. Something is off by a systematic factor of about
        #1000. This factor corrects for it:
        flux_sky_scale = 1.e-3

        self.side = _side
        if self.side!='blue' and self.side!='red':
            print("Please declare red or blue side")
            return

        self.set_side(gmag,rmag)
        super(Fake_LRIS_Spec,self).__init__(z=z,flux_sky_scale=flux_sky_scale)

        return

    
    def set_side(self,gmag,rmag):

        if self.side=='blue':
            self.lam_spec_min = 3200.*u.AA
            self.lam_spec_max = 5500.*u.AA
            self.RON          = 3.82
            self.mag_norm     = gmag
            self.lam_eff_norm = 4770.*u.AA
            self.sky_fname    = "../Spec_pipeline/Sky_Templates/"+\
                                "template_sky_LRIS_b.dat"
            self.sens_fname   = "../Spec_pipeline/Sensitivity_Files/"+\
                                "Sens_LRIS_B600.txt"
        else:
            self.lam_spec_min =  5500.*u.AA
            self.lam_spec_max = 10000.*u.AA
            self.RON          = 4.64
            self.mag_norm     = rmag
            self.lam_eff_norm = 6231.*u.AA
            self.sky_fname    = "../Spec_pipeline/Sky_Templates/"+\
                                "template_sky_LRIS_r.dat"
            self.sens_fname   = "../Spec_pipeline/Sensitivity_Files/"+\
                                "Sens_LRIS_R400.txt"
        return

