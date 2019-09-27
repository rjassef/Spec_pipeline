#!/usr/bin/env python

import warnings
warnings.simplefilter("ignore")

import numpy as np
import astropy.units as u
from astropy.constants import h,c

import sys
import os
sys.path.append(os.environ['SPEC_PIPE_LOC'])
from Spec_pipeline.Spec_Reader.Spec import Spec

#Start by creating a fake spectrum.
fake = Spec("Fake",3.0)

#Set the wavelength range.
fake.lam_obs = np.arange(3800.,7000.,2.)*u.AA


