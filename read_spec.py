#This is just a wrapper function for reading the spectra. 

import sys
sys.path.append("Spec_Reader/")
from GMOS_Spec import GMOS_Spec
from LRIS_Spec import LRIS_Spec
from DBSP_Spec import DBSP_Spec

def read_spec(name, zspec, instrument, fits_files):
    
    if instrument=="GMOS":
        spec = GMOS_Spec(name,zspec,fits_files)
    elif instrument=="LRIS":
        spec = LRIS_Spec(name,zspec,fits_files)
    elif instrument=="DBSP":
        spec = DBSP_Spec(name,zspec,fits_files)
    else:
        print("Unknown instrument {0:s}".format(instrument))
        return None
    return spec

