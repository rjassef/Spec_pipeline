#This is just a wrapper function for reading the spectra. 

import sys
sys.path.append("Spec_Reader/")
from GMOS_Spec import GMOS_Spec
from LRIS_Spec import LRIS_Spec
from DBSP_Spec import DBSP_Spec

def read_spec(name, zspec, instrument, fits_files, line_center=None, 
              blue=False,red=False):
    
    if instrument=="GMOS":
        spec = GMOS_Spec(name,zspec,fits_files,line_center)
    elif instrument=="LRIS":
        if line_center is None and not blue and not red:
            print("Cannot read spectrum.")
            print("Provide a side or a wavelength of interest")
            return None
        spec = LRIS_Spec(name,zspec,fits_files,line_center,blue,red)
    elif instrument=="DBSP":
        if line_center is None and not blue and not red:
            print("Cannot read spectrum.")
            print("Provide a side or a wavelength of interest")
            return None
        spec = DBSP_Spec(name,zspec,fits_files,line_center,blue,red)
    else:
        print("Unknown instrument {0:s}".format(instrument))
        return None
    return spec

