name = "Spec_pipeline"

#First, write the standard Spec_Reader modules we'll be using.
from .Spec_Reader.LRIS_Spec import LRIS_Spec
from .Spec_Reader.GMOS_Spec import GMOS_Spec
from .Spec_Reader.DBSP_Spec import DBSP_Spec
from .Spec_Reader.SDSS_Spec import SDSS_Spec
from .Spec_Reader.DEIMOS_Spec import DEIMOS_Spec

#Now, from the Line_Fitter part, import only the multi_line class modules, the rest are deprecated.
from .Line_Fitter.multi_line import Multi_Line_fit
from .Line_Fitter.complex_line import Complex_Line_fit
