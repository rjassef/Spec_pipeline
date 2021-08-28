name = "Spec_pipeline"

#First, write the standard Spec_Reader modules we'll be using.
from .Spec_Reader.LRIS_Spec import LRIS_Spec
from .Spec_Reader.GMOS_Spec import GMOS_Spec
from .Spec_Reader.DBSP_Spec import DBSP_Spec
from .Spec_Reader.SDSS_Spec import SDSS_Spec
from .Spec_Reader.DEIMOS_Spec import DEIMOS_Spec
#from .Spec_Reader.stern_plot import stern_plot
#from .Stern_plots.stern_plot import stern_plot
from .Stern_plots.stern_plot_v2 import stern_plot #as stern_plot_new

#Now, from the Line_Fitter part, import only the multi_line class modules, the rest are deprecated.
from .Line_Fitter.multi_line import Multi_Line_fit
from .Line_Fitter.complex_line import Complex_Line_fit
from .Line_Fitter.Powc_Line_fit import Powc_Line_fit
