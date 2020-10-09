Quick Start
***********

Reading a Spectrum
==================

First, create a folder in which you will be working. We will call this folder `work-folder`. Inside that folder, create a directory called `data` and copy your fits files there.

Once you have done that, you can load your spectrum in a python shell. For this example, we will assume your object is called `W0123+4567` with redshift 1.234, and that the spectrum you want load corresponds to the blue arm of the LRIS spectrograph in the Keck telescope. Assuming the fits file with the spectrum is called `w0123p4567_b.f.fits`, you would load it as::

    >>> from Spec_pipeline import LRIS_Spec

    >>> spec = LRIS_Spec("W0123+4567", 1.234, "w0123p4567_b.f.fits", blue=True)

You can find a list of all the currently implemented instrument classes in :ref:`implemented-spec-classes`. Once the spectrum is loaded, an error spectrum is immediately created by following the process described in Eisenhardt et al. (2020).

Fitting an Emission Line to the Spectrum
========================================

One you have loaded an spectrum, you can fit an emission line to it. For this example, we will use the Multi_Line_fit module, which is capable of fitting simultaneous single Gaussian emission lines and a linear continuum to an spectrum. You can find details of the currently implemented fitting models at :ref:`implemented-line-fitters`.

To fit the `[OII]` doublet and plot the results, you would type::

    >>> from Spec_pipeline import Multi_Line_fit

    >>> oii_fitter = Multi_Line_fit("[OII]", spec=spec)

    >>> oii_fitter.run_fit()

    >>> oii_fitter.plot()

To obtain uncertainties in the parameters, you can also run a Monte Carlo process to produce nrep re-sampled spectra, fit them, and get the uncertainties from the dispersion of the chain of best-fit parameters. This is done by typing::

    >>> oii_fitter.run_MC(nrep=nrep)

    >>> oii_fitter.plot()

Note that this MC process runs by default in parallel mode using all available cores.
