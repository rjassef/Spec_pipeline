Installation
************

Requirements
============

``Spec_pipeline`` requires a number of packages. The version number is the one the code has been tested with.

- `Python <https://www.python.org/>`_ 3.6.9
- `Numpy <https://numpy.org/>`_ 1.18.5
- `Astropy <https://www.astropy.org/>`_ 4.0.1
- `Matplotlib <https://matplotlib.org/>`_ 4.1.1

You should be able to obtain all of these packages trough `Astroconda <https://astroconda.readthedocs.io/en/latest/>`_ .

This documentation was written using `astropy-sphinx`. 

Installing Spec_pipeline
========================

Download the source from the `github <https://github.com/rjassef/Spec_pipeline/tree/v0.5>`_ repository and unpack.

Before using, you need to add the path to the folder to your `PYTHONPATH`. In bash and zsh, you can do this by executing the included script setup_bash in the root folder of the package run::

    source setup_bash

You can also add it to your .bashrc, .zshrc or similar by inluding the following lines run::

    SPEC_PIPE_LOC=path_to_Spec_pipeline
    export SPEC_PIPE_LOC

    PYTHONPATH=$PYTHONPATH:$SPEC_PIPE_LOC
    export PYTHONPATH


Testing the Installation
========================

TBA
