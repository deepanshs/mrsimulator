# -*- coding: utf-8 -*-
"""
Solid-state line-shape simulation module for Python
===================================================

mrsimulator is a very fast, real-time solid-state NMR line-shape simulation package,
capable of simulating line-shapes from both crystalline and amorphous-like materials.


It aims to provide simple and efficient solutions to the solid-state NMR line-shape
simulation problem. It includes tools for users to create their model spin-systems,
simulate and compare the line-shapes with the measurement, and perform least-squares
minimization, using a collection of pre-defined NMR methods.

See https://mrsimulator.readthedocs.io/en/stable/ for complete documentation.
"""
# version has to be specified at the start.
__version__ = "0.3.0b1"

from .spin_system import Site  # lgtm [py/import-own-module]
from .spin_system import SpinSystem  # lgtm [py/import-own-module]
from .simulator import Simulator  # lgtm [py/import-own-module]
from .transition import Transition  # lgtm [py/import-own-module]
from .method.event import Event  # lgtm [py/import-own-module]
from .method.spectral_dimension import SpectralDimension  # lgtm [py/import-own-module]
from .method import Method  # lgtm [py/import-own-module]
