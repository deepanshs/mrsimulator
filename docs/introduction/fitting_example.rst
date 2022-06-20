.. _fitting_example:

Least-Squares Fitting Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


In this example, we illustrate the use of the mrsimulator objects to

- create a quadrupolar fitting model using Simulator and SignalProcessor objects,
- use the fitting model to perform a least-squares analysis, and
- extract the fitting parameters from the model.

We use the `LMFIT <https://lmfit.github.io/lmfit-py/>`_ library to fit the spectrum.
The following example shows the least-squares fitting procedure applied to the
:math:`^{17}\text{O}` MAS NMR spectrum of :math:`\text{Na}_{2}\text{SiO}_{3}` [#f5]_.

Start by importing the relevant modules.

.. plot::
    :context: close-figs

    import csdmpy as cp
    import matplotlib.pyplot as plt
    from lmfit import Minimizer

    from mrsimulator import Simulator, SpinSystem, Site
    from mrsimulator.method.lib import BlochDecayCTSpectrum
    from mrsimulator import signal_processing as sp
    from mrsimulator.utils import spectral_fitting as sf
    from mrsimulator.utils import get_spectral_dimensions
    from mrsimulator.spin_system.tensors import SymmetricTensor






