#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Lorentzian (Exponential) Apodization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# sphinx_gallery_thumbnail_number = 1
# %%
# In this example, we will use an exponential function to preform a pointwise Lorentzian
# apodization on the Fourier transform an example dataset. The exponential function
# used for this apodization is defined as follows
#
# .. math::
#
#   \usepackage{asmmath}
#   \begin{align*}
#   f(x) = e^{-\sigma \pi |x|}
#   \end{align*}
#
# where :math:`\sigma` is parametrized by the the full width at half maximum as follows
#
# .. math::
#
#   \sigma = \frac{\text{FWHM}}{2 \sqrt{2 \ln{2}}}
#
# Below we import the necessary modules
import csdmpy as cp
import numpy as np
from mrsimulator import signal_processing as sp

# First we create ``processor``, an instance of the
# :py:class:`~mrsimulator.signal_processing.SignalProcessor` class. The required
# attribute of the SignalProcessor class, *operations*, is a list of operations to which
# we add a :py:class:`~mrsimulator.signal_processing.apodization.Exponential` object
# sandwitched between two Fourier transformations.
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.Exponential(FWHM="100 s"),
        sp.FFT(),
    ]
)

# %%
# Next we create a CSDM object with a test dataset which our signal processor will
# operate on. Here, the dataset spans 500 seconds with a delta function centered at
# 250 seconds.
test_data = np.zeros(500)
test_data[250] = 1
csdm_object = cp.CSDM(
    dependent_variables=[cp.as_dependent_variable(test_data)],
    dimensions=[cp.LinearDimension(count=500, increment="1 s")],
)

# %%
# Now to apply the processor to the CSDM object, use the
# :py:meth:`~mrsimulator.signal_processing.SignalProcessor.apply_operations` method as
# follows
processed_data = processor.apply_operations(data=csdm_object.copy())

# %%
# To see the results of the exponential apodization, we create a simple plot using the
# ``matplotlib`` library.
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=(8, 3), subplot_kw={"projection": "csdm"})
ax[0].plot(csdm_object, color="black", linewidth=1)
ax[0].set_title("Before")
ax[1].plot(processed_data.real, color="black", linewidth=1)
ax[1].set_title("After")
plt.tight_layout()
plt.show()
