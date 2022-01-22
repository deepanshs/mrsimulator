#!/usr/bin/env python
# -*- coding: utf-8 -*-
# NOTE: Example removed until skewed gaussian apodization resolved
"""
Skewed Gaussian Apodization
^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# In this example, we will use the
# :py:class:`~mrsimulator.signal_processing.apodization.SkewedGaussian` class to
# apply a Skewed Gaussian convolution on an example dataset. The
# skewed Gaussian function is defined as follows
#
# .. math::
#
#   f(x, \sigma, \alpha) = 2\phi(x, \sigma)\Phi(\alpha x, \sigma)
#
# where :math:`\phi` is a normal PDF and :math:`\Phi` is
# a normal CDF both with standard deviation
# :math:`\sigma = \frac{\text{FWHM}}{2 \sqrt{2 \ln{2}}}` and
# :math:`\alpha` is the skewness parameter.
#
# Below we import the necessary modules
import csdmpy as cp
import numpy as np
from mrsimulator import signal_processing as sp

# sphinx_gallery_thumbnail_number = 1

# %%
# First we create ``processor``, and instance of the
# :py:class:`~mrsimulator.signal_processing.SignalProcessor` class. The required
# attribute of the SignalProcessor class, *operations*, is a list of operations to which
# we add a :py:class:`~mrsimulator.signal_processing.apodization.SkewedGaussian` object
# sandwitched between two Fourier transformations. Here :math:`\text{skew} = \alpha = 2`
# for the first dependent variable and :math:`\text{skew} = \alpha = -3.5`
# and the FWHM is 100 seconds.
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.SkewedGaussian(skew=2, FWHM="100 Hz", dv_index=0),
        sp.apodization.SkewedGaussian(skew=-3.5, FWHM="100 Hz", dv_index=1),
        sp.FFT(),
    ]
)

# %%
# Next we create a CSDM object with a test dataset which our signal processor will
# operate on. Here, the dataset consists of two dependent variables each
# spaning 500 seconds with delta functions centered at -100 and
# 100 Hz, respectively.
test_data_0 = np.zeros(500)
test_data_0[150] = 1
test_data_1 = np.zeros(500)
test_data_1[350] = 0.5
csdm_object = cp.CSDM(
    dependent_variables=[
        cp.as_dependent_variable(test_data_0),
        cp.as_dependent_variable(test_data_1),
    ],
    dimensions=[cp.LinearDimension(count=500, increment="1 Hz", complex_fft=True)],
)

# %%
# To apply the previously defined signal processor, we use the
# :py:meth:`~mrsimulator.signal_processing.SignalProcessor.apply_operations` method as
# as follows
processed_data = processor.apply_operations(data=csdm_object)

# %%
# To see the results of the skewed Gaussian apodization, we create a simple plot using
# the ``matplotlob`` library.
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=(8, 3.5), subplot_kw={"projection": "csdm"})
ax[0].plot(csdm_object, linewidth=1)
ax[0].set_title("Before")
ax[1].plot(processed_data.real, linewidth=1)
ax[1].set_title("After")
ax[1].legend(["skew=2", "skew=-3.5"])
plt.tight_layout()
plt.show()
