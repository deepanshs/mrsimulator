#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Skewed Gaussian Apodization
^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# In this example, we will use the
# :py:class:`~mrsimulator.signal_processing.apodization.SkewedGaussian` class to
# apply a apodization on the Foruier transform of an example dataset. The
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
        sp.apodization.SkewedGaussian(skew=2, FWHM="100 s", dv_index=0),
        sp.apodization.SkewedGaussian(skew=-3.5, FWHM="100 s", dv_index=1),
        sp.FFT(),
    ]
)

# %%
# Next we create a CSDM object with a test dataset which our signal processor will
# operate on. Here, the dataset consists of two dependent variables each
# spaning 500 seconds with a delta function centered at
# 250 seconds.
test_data = np.zeros(500)
test_data[250] = 1
csdm_object = cp.CSDM(
    dependent_variables=[
        cp.as_dependent_variable(test_data),
        cp.as_dependent_variable(test_data),
    ],
    dimensions=[cp.LinearDimension(count=500, increment="1 s")],
)
# set the labels for dependent variables
csdm_object.y[0].name = "skew=2"
csdm_object.y[1].name = "skew=-3.5"

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
ax[0].plot(csdm_object, color="black", linewidth=1)
ax[0].set_title("Before")
ax[1].plot(processed_data.real, linewidth=1)
ax[1].set_title("After")
plt.tight_layout()
plt.show()
