#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step Apodization
^^^^^^^^^^^^^^^^
"""
# %%
# In this example, we will use the
# :py:class:`~mrsimulator.signal_processing.apodization.Step` class to apply a pointwise
# step apodization on the Fourier transform of an example dataset. The step function
# is defined as follows
#
# .. math::
#
#   f(x) = \begin{cases}
#     1, \texttt{rising_edge} \leq x \leq \texttt{falling_edge} \\
#     0, \text{otherwise}
#     \end{cases}
#
# where ``rising_edge`` is the start of the step function window and ``falling_edge``
# is the end of the step function window
#
# Below we import the necessary modules
import csdmpy as cp
import matplotlib.pyplot as plt
import numpy as np
from mrsimulator import signal_processing as sp

# sphinx_gallery_thumbnail_number = 1

# %%
# First we create ``processor``, and instance of the
# :py:class:`~mrsimulator.signal_processing.SignalProcessor` class. The required
# attribute of the SignalProcessor class, *operations*, is a list of operations to which
# we add a :py:class:`~mrsimulator.signal_processing.apodization.Step` object
# sandwitched between two Fourier transformations. Here the step window is between
# -0.01 and 0.01 seconds.
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.Step(rising_edge="-0.01 s", falling_edge="0.01 s"),
        sp.FFT(),
    ]
)

# %%
for i in range(10):
    print(i)

# %%
# Next we create a CSDM object with a test dataset which our signal processor will
# operate on. Here, the dataset is a delta function centered at 250 Hz with a some
# applied Gaussian line broadening.
test_data = np.zeros(500)
test_data[250] = 1
csdm_object = cp.CSDM(
    dependent_variables=[cp.as_dependent_variable(test_data)],
    dimensions=[cp.LinearDimension(count=500, increment="1 Hz")],
)

# Create processor to apply line broadening
pre_processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.Gaussian(FWHM="50 Hz"),
        sp.FFT(),
    ]
)

# Apply Gaussian line broadening
pre_processed_data = pre_processor.apply_operations(data=csdm_object)

# %%
# To apply the previously defined signal processor, we use the
# :py:meth:`~mrsimulator.signal_processing.SignalProcessor.apply_operations` method as
# as follows
processed_data = processor.apply_operations(data=pre_processed_data.copy())

# %%
# To see the results of the step apodization, we create a simple plot using the
# ``matplotlob`` library.
fig, ax = plt.subplots(1, 2, figsize=(8, 3), subplot_kw={"projection": "csdm"})
ax[0].plot(pre_processed_data, color="black", linewidth=1)
ax[0].set_title("Before")
ax[1].plot(processed_data.real, color="black", linewidth=1)
ax[1].set_title("After")
plt.tight_layout()
plt.show()
