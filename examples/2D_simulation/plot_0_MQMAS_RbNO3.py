#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
RbNO3, 87Rb (I=3/2) MQMAS
^^^^^^^^^^^^^^^^^^^^^^^^^

87Rb (I=3/2) multi-quantum magic-angle (MQ-MAS) simulation.
"""
# %%
# The following is an MQMAS simulation of :math:`\text{RbNO}_3`, which has three
# distinct :math:`^{87}\text{Rb}` sites. The :math:`^{87}\text{Rb}` tensor parameters
# were obtained from Massiot `et. al.` [#f1]_
import matplotlib as mpl
import matplotlib.pyplot as plt
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.affine_transformation as af
import mrsimulator.signal_processing.apodization as apo
import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import MQVAS

# global plot configuration
font = {"size": 9}
mpl.rc("font", **font)
mpl.rcParams["figure.figsize"] = [4.25, 3.0]
# sphinx_gallery_thumbnail_number = 2

# %%
# **Step 1:** Create sites and spin systems.
Rb87_1 = Site(
    isotope="87Rb-1",
    isotropic_chemical_shift=-27.4,  # in ppm
    quadrupolar={"Cq": 1.68e6, "eta": 0.2},  # Cq is in Hz
)
Rb87_2 = Site(
    isotope="87Rb-2",
    isotropic_chemical_shift=-28.5,  # in ppm
    quadrupolar={"Cq": 1.94e6, "eta": 1.0},  # Cq is in Hz
)
Rb87_3 = Site(
    isotope="87Rb-3",
    isotropic_chemical_shift=-31.3,  # in ppm
    quadrupolar={"Cq": 1.72e6, "eta": 0.5},  # Cq is in Hz
)

sites = [Rb87_1, Rb87_2, Rb87_3]  # all sites
spin_systems = [SpinSystem(sites=[s]) for s in sites]

# %%
# **Step 2:** Create a Multi Quantum variable angle spinning method.
method = MQVAS(
    channels=["87Rb"],
    magnetic_flux_density=7,  # in T
    rotor_angle=54.735 * np.pi / 180,  # in radians
    spectral_dimensions=[
        {
            "count": 1024,
            "spectral_width": 2e4,  # in Hz
            "reference_offset": -5e3,  # in Hz
            "label": "3Q dimension - 1",
        },
        {
            "count": 1024,
            "spectral_width": 1e4,  # in Hz
            "reference_offset": -4e3,  # in Hz
            "label": "MAS dimension - 0",
        },
    ],
)

# %%
# **Step 3:** Create the Simulator object, add the method and spin system objects, and
# run the simulation
sim = Simulator()
sim.spin_systems += spin_systems  # add the spin systems
sim.methods += [method]  # add the method.

sim.run()

# %%
# **Step 4:**
# The plot of the simulation.
data = sim.methods[0].simulation
ax = plt.subplot(projection="csdm")
ax.imshow(data / data.max(), aspect="auto", cmap="gist_ncar_r")
ax.invert_xaxis()
ax.invert_yaxis()
plt.tight_layout()
plt.show()

# %%
# **Step 5:** Add post-simulation signal processing.
#
# Note, the above spectrum is a correlation between the triple quantum and MAS
# dimension, that is, the spectrum as acquired. To obtain an isotropic vs. MAS spectrum,
# we need to apply skew, and scaling transformations. For spin 3/2 undergoing a triple
# quantum excitation, apply a shear transformation parallel to the MAS dimension with
# a shear factor of 21/27 and then scale the 3Q dimension by 27/48. The following signal
# processing scheme first applies a shear and scaling transformation, followed by the
# line-broadening convolutions. Note, the MAS and 3Q dimensions are at index 0 and 1,
# respectively.
processor = sp.SignalProcessor(
    operations=[
        # shear parallel to MAS dimension.
        sp.IFFT(dim_index=1),
        af.Shear(factor=-21 / 27, dim_index=1, normal=0),
        sp.FFT(dim_index=1),
        # scale the 3Q-dimension
        af.Scale(factor=27 / 48, dim_index=1),
        # Gaussian convolution along both dimensions.
        sp.IFFT(dim_index=(0, 1)),
        apo.Gaussian(FWHM="0.2 kHz", dim_index=0),
        apo.Gaussian(FWHM="0.2 kHz", dim_index=1),
        sp.FFT(dim_index=(0, 1)),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation)
processed_data /= processed_data.max()

# %%
# **Step 6:** The plot of the simulation after signal processing.
ax = plt.subplot(projection="csdm")
ax.imshow(processed_data.real, vmax=1, cmap="gist_ncar_r", aspect="auto")
ax.set_ylabel("Isotropic dimension (ppm from 1M RbNO$_3$)")
ax.set_xlabel("MAS dimension (ppm from 1M RbNO$_3$)")
ax.set_xlim(-15, -70)
ax.set_ylim(-35, -65)
plt.tight_layout()
plt.show()

# %%
# .. [#f1] Massiot, D., Touzoa, B., Trumeaua, D., Coutures, J.P., Virlet, J., Florian,
#       P., Grandinetti, P.J. Two-dimensional magic-angle spinning isotropic
#       reconstruction sequences for quadrupolar nuclei, ssnmr, (1996), **6**, *1*,
#       73-83. `DOI: 10.1016/0926-2040(95)01210-9
#       <https://doi.org/10.1016/0926-2040(95)01210-9>`_
