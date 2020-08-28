#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Albite, 27Al (I=5/2) MQMAS
^^^^^^^^^^^^^^^^^^^^^^^^^^

27Al (I=5/2) triple-quantum magic-angle (3Q-MAS) simulation.
"""
# %%
# The following is an :math:`^{27}\text{Al}` MQMAS simulation of albite
# :math:`\text{NaSi}_3\text{AlO}_8`. The :math:`^{87}\text{Rb}` tensor parameters were
# obtained from Massiot `et. al.` [#f1]_.
import matplotlib as mpl
import matplotlib.pyplot as plt
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import ThreeQ_MAS

# global plot configuration
font = {"size": 9}
mpl.rc("font", **font)
mpl.rcParams["figure.figsize"] = [4.25, 3.0]
# sphinx_gallery_thumbnail_number = 2

# %%
# **Step 1:** Create the site and spin system.
site = Site(
    isotope="27Al",
    isotropic_chemical_shift=64.7,  # in ppm
    quadrupolar={"Cq": 3.25e6, "eta": 0.68},  # Cq is in Hz
)

spin_systems = [SpinSystem(sites=[site])]

# %%
# **Step 2:** Create a Triple Quantum magic-angle spinning method.
method = ThreeQ_MAS(
    channels=["27Al"],
    magnetic_flux_density=7,  # in T
    spectral_dimensions=[
        {
            "count": 512,
            "spectral_width": 1e4,  # in Hz
            "reference_offset": -2.5e3,  # in Hz
            "label": "Isotropic dimension",
        },
        {
            "count": 512,
            "spectral_width": 1e4,  # in Hz
            "reference_offset": 4e3,  # in Hz
            "label": "MAS dimension",
        },
    ],
)

# %%
# **Step 3:** Create the Simulator object, add the method and spin system objects, and
# run the simulation
sim = Simulator()
sim.spin_systems = spin_systems  # add the spin systems
sim.methods = [method]  # add the method.
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
ax.imshow(processed_data.real, cmap="gist_ncar_r", aspect="auto")
ax.set_xlim(75, 25)
ax.set_ylim(-15, -65)
plt.tight_layout()
plt.show()

# %%
# .. [#f1] Massiot, D., Touzoa, B., Trumeaua, D., Coutures, J.P., Virlet, J., Florian,
#       P., Grandinetti, P.J. Two-dimensional magic-angle spinning isotropic
#       reconstruction sequences for quadrupolar nuclei, ssnmr, (1996), **6**, *1*,
#       73-83. `DOI: 10.1016/0926-2040(95)01210-9
#       <https://doi.org/10.1016/0926-2040(95)01210-9>`_
