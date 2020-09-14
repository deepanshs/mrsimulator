#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Rb2CrO4, 87Rb (I=3/2) SAS
^^^^^^^^^^^^^^^^^^^^^^^^^

87Rb (I=3/2) Switched-angle spinning (SAS) simulation.
"""
# %%
# The following is a switched-angle spinning (SAS) simulation of
# :math:`\text{Rb}_2\text{CrO}_4`. While :math:`\text{Rb}_2\text{CrO}_4` has two
# rubidium sites, the site with the smaller quadrupolar interaction was selectively
# observed and reported by Shore `et. al.` [#f1]_. The following is the simulation
# based on the published tensor parameters.
import matplotlib as mpl
import matplotlib.pyplot as plt
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo
import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import Method2D

# global plot configuration
font = {"size": 9}
mpl.rc("font", **font)
mpl.rcParams["figure.figsize"] = [4.25, 3.0]
# sphinx_gallery_thumbnail_number = 2

# %%
# **Step 1:** Create the site and the spin systems.
site = Site(
    isotope="87Rb",
    isotropic_chemical_shift=-7,  # in ppm
    shielding_symmetric={"zeta": 110, "eta": 0},
    quadrupolar={
        "Cq": 3.5e6,  # in Hz
        "eta": 0.3,
        "alpha": 0,  # in rads
        "beta": 70 * np.pi / 180,  # in rads
        "gamma": 0,  # in rads
    },
)
spin_system = SpinSystem(sites=[site])

# %%
# **Step 2:** Create a generic 2D method, `Method2D`, method. To simulate a SAS
# spectrum, customize the method parameters, as shown below. Note, the Method2D method
# simulates an infinite spinning speed spectrum.
method = Method2D(
    channels=["87Rb"],
    magnetic_flux_density=4.2,  # in T
    spectral_dimensions=[
        {
            "count": 256,
            "spectral_width": 1.5e4,  # in Hz
            "reference_offset": -5e3,  # in Hz
            "label": "70.12 dimension",
            "events": [{"rotor_angle": 70.12 * np.pi / 180}],  # in radians
        },
        {
            "count": 512,
            "spectral_width": 15e3,  # in Hz
            "reference_offset": -7e3,  # in Hz
            "label": "MAS dimension",
            "events": [{"rotor_angle": 54.74 * np.pi / 180}],  # in radians
        },
    ],
)

# %%
# **Step 3:** Create the Simulator object, add the method and spin system objects, and
# run the simulation.
sim = Simulator()
sim.spin_systems = [spin_system]  # add the spin systems
sim.methods = [method]  # add the method.

# configure the simulator object
sim.config.integration_volume = "hemisphere"
sim.run()

# %%
# **Step 4:**
# The plot of the simulation.
data = sim.methods[0].simulation
ax = plt.subplot(projection="csdm")
ax.imshow(data / data.max(), aspect="auto", cmap="gist_ncar_r")
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# **Step 5:** Add post-simulation signal processing.
#
processor = sp.SignalProcessor(
    operations=[
        # Gaussian convolution along both dimensions.
        sp.IFFT(dim_index=(0, 1)),
        apo.Gaussian(FWHM="0.2 kHz", dim_index=0),
        apo.Gaussian(FWHM="0.2 kHz", dim_index=1),
        sp.FFT(dim_index=(0, 1)),
    ]
)
processed_data = processor.apply_operations(data=data)
processed_data /= processed_data.max()

# %%
# **Step 6:** The plot of the simulation after signal processing.
ax = plt.subplot(projection="csdm")
ax.imshow(processed_data.real, cmap="gist_ncar_r", aspect="auto")
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# .. [#f1] Shore, J.S., Wang, S.H., Taylor, R.E., Bell, A.T., Pines, A. Determination of
#       quadrupolar and chemical shielding tensors using solid-state two-dimensional NMR
#       spectroscopy, J. Chem. Phys. (1996)  **105** *21*, 9412.
#       `DOI: 10.1063/1.472776 <https://doi.org/10.1063/1.472776>`_
