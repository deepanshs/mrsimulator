#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Rb₂CrO₄, ⁸⁷Rb (I=3/2) COASTER
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

⁸⁷Rb (I=3/2) Correlation of anisotropies separated through echo refocusing (COASTER)
simulation.
"""
# %%
# The following is a Correlation of Anisotropies Separated Through Echo Refocusing
# (COASTER) simulation of :math:`\text{Rb}_2\text{CrO}_4`. The Rb site with the smaller
# quadrupolar interaction is selectively observed and reported by Ash `et al.` [#f1]_.
# The following is the simulation based on the published tensor parameters.
import matplotlib.pyplot as plt
import numpy as np

from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator import signal_processing as sp
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator.method import Method, SpectralDimension, SpectralEvent, MixingEvent

# sphinx_gallery_thumbnail_number = 3

# %%
# Generate the site and spin system objects.
site = Site(
    isotope="87Rb",
    isotropic_chemical_shift=-9,  # in ppm
    shielding_symmetric=SymmetricTensor(zeta=110, eta=0),
    quadrupolar=SymmetricTensor(
        Cq=3.5e6,  # in Hz
        eta=0.36,
        alpha=0,  # in rads
        beta=70 * 3.14159 / 180,  # in rads
        gamma=0,  # in rads
    ),
)
spin_system = SpinSystem(sites=[site])

# %%
# Use the generic method, `Method`, to simulate a 2D COASTER spectrum by customizing
# the method parameters, as shown below.
#
# By default, all transitions selected from one `SpectralEvent` will connect to all
# transitions selected by another `SpectralEvent` if there is no `MixingEvent` between
# defined them. Here we define a `MixingEvent` with an angle of 109.5 degrees to 
# selectively connect the three quantum and 1 quantum coherences.
coaster = Method(
    name="COASTER",
    channels=["87Rb"],
    magnetic_flux_density=9.4,  # in T
    rotor_angle=70.12 * np.pi / 180,  # in rads
    rotor_frequency=np.inf,
    spectral_dimensions=[
        SpectralDimension(
            count=512,
            spectral_width=4e4,  # in Hz
            reference_offset=-8e3,  # in Hz
            label="$\\omega_1$ (CSA)",
            events=[
                SpectralEvent(transition_query=[{"ch1": {"P": [3], "D": [0]}}]),
                MixingEvent(query={"ch1": {"angle": np.pi * 109.5 / 180, "phase": 0}}),
            ],
        ),
        # The last spectral dimension block is the direct-dimension
        SpectralDimension(
            count=512,
            spectral_width=8e3,  # in Hz
            reference_offset=-4e3,  # in Hz
            label="$\\omega_2$ (Q)",
            events=[SpectralEvent(transition_query=[{"ch1": {"P": [-1], "D": [0]}}])],
        ),
    ],
    affine_matrix=[[1, 0], [1 / 4, 3 / 4]],
)

# A graphical representation of the method object.
plt.figure(figsize=(5, 2.75))
coaster.plot()
plt.show()

# %%
# Create the Simulator object, add the method and spin system objects, and
# run the simulation.
sim = Simulator()
sim.spin_systems = [spin_system]  # add the spin systems
sim.methods = [coaster]  # add the method.

# configure the simulator object. For non-coincidental tensors, set the value of the
# `integration_volume` attribute to `hemisphere`.
sim.config.integration_volume = "hemisphere"
sim.run()

# %%
# The plot of the simulation.
data = sim.methods[0].simulation

plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
cb = ax.imshow(data.real / data.real.max(), aspect="auto", cmap="gist_ncar_r")
plt.colorbar(cb)
ax.set_xlim(-0, -55)
ax.set_ylim(80, -180)
plt.tight_layout()
plt.show()

# %%
# Add post-simulation signal processing.
processor = sp.SignalProcessor(
    operations=[
        # Gaussian convolution along both dimensions.
        sp.IFFT(dim_index=(0, 1)),
        sp.apodization.Gaussian(FWHM="0.15 kHz", dim_index=0),
        sp.apodization.Gaussian(FWHM="0.15 kHz", dim_index=1),
        sp.FFT(dim_index=(0, 1)),
    ]
)
processed_data = processor.apply_operations(data=data)
processed_data /= processed_data.max()

# %%
# The plot of the simulation after signal processing.
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
cb = ax.imshow(processed_data.real, cmap="gist_ncar_r", aspect="auto")
plt.colorbar(cb)
ax.set_xlim(-0, -55)
ax.set_ylim(80, -180)
plt.tight_layout()
plt.show()

# %%
# .. [#f1] Jason T. Ash, Nicole M. Trease, and Philip J. Grandinetti. Separating
#       Chemical Shift and Quadrupolar Anisotropies via Multiple-Quantum NMR
#       Spectroscopy, J. Am. Chem. Soc. (2008) **130**, 10858-10859.
#       `DOI: 10.1021/ja802865x <https://doi.org/10.1021/ja802865x>`_
