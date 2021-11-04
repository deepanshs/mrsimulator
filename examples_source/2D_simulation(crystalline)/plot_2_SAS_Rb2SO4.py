#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Rb₂SO₄, ⁸⁷Rb (I=3/2) SAS
^^^^^^^^^^^^^^^^^^^^^^^^

⁸⁷Rb (I=3/2) Switched-angle spinning (SAS) simulation.
"""
# %%
# The following is an example of Switched-Angle Spinning (SAS) simulation of
# :math:`\text{Rb}_2\text{SO}_4`, which has two distinct rubidium sites. The NMR
# tensor parameters for these sites are taken from Shore `et al.` [#f1]_.
import matplotlib.pyplot as plt

from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.methods import Method2D
from mrsimulator import signal_processing as sp

# sphinx_gallery_thumbnail_number = 3

# %%
# Generate the site and spin system objects.
sites = [
    Site(
        isotope="87Rb",
        isotropic_chemical_shift=16,  # in ppm
        quadrupolar={"Cq": 5.3e6, "eta": 0.1},  # Cq in Hz
    ),
    Site(
        isotope="87Rb",
        isotropic_chemical_shift=40,  # in ppm
        quadrupolar={"Cq": 2.6e6, "eta": 1.0},  # Cq in Hz
    ),
]
spin_systems = [SpinSystem(sites=[s]) for s in sites]

# %%
# Use the generic 2D method, `Method2D`, to simulate a SAS spectrum by customizing the
# method parameters, as shown below. Note, the Method2D method simulates an infinite
# spinning speed spectrum.
sas = Method2D(
    name="Switched Angle Spinning",
    channels=["87Rb"],
    magnetic_flux_density=9.4,  # in T
    spectral_dimensions=[
        {
            "count": 256,
            "spectral_width": 3.5e4,  # in Hz
            "reference_offset": 1e3,  # in Hz
            "label": "90 dimension",
            "events": [
                {
                    "rotor_angle": 90 * 3.14159 / 180,
                    "transition_query": [{"P": [-1], "D": [0]}],
                }
            ],  # in radians
        },
        {
            "count": 256,
            "spectral_width": 22e3,  # in Hz
            "reference_offset": -4e3,  # in Hz
            "label": "MAS dimension",
            "events": [
                {
                    "rotor_angle": 54.74 * 3.14159 / 180,
                    "transition_query": [{"P": [-1], "D": [0]}],
                }
            ],  # in radians
        },
    ],
)

# A graphical representation of the method object.
plt.figure(figsize=(5, 3.5))
sas.plot()
plt.show()

# %%
# Create the Simulator object, add the method and spin system objects, and
# run the simulation.
sim = Simulator()
sim.spin_systems = spin_systems  # add the spin systems
sim.methods = [sas]  # add the method.
sim.run()

# %%
# The plot of the simulation.
data = sim.methods[0].simulation

plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
cb = ax.imshow(data.real / data.real.max(), aspect="auto", cmap="gist_ncar_r")
plt.colorbar(cb)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# Add post-simulation signal processing.
processor = sp.SignalProcessor(
    operations=[
        # Gaussian convolution along both dimensions.
        sp.IFFT(dim_index=(0, 1)),
        sp.apodization.Gaussian(FWHM="0.4 kHz", dim_index=0),
        sp.apodization.Gaussian(FWHM="0.4 kHz", dim_index=1),
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
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# .. [#f1] Shore, J.S., Wang, S.H., Taylor, R.E., Bell, A.T., Pines, A. Determination of
#       quadrupolar and chemical shielding tensors using solid-state two-dimensional NMR
#       spectroscopy, J. Chem. Phys. (1996)  **105** *21*, 9412.
#       `DOI: 10.1063/1.472776 <https://doi.org/10.1063/1.472776>`_
