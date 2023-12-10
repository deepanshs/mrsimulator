#!/usr/bin/env python
"""
Coupled spin-1/2 (MAS Quadrupolar-dipolar spectrum)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

¹³C-14N static dipolar coupling simulation.
"""
# %%
import matplotlib.pyplot as plt
import numpy as np

from mrsimulator import Simulator, SpinSystem, Site, Coupling
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator import signal_processor as sp
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator.method import SpectralDimension

# sphinx_gallery_thumbnail_number = 1

# %%
# Create a 13C-14N coupled spin system.
spin_system = SpinSystem(
    sites=[
        Site(isotope="13C", isotropic_chemical_shift=0.0),
        Site(
            isotope="14N",
            isotropic_chemical_shift=0,  # in ppm
            quadrupolar=SymmetricTensor(
                Cq=1.18e6,  # in Hz
                eta=0.54,
                alpha=0,
                beta=5 * np.pi / 180,
                gamma=0,
            ),
        ),
    ],
    couplings=[Coupling(site_index=[0, 1], dipolar=SymmetricTensor(D=-370))],
)
# %%
# Create a BlochDecaySpectrum method.
method = BlochDecaySpectrum(
    channels=["13C"],
    magnetic_flux_density=3.5,  # in T
    rotor_frequency=12000,  # in Hz
    spectral_dimensions=[SpectralDimension(count=2048, spectral_width=400)],
)

# %%
# Create the Simulator object and add the method and the spin system object.
sim = Simulator(spin_systems=[spin_system], methods=[method])
sim.run()

# %%
# Add post-simulation signal processing.
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.Gaussian(FWHM="1 Hz"),
        sp.FFT(),
    ]
)
processed_dataset = processor.apply_operations(dataset=sim.methods[0].simulation)

# %%
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(processed_dataset.real, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()
