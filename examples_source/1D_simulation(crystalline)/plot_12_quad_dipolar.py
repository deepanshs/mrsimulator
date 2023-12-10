#!/usr/bin/env python
"""
Influence of 14N on 13C NMR MAS spectra of glycine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The alpha-carbon resonance of glycine, 13C (I=1/2), attached to 14N (I=1).
The 14N quadrupolar tensor parameters were obtained from Hexem `et al.` [#f1]_
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
processed_dataset.dimensions[0].to("Hz")

# %%
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(processed_dataset.real, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()


# %%
#  .. [#f1] Hexem, J. G., Frey, M. H., and Opella, S. J., Influence of $^{14}\text{N}$ on $^{13}\text{C}$ NMR Spectra of Solids, J. Am. Chem. Soc., 1981, **103**, 224-226. [DOI: 10.1021/ja00391a057](https://doi.org/10.1021/ja00391a057)
#
