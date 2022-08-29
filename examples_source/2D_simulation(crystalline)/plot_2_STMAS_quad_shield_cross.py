#!/usr/bin/env python
"""
RbNO₃, ⁸⁷Rb (I=3/2) STMAS
^^^^^^^^^^^^^^^^^^^^^^^^^

⁸⁷Rb (I=3/2) satellite-transition off magic-angle spinning simulation.
"""
# %%
# The following is an example of the STMAS simulation of :math:`\text{RbNO}_3`. The
# :math:`^{87}\text{Rb}` tensor parameters were obtained from Massiot `et al.` [#f1]_.
import matplotlib.pyplot as plt

from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.method.lib import ST1_VAS
from mrsimulator import signal_processor as sp
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator.method import SpectralDimension

# sphinx_gallery_thumbnail_number = 2

# %%
# Generate the site and spin system objects.
Co_site = Site(
    isotope="59Co",  # 59Co
    isotropic_chemical_shift=0,  # in ppm
    shielding_symmetric=SymmetricTensor(zeta=-1750, eta=0.2),
    quadrupolar=SymmetricTensor(Cq=3.1e6, eta=0),  # Cq is in Hz
)

spin_system = SpinSystem(sites=[Co_site])

# %%
# Select a satellite-transition variable-angle spinning method. The
# following `ST1_VAS` method correlates the frequencies from the two inner-satellite
# transitions to the central transition. Note, STMAS measurements are highly suspectable
# to rotor angle mismatch. In the following, we show two methods, the first at the
# magic angle and second deliberately miss-sets by approximately 0.0059 degrees.

method = ST1_VAS(
    channels=["59Co"],
    magnetic_flux_density=4.684,  # in T
    rotor_angle=54.7359 * 3.14159 / 180,  # in rad (magic angle)
    spectral_dimensions=[
        SpectralDimension(
            count=256,
            spectral_width=2e3,  # in Hz
            # reference_offset=-1e3,  # in Hz
            label="Isotropic dimension",
        ),
        SpectralDimension(
            count=512,
            spectral_width=2e3,  # in Hz
            reference_offset=-1e3,  # in Hz
            label="MAS dimension",
        ),
    ],
)

# A graphical representation of the method object.
plt.figure(figsize=(5, 2.5))
method.plot()
plt.show()

# %%
# Create the Simulator object, add the method and spin system objects, and
# run the simulation.
sim = Simulator(spin_systems=[spin_system], methods=[method])
sim.run()

# %%
# Add post-simulation signal processing.
dataset = sim.methods[0].simulation
processor = sp.SignalProcessor(
    operations=[
        # Gaussian convolution along both dimensions.
        sp.IFFT(dim_index=(0, 1)),
        sp.apodization.Gaussian(FWHM="20 Hz", dim_index=0),
        sp.apodization.Gaussian(FWHM="20 Hz", dim_index=1),
        sp.FFT(dim_index=(0, 1)),
    ]
)
processed_dataset = processor.apply_operations(dataset=dataset)
processed_dataset /= processed_dataset.max()

# %%
# The plot of the simulation.
_ = [item.to("kHz", "nmr_frequency_ratio") for item in processed_dataset.x]

plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
cb = ax.imshow(processed_dataset.real, cmap="gist_ncar_r", aspect="auto")
plt.colorbar(cb)
ax.invert_xaxis()
ax.invert_yaxis()
plt.tight_layout()
plt.show()
