#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
RbNO3, 87Rb (I=3/2) 3QMAS
^^^^^^^^^^^^^^^^^^^^^^^^^

87Rb (I=3/2) triple-quantum magic-angle spinning (3Q-MAS) simulation.
"""
# %%
# The following is an example of the 3QMAS simulation of :math:`\text{RbNO}_3`, which
# has three distinct :math:`^{87}\text{Rb}` sites. The :math:`^{87}\text{Rb}` tensor
# parameters were obtained from Massiot `et. al.` [#f1]_.
import matplotlib as mpl
import matplotlib.pyplot as plt
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import ThreeQ_VAS

# global plot configuration
font = {"size": 9}
mpl.rc("font", **font)
mpl.rcParams["figure.figsize"] = [4.25, 3.0]
# sphinx_gallery_thumbnail_number = 2

# %%
# Generate the site and spin system objects.
Rb87_1 = Site(
    isotope="87Rb",
    isotropic_chemical_shift=-27.4,  # in ppm
    quadrupolar={"Cq": 1.68e6, "eta": 0.2},  # Cq is in Hz
)
Rb87_2 = Site(
    isotope="87Rb",
    isotropic_chemical_shift=-28.5,  # in ppm
    quadrupolar={"Cq": 1.94e6, "eta": 1.0},  # Cq is in Hz
)
Rb87_3 = Site(
    isotope="87Rb",
    isotropic_chemical_shift=-31.3,  # in ppm
    quadrupolar={"Cq": 1.72e6, "eta": 0.5},  # Cq is in Hz
)

sites = [Rb87_1, Rb87_2, Rb87_3]  # all sites
spin_systems = [SpinSystem(sites=[s]) for s in sites]

# %%
# Select a Triple Quantum variable-angle spinning method. You may optionally
# provide a `rotor_angle` to the method. The default `rotor_angle` is the magic-angle.
method = ThreeQ_VAS(
    channels=["87Rb"],
    magnetic_flux_density=7,  # in T
    spectral_dimensions=[
        {
            "count": 256,
            "spectral_width": 4e3,  # in Hz
            "reference_offset": -5e3,  # in Hz
            "label": "Isotropic dimension",
        },
        {
            "count": 512,
            "spectral_width": 1e4,  # in Hz
            "reference_offset": -4e3,  # in Hz
            "label": "MAS dimension",
        },
    ],
)

# %%
# Create the Simulator object, add the method and spin system objects, and
# run the simulation.
sim = Simulator()
sim.spin_systems = spin_systems  # add the spin systems
sim.methods = [method]  # add the method.
sim.run()

# %%
# The plot of the simulation.
data = sim.methods[0].simulation
ax = plt.gca(projection="csdm")
ax.imshow(data / data.max(), aspect="auto", cmap="gist_ncar_r")
ax.invert_xaxis()
ax.invert_yaxis()
plt.tight_layout()
plt.show()

# %%
# Add post-simulation signal processing.
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
# The plot of the simulation after signal processing.
ax = plt.subplot(projection="csdm")
ax.imshow(processed_data.real, cmap="gist_ncar_r", aspect="auto")
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
