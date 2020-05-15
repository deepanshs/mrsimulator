#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Wollastonite
^^^^^^^^^^^^

29Si (I=1/2) spinning sideband simulation.
"""
#%%
# Wollastonite is a high-temperature calcium-silicate,
# :math:`\beta−\text{Ca}_3\text{Si}_3\text{O}_9`, with three distinct
# :math:`^{29}\text{Si}` sites. The :math:`^{29}\text{Si}` tensor parameters
# were obtained from Hansen et. al. [#f1]_
import matplotlib.pyplot as plt
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site

#%%
# **Step 1** Create the sites.

S29_1 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-89.0,  # in ppm
    shielding_symmetric={"zeta": 59.8, "eta": 0.62},  # zeta in ppm
)
S29_2 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-89.5,  # in ppm
    shielding_symmetric={"zeta": 52.1, "eta": 0.68},  # zeta in ppm
)
S29_3 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-87.8,  # in ppm
    shielding_symmetric={"zeta": 69.4, "eta": 0.60},  # zeta in ppm
)

sites = [S29_1, S29_2, S29_3]  # all sites

#%%
# **Step 2** Create the isotopomers from these sites. Again, we create three
# single-site isotopomers for better performance.

isotopomers = [Isotopomer(sites=[site]) for site in sites]

#%%
# **Step 3** Create a Bloch decay spectrum method.

from mrsimulator.methods import BlochDecaySpectrum

method = BlochDecaySpectrum(
    channels=["29Si"],
    magnetic_flux_density=14.1,  # in T
    rotor_frequency=1500,  # in Hz
    spectral_dimensions=[
        {
            "count": 2046,
            "spectral_width": 25000,  # in Hz
            "reference_offset": -10000,  # in Hz
        }
    ],
)


#%%
# **Step 4** Create the Simulator object and add the method and isotopomer objects.

sim_wollastonite = Simulator()
sim_wollastonite.isotopomers += isotopomers  # add isotopomers
sim_wollastonite.methods += [method]  # add method

#%%
# **Step 5** Simulate the spectrum.

sim_wollastonite.run()

#%%
# **Step 6** The plot of the simulation.

x, y = sim_wollastonite.methods[0].simulation.to_list()
plt.figure(figsize=(4, 3))
plt.plot(x, y, color="black", linewidth=1)
plt.xlabel("frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()

#%%
# .. [#f1] Hansen, M. R., Jakobsen, H. J., Skibsted, J., :math:`^{29}\text{Si}`
#       Chemical Shift Anisotropies in Calcium Silicates from High-Field
#       :math:`^{29}\text{Si}` MAS NMR Spectroscopy, Inorg. Chem. 2003,
#       **42**, *7*, 2368-2377.
#       `DOI: 10.1021/ic020647f <https://doi.org/10.1021/ic020647f>`_
