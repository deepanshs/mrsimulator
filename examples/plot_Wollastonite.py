#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Wollastonite
^^^^^^^^^^^^
"""
#%%
# Wollastonite is a high-temperature calcium-silicate,
# :math:`\beta−\text{Ca}_3\text{Si}_3\text{O}_9`,
# with three distinct :math:`^{29}\text{Si}` sites. The :math:`^{29}\text{Si}`
# tensor parameters
# were obtained from Hansen et. al. [#f1]_
import matplotlib.pyplot as plt
from mrsimulator import Dimension
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator.methods import one_d_spectrum

#%%
# **Step 1** Create sites.

S29_1 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-89.0,
    shielding_symmetric={"zeta": 59.8, "eta": 0.62},
)
S29_2 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-89.5,
    shielding_symmetric={"zeta": 52.1, "eta": 0.68},
)
S29_3 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-87.8,
    shielding_symmetric={"zeta": 69.4, "eta": 0.60},
)

sites = [S29_1, S29_2, S29_3]

#%%
# **Step 2** Create isotopomers from these sites.

isotopomers = [Isotopomer(sites=[site]) for site in sites]

#%%
# **Step 3** Create a dimension.

dimension = Dimension(
    isotope="29Si",
    magnetic_flux_density=14.1,  # in T
    number_of_points=2046,
    spectral_width=25000,  # in Hz
    reference_offset=-10000,  # in Hz
    rotor_frequency=1500,  # in Hz
)

#%%
# **Step 4** Create the Simulator object and add dimension and isotopomer objects.

sim_wollastonite = Simulator()

# add isotopomers
sim_wollastonite.isotopomers += isotopomers

# add dimensions
sim_wollastonite.dimensions += [dimension]

#%%
# **Step 5** Simulate the spectrum.

x, y = sim_wollastonite.run(method=one_d_spectrum)

#%%
# **Step 6** Plot.

_, ax = plt.subplots(1, 1, figsize=(4.5, 3))
ax.plot(x, y, color="black", linewidth=0.5)
ax.set_xlabel("frequency / ppm")
ax.set_xlim(x.value.max(), x.value.min())
ax.grid(color="gray", linestyle="--", linewidth=0.25)
plt.tight_layout()
plt.show()

#%%
# .. [#f1] Hansen, M. R., Jakobsen, H. J., Skibsted, J., :math:`^{29}\text{Si}`
#       Chemical Shift Anisotropies in Calcium Silicates from High-Field
#       :math:`^{29}\text{Si}` MAS NMR Spectroscopy, Inorg. Chem. 2003,
#       **42**, *7*, 2368-2377.
#       `DOI: 10.1021/ic020647f <https://doi.org/10.1021/ic020647f>`_
