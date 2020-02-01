#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Coesite
^^^^^^^

17O (I=5/2) quadrupolar line-shape simulation.
"""
#%%
# Coesite is a high-pressure (2-3 GPa) and high-temperature (700°C) polymorph of silicon
# dioxide :math:`\text{SiO}_2`. Coesite has five distinct :math:`^{17}\text{O}` sites.
# We use the :math:`^{17}\text{O}` tensor information from Grandinetti et. al. [#f2]_
import matplotlib.pyplot as plt
from mrsimulator import Dimension
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator.methods import one_d_spectrum

#%%
# **Step 1** Create sites.

O17_1 = Site(
    isotope="17O", isotropic_chemical_shift=29, quadrupolar={"Cq": 6.05e6, "eta": 0.000}
)
O17_2 = Site(
    isotope="17O", isotropic_chemical_shift=41, quadrupolar={"Cq": 5.43e6, "eta": 0.166}
)
O17_3 = Site(
    isotope="17O", isotropic_chemical_shift=57, quadrupolar={"Cq": 5.45e6, "eta": 0.168}
)
O17_4 = Site(
    isotope="17O", isotropic_chemical_shift=53, quadrupolar={"Cq": 5.52e6, "eta": 0.169}
)
O17_5 = Site(
    isotope="17O", isotropic_chemical_shift=58, quadrupolar={"Cq": 5.16e6, "eta": 0.292}
)

sites = [O17_1, O17_2, O17_3, O17_4, O17_5]

#%%
# **Step 2** Create isotopomers from these sites.

abundance = [0.83, 1.05, 2.16, 2.05, 1.90]  # abundance of each isotopomer
isotopomers = [Isotopomer(sites=[s], abundance=a) for s, a in zip(sites, abundance)]

#%%
# **Step 3** Create a dimension.

dimension = Dimension(
    isotope="17O", number_of_points=2046, spectral_width=50000, rotor_frequency=14000
)

#%%
# The above dimension is set up to record the :math:`^{17}\text{O}` resonances at the
# magic angle, spinning at 14 kHz and 9.4 T external magnetic flux density. The
# resonances are recorded over 50 kHz using 2046 points. You may also request a full
# description of the dimension object using the
# :meth:`mrsimulator.Dimension.to_dict_with_units` method.
from pprint import pprint

pprint(dimension.to_dict_with_units())

#%%
# **Step 4** Create the Simulator object and add dimension and isotopomer objects.

sim_coesite = Simulator()

# add isotopomers
sim_coesite.isotopomers += isotopomers

# add dimensions
sim_coesite.dimensions += [dimension]

#%%
# **Step 5** Simulate the spectrum.


x, y = sim_coesite.run(method=one_d_spectrum)

#%%
# **Step 6** Plot.

plt.figure(figsize=(4, 3))
plt.plot(x, y, color="black", linewidth=1)
plt.xlabel("frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()

#%%
# .. [#f2] Grandinetti, P. J., Baltisberger, J. H., Farnan, I., Stebbins, J. F.,
#       Werner, U. and Pines, A.
#       Solid-State :math:`^{17}\text{O}` Magic-Angle and Dynamic-Angle Spinning NMR
#       Study of the :math:`\text{SiO}_2` Polymorph Coesite, J. Phys. Chem. 1995,
#       **99**, *32*, 12341-12348.
#       `DOI: 10.1021/j100032a045 <https://doi.org/10.1021/j100032a045>`_
