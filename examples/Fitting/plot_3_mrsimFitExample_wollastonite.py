#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fitting Wollastonite.
^^^^^^^^^^^^^^^^^^^^^
"""
# sphinx_gallery_thumbnail_number = 4
#%%
# Often, after obtaining an NMR measurement we must fit tensors to our data so we can
# obtain the tensor parameters. In this example, we will illustrate the use of the *mrsimulator*
# method to simulate the experimental spectrum and fit the simulation to the data allowing us to
# extract the tensor parameters for our isotopomers. We will be using the `LMFIT <https://lmfit.github.io/lmfit-py/>`_
# methods to establish fitting parameters and fit the spectrum. The following examples will show fitting with
# two synthetic :math:`^{29}\text{Si}` spectra--cuspidine and wollastonite--as well as the
# measurements from an :math:`^{17}\text{O}` experiment on :math:`\text{Na}_{2}\text{SiO}_{3}`.
# The *mrsimulator* library and data make use of CSDM compliant files.
#
# In this example we will be creating a synthetic spectrum of wollastonite from reported tensor
# parameters and then fit a simulation to the spectrum to demonstrate a fitting procedure.
# The :math:`^{29}\text{Si}` tensor parameters were obtained from Hansen et. al. [#f1]_
#
# We will begin by importing *matplotlib* and establishing figure size.
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt

params = {"figure.figsize": (4.5, 3), "font.size": 9}
pylab.rcParams.update(params)

# Next, we will create a synthetic spectrum. In order to do so we must create a simulation.

from mrsimulator import Dimension, Isotopomer, Simulator, Site
from mrsimulator import SymmetricTensor as st
from mrsimulator.methods import one_d_spectrum

params = {"figure.figsize": (4.5, 3)}
pylab.rcParams.update(params)

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

dimension = Dimension(
    isotope="29Si",
    magnetic_flux_density=14.1,  # in T
    number_of_points=2046,
    spectral_width=25000,  # in Hz
    reference_offset=-10000,  # in Hz
    rotor_frequency=0,  # in Hz
)


isotopomers = [Isotopomer(sites=[site]) for site in [S29_1, S29_2, S29_3]]

synth_wollastonite = Simulator()
synth_wollastonite.isotopomers += isotopomers
synth_wollastonite.dimensions += [dimension]

x1, y1 = synth_wollastonite.run(method=one_d_spectrum)

plt.plot(x1, y1)
plt.xlabel("$^{29}$Si frequency / ppm")
plt.xlim(x1.value.max(), x1.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()

#%%
# In order to create a synthetic spectrum we will add noise and line broadening to the
# spectrum. The line broadening is performed using :meth:`mrsimulator.Apodization.apodize` which is a
# member of a simulator object.

#%%


from mrsimulator.apodization import Apodization
import numpy as np

synth_data = synth_wollastonite.as_csdm_object()
synth_data.dimensions[0].to("ppm", "nmr_frequency_ratio")

x2 = synth_data.dimensions[0].coordinates
y2 = synth_data.dependent_variables[0].components[0]
# Here we apodize the spectrum using a Lorentzian line broadening
synth_data.dependent_variables[0].components[0] = synth_wollastonite.apodize(
    Apodization.Lorentzian, sigma=200
)

# Here we simulate noise with RNG and add the noise to the spectrum
noise = np.random.normal(-0.025, 0.025, y2.shape)
synthetic = y2 + noise
synth_data.dependent_variables[0].components[0] = synthetic


plt.plot(x2, synthetic)
plt.xlabel("$^{29}$Si frequency / ppm")
plt.xlim(x2.value.max(), x2.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()

#%%
# We have just created a synthetic experimental data set and now we must create the simulation
# to fit to this synthetic data. We will use an initial approximation to generate our fitting simulation:

#%%


S29_1 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-91.0,
    shielding_symmetric={"zeta": 60, "eta": 0.6},
)
S29_2 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-90,
    shielding_symmetric={"zeta": 50, "eta": 0.70},
)
S29_3 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-88,
    shielding_symmetric={"zeta": 70, "eta": 0.55},
)

isotopomers = [Isotopomer(sites=[site]) for site in [S29_1, S29_2, S29_3]]

sim = Simulator()
sim.isotopomers += isotopomers
sim.dimensions += [dimension]
x3, y3 = sim.run(method=one_d_spectrum)

plt.plot(x3, y3)
plt.xlabel("$^{29}$Si frequency / ppm")
plt.xlim(x3.value.max(), x3.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()

#%%
# Again, we will need a list of parameters to change during the fitting. We create this list using the
# :func:`~mrsimulator.spectral_fitting.make_fitting_parameters` function.

#%%


from mrsimulator.spectral_fitting import make_fitting_parameters

params = make_fitting_parameters(sim)
params.add(name="sigma", value=0)
params.add(name="factor", value=1, min=0)

#%%
# With the synthetic data, simulation, and the parameters we are ready to perform the fit

#%%
from mrsimulator.spectral_fitting import min_function
from lmfit import Minimizer

minner = Minimizer(min_function, params, fcn_args=(synth_data, sim, "Lorentzian"))
result = minner.minimize()

#%%
# After the fit, we can plot the new simulated spectrum using the *matplotlib* library.

plt.figsize = (4, 3)
residual = synth_data.copy()
residual[:] = result.residual
plt.plot(*synth_data.to_list(), label="Spectrum")
plt.plot(*(synth_data - residual).to_list(), "r", alpha=0.5, label="Fit")
plt.plot(*residual.to_list(), alpha=0.5, label="Residual")

plt.xlim(x3.value.max(), x3.value.min())
plt.xlabel("$^{29}$Si frequency / ppm")
plt.grid(which="major", axis="both", linestyle="--")
plt.legend()

plt.tight_layout()
plt.show()

#%%
# .. [#f1] Hansen, M. R., Jakobsen, H. J., Skibsted, J., :math:`^{29}\text{Si}`
#       Chemical Shift Anisotropies in Calcium Silicates from High-Field
#       :math:`^{29}\text{Si}` MAS NMR Spectroscopy, Inorg. Chem. 2003,
#       **42**, *7*, 2368-2377.
#       `DOI: 10.1021/ic020647f <https://doi.org/10.1021/ic020647f>`_
