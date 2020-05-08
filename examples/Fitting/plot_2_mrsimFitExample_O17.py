#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fitting Sodium Silicate.
^^^^^^^^^^^^^^^^^^^^^^^^
"""
# sphinx_gallery_thumbnail_number = 3
#%%
# Often, after obtaining an NMR measurement we must fit tensors to our data so we can
# obtain the tensor parameters. In this example, we will illustrate the use of the *mrsimulator*
# method to simulate the experimental spectrum and fit the simulation to the data allowing us to
# extract the tensor parameters for our isotopomers. We will be using the `LMFIT <https://lmfit.github.io/lmfit-py/>`_
# methods to establish fitting parameters and fit the spectrum. The following examples will show fitting with
# two synthetic :math:`^{29}\text{Si}` spectra--cuspidine and wollastonite--as well as the
# measurements from an :math:`^{17}\text{O}` experiment on :math:`\text{Na}_{2}\text{SiO}_{3}`.
# The *mrsimulator* library and data make use of CSDM compliant files.
# In this example we will fit a simulation to an experimentally obtained :math:`^{17}\text{O}` spectrum.
# We use the :math:`^{17}\text{O}` tensor information from Grandinetti et. al. [#f5]_
#
# We will begin by importing *matplotlib* and establishing figure size.
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt

params = {"figure.figsize": (4.5, 3), "font.size": 9}
pylab.rcParams.update(params)

#%%
# Next we will import `csdmpy <https://csdmpy.readthedocs.io/en/latest/index.html>`_ and loading the data file.
import csdmpy as cp
import numpy as np


filename = "https://osu.box.com/shared/static/kfgt0jxgy93srsye9pofdnoha6qy58qf.csdf"
oxygen_experiment = cp.load(filename).real
oxygen_experiment.dimensions[0].to("ppm", "nmr_frequency_ratio")

x1, y1 = oxygen_experiment.to_list()

plt.plot(x1, y1)
plt.xlabel("$^{17}$O frequency / ppm")
plt.xlim(x1.value.max(), x1.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()


#%%
# Next, we will want to create a ``simulator`` object that we will use to fit to our
# spectrum. We will need to import the necessary libraries for the *mrsimulator*
# methods. We will then create ``isotopomer`` objects.

#%%

from mrsimulator import Simulator, Isotopomer, Site
from mrsimulator import SymmetricTensor as st
from mrsimulator.methods import BlochDecayFT

sim = Simulator()
O17_1 = Site(
    isotope="17O",
    isotropic_chemical_shift=60.0,
    quadrupolar={"Cq": 4200000, "eta": 0.5},
)
O17_2 = Site(
    isotope="17O", isotropic_chemical_shift=40, quadrupolar={"Cq": 2400000, "eta": 0.18}
)

isotopomers = [Isotopomer(sites=[site]) for site in [O17_1, O17_2]]


count = oxygen_experiment.dimensions[0].count
increment = oxygen_experiment.dimensions[0].increment.to("Hz").value
offset = oxygen_experiment.dimensions[0].coordinates_offset.to("Hz").value

method = BlochDecayFT(
    channel="17O",
    magnetic_flux_density=9.4,
    rotor_frequency=14000,
    dimensions=[
        {
            "count": count,
            "spectral_width": count * increment,
            "reference_offset": offset,
        }
    ],
)

sim.isotopomers += isotopomers
sim.methods += [method]

# To avoid querying at every iteration we will save the relevant transition pathways
for iso in isotopomers:
    iso.transition_pathways = method.get_transition_pathways(iso).tolist()
sim.run()

sim.methods[0].simulation.dimensions[0].to("ppm", "nmr_frequency_ratio")

x, y = sim.methods[0].simulation.to_list()

plt.plot(x, y)
plt.xlabel("$^{17}$O frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()

#%%
# Once we have our simulation we must create our list of parameters to use in our
# fitting. We will be using the `Parameters <https://lmfit.github.io/lmfit-py/parameters.html>`_ class from *LMFIT*.
#
# In each experiment the number of isotopomers and sites present as well as their
# attributes may vary. To simplify the parameter list creation we will use the
# :func:`~mrsimulator.spectral_fitting.make_fitting_parameters`

#%%


from mrsimulator.spectral_fitting import make_fitting_parameters

params = make_fitting_parameters(sim)
params.add(name="sigma", value=100, min=0)
params.add(name="factor", value=oxygen_experiment.max() / 4, min=0)

#%%
# With an experimental spectrum, a simulaton, and a list of parameters we are now
# ready to perform a fit. This fit will be performed using the *LMFIT* library as
# well as our error function, :func:`~mrsimulator.spectral_fitting.min_function`. The arguments
# for ``min_function`` are the intensities from the experimental data and the simulation
# CSDM object. Reporting the results of the fit will give us our tensor parameters.
#
# One thing to note is that the names of our parameters must correspond to their addresses within the simulation object
# in order to update the simulation during the fit. The *LMFIT* library does not allow for the use of special characters
# such as "\[", "\]", or "." so our current workaround is converting the special characters to their corresponding HTML
# character code numbers and converting back to the special character when updating the simulation.


#%%

from mrsimulator.spectral_fitting import min_function
from lmfit import Minimizer, report_fit

minner = Minimizer(
    min_function, params, fcn_args=(oxygen_experiment, sim, "Lorentzian")
)
result = minner.minimize()

#%%
# Next, we can compare the fit to the experimental data:

#%%


plt.figsize = (4, 3)
residual = oxygen_experiment.copy()
residual[:] = result.residual
plt.plot(*oxygen_experiment.to_list(), label="Spectrum")
plt.plot(*(oxygen_experiment - residual).to_list(), "r", alpha=0.5, label="Fit")
plt.plot(*residual.to_list(), alpha=0.5, label="Residual")

plt.xlim(x.value.max(), x.value.min())
plt.xlabel("$^{17}$O frequency / ppm")
plt.grid(which="major", axis="both", linestyle="--")
plt.legend()

plt.tight_layout()
plt.show()

#%%
# .. [#f5] T. M. Clark, P. Florian, J. F. Stebbins, and P. J. Grandinetti,
#       An :math:`^{17}\text{O}` NMR Investigation of Crystalline Sodium Metasilicate:
#       Implications for the Determination of Local Structure in Alkali Silicates,
#       J. Phys. Chem. B. 2001, **105**, 12257-12265.
#       `DOI: 10.1021/jp011289p  <https://doi.org/10.1021/jp011289p>`_
