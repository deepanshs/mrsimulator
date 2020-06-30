#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fitting Sodium Silicate.
^^^^^^^^^^^^^^^^^^^^^^^^
.. sectionauthor:: Maxwell C. Venetos <maxvenetos@gmail.com>
"""
# sphinx_gallery_thumbnail_number = 3
#%%
# Often, after obtaining an NMR measurement we must fit tensors to our data so we can
# obtain the tensor parameters. In this example, we will illustrate the use of the *mrsimulator*
# method to simulate the experimental spectrum and fit the simulation to the data allowing us to
# extract the tensor parameters for our spin systems. We will be using the `LMFIT <https://lmfit.github.io/lmfit-py/>`_
# methods to establish fitting parameters and fit the spectrum. The following examples will show fitting with
# two synthetic :math:`^{29}\text{Si}` spectra--cuspidine and wollastonite--as well as the
# measurements from an :math:`^{17}\text{O}` experiment on :math:`\text{Na}_{2}\text{SiO}_{3}`.
# The *mrsimulator* library and data make use of CSDM compliant files.
# In this example we will fit a simulation to an experimentally obtained :math:`^{17}\text{O}` spectrum.
# We use the :math:`^{17}\text{O}` tensor information from Grandinetti `et. al.` [#f5]_
#
# We will begin by importing *matplotlib* and establishing figure size.
import matplotlib as mpl
import matplotlib.pyplot as plt

font = {"weight": "light", "size": 9}
mpl.rc("font", **font)
mpl.rcParams["figure.figsize"] = [4.25, 3.0]

#%%
# Next we will import `csdmpy <https://csdmpy.readthedocs.io/en/latest/index.html>`_ and loading the data file.
import csdmpy as cp


filename = "https://osu.box.com/shared/static/kfgt0jxgy93srsye9pofdnoha6qy58qf.csdf"
oxygen_experiment = cp.load(filename).real
oxygen_experiment.dimensions[0].to("ppm", "nmr_frequency_ratio")

ax = plt.subplot(projection="csdm")
ax.plot(oxygen_experiment, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()


#%%
# Next, we will want to create a ``simulator`` object that we will use to fit to our
# spectrum. We will need to import the necessary libraries for the *mrsimulator*
# methods. We will then create ``SpinSystem`` objects.


from mrsimulator import SpinSystem
from mrsimulator import Simulator

# %%
# **Step 1** Create the site.

O17_1 = {
    "isotope": "17O",
    "isotropic_chemical_shift": "60.0 ppm",
    "quadrupolar": {"Cq": "4200000 Hz", "eta": 0.5},
}

O17_2 = {
    "isotope": "17O",
    "isotropic_chemical_shift": "40.0 ppm",
    "quadrupolar": {"Cq": "2400000 Hz", "eta": 0.5},
}

# %%
# **Step 2** Create the spin system for the site.

spin_system = {"sites": [O17_1, O17_2], "abundance": "100%"}  # from the above code

system_object = SpinSystem.parse_dict_with_units(spin_system)

# %%
# **Step 3** Create the Bloch Decay method.

from mrsimulator.methods import BlochDecayCentralTransitionSpectrum

count = oxygen_experiment.dimensions[0].count
increment = oxygen_experiment.dimensions[0].increment.to("Hz").value
offset = oxygen_experiment.dimensions[0].coordinates_offset.to("Hz").value

method = BlochDecayCentralTransitionSpectrum(
    channels=["17O"],
    magnetic_flux_density=9.4,  # in T
    rotor_frequency=14000,  # in Hz
    spectral_dimensions=[
        {
            "count": count,
            "spectral_width": count * increment,  # in Hz
            "reference_offset": offset,  # in Hz
        }
    ],
)

# %%
# **Step 4** Create the Simulator object and add the method and spin-system objects.

sim = Simulator()
sim.spin_systems += [system_object]
sim.methods += [method]

sim.methods[0].experiment = oxygen_experiment

# %%
# **Step 5** simulate the spectrum.

for iso in sim.spin_systems:
    # To avoid querying at every iteration we will save the relevant transition pathways
    iso.transition_pathways = method.get_transition_pathways(iso).tolist()
sim.run()

# %%
# **Step 6** Create a SignalProcessor

import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo

factor = oxygen_experiment.dependent_variables[0].components[0].max().real / 4

op_list = [
    sp.IFFT(dim_indx=0),
    apo.Exponential(Lambda=100, dim_indx=0, dep_var_indx=0),
    sp.FFT(dim_indx=0),
    sp.Scale(factor=factor),
]

post_sim = sp.SignalProcessor(data=sim.methods[0].simulation, operations=op_list)

# %%
# ** Step 7** Process and plot the spectrum.

processed_data = post_sim.apply_operations()

ax = plt.subplot(projection="csdm")
ax.plot(processed_data, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

#%%
# Once we have our simulation we must create our list of parameters to use in our
# fitting. We will be using the `Parameters <https://lmfit.github.io/lmfit-py/parameters.html>`_ class from *LMFIT*.
#
# In each experiment the number of spin systems and sites present as well as their
# attributes may vary. To simplify the parameter list creation we will use the
# :func:`~mrsimulator.spectral_fitting.make_fitting_parameters`

# %%
# **Step 8** Create a list of parameters to vary during fitting.


from mrsimulator.spectral_fitting import make_LMFIT_parameters

params = make_LMFIT_parameters(sim, post_sim)
params.pretty_print()
# %%
# **Step 9** Perform minimization.

from mrsimulator.spectral_fitting import LMFIT_min_function
from lmfit import Minimizer, report_fit

minner = Minimizer(LMFIT_min_function, params, fcn_args=(sim, post_sim))
result = minner.minimize()
report_fit(result)

# %%
# **Step 10** Plot the fitted spectrum.


plt.figsize = (4, 3)
residual = oxygen_experiment.copy()
residual[:] = result.residual
plt.plot(*oxygen_experiment.to_list(), label="Spectrum")
plt.plot(*(oxygen_experiment - residual).to_list(), "r", alpha=0.5, label="Fit")
plt.plot(*residual.to_list(), alpha=0.5, label="Residual")

plt.xlabel("$^{17}$O frequency / ppm")
plt.gca().invert_xaxis()
plt.grid(which="major", axis="both", linestyle="--")
plt.legend()

plt.tight_layout()
plt.show()

#%%
#
# .. [#f5] T. M. Clark, P. Florian, J. F. Stebbins, and P. J. Grandinetti,
#       An :math:`^{17}\text{O}` NMR Investigation of Crystalline Sodium Metasilicate:
#       Implications for the Determination of Local Structure in Alkali Silicates,
#       J. Phys. Chem. B. 2001, **105**, 12257-12265.
#       `DOI: 10.1021/jp011289p  <https://doi.org/10.1021/jp011289p>`_
