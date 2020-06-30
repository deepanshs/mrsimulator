#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fitting Cusipidine.
^^^^^^^^^^^^^^^^^^^
.. sectionauthor:: Maxwell C. Venetos <maxvenetos@gmail.com>
"""
# sphinx_gallery_thumbnail_number = 3
# %%
# Often, after obtaining an NMR measurement we must fit tensors to our data so we can
# obtain the tensor parameters. In this example, we will illustrate the use of the *mrsimulator*
# method to simulate the experimental spectrum and fit the simulation to the data allowing us to
# extract the tensor parameters for our spin systems. We will be using the `LMFIT <https://lmfit.github.io/lmfit-py/>`_
# methods to establish fitting parameters and fit the spectrum. The following examples will show fitting with
# two synthetic :math:`^{29}\text{Si}` spectra--cuspidine and wollastonite--as well as the
# measurements from an :math:`^{17}\text{O}` experiment on :math:`\text{Na}_{2}\text{SiO}_{3}`.
# The *mrsimulator* library and data make use of CSDM compliant files.
# In this example we will be creating a synthetic spectrum of cuspidine from reported tensor
# parameters and then fit a simulation to the spectrum to demonstrate a simple fitting procedure.
# The :math:`^{29}\text{Si}` tensor parameters were obtained from Hansen `et. al.` [#f1]_
#
# We will begin by importing *matplotlib* and establishing figure size.
import matplotlib as mpl
import matplotlib.pyplot as plt

font = {"weight": "light", "size": 9}
mpl.rc("font", **font)
mpl.rcParams["figure.figsize"] = [4.25, 3.0]

# %%
# Next we will import `csdmpy <https://csdmpy.readthedocs.io/en/latest/index.html>`_ and loading the data file.
import csdmpy as cp

filename = "https://osu.box.com/shared/static/a45xj96iekdjrs2beri0nkrow4vjewdh.csdf"
synthetic_experiment = cp.load(filename)
synthetic_experiment.dimensions[0].to("ppm", "nmr_frequency_ratio")

ax = plt.subplot(projection="csdm")
ax.plot(synthetic_experiment, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# In order to to fit a simulation to the data we will need to establish a ``Simulator`` object. We will
# use approximate initial parameters to generate our simulation:

from mrsimulator import SpinSystem
from mrsimulator import Simulator

# %%
# **Step 1** Create the site.

site = {
    "isotope": "29Si",
    "isotropic_chemical_shift": "-80.0 ppm",
    "shielding_symmetric": {"zeta": "-60 ppm", "eta": 0.6},
}

# %%
# **Step 2** Create the spin system for the site.

spin_system = {
    "name": "Si Site",
    "description": "A 29Si site in cuspidine",
    "sites": [site],  # from the above code
    "abundance": "100%",
}
system_object = SpinSystem.parse_dict_with_units(spin_system)

# %%
# **Step 3** Create the Bloch Decay method.

from mrsimulator.methods import BlochDecaySpectrum

method = BlochDecaySpectrum(
    channels=["29Si"],
    magnetic_flux_density=7.1,  # in T
    rotor_frequency=780,  # in Hz
    spectral_dimensions=[
        {
            "count": 2048,
            "spectral_width": 25000,  # in Hz
            "reference_offset": -5000,  # in Hz
        }
    ],
)

# %%
# The above method is set up to record the :math:`^{29}\text{Si}` resonances at the
# magic angle, spinning at 780 Hz and 7.1 T external magnetic flux density.
# The resonances are recorded over 25 kHz spectral width ofset by -5000 Hz and
# using 2046 points.

# %%
# **Step 4** Create the Simulator object and add the method and spin-system objects.

sim = Simulator()
sim.spin_systems += [system_object]
sim.methods += [method]

sim.methods[0].experiment = synthetic_experiment

# %%
# **Step 5** simulate the spectrum.

sim.run()

# %%
# **Step 6** Create a SignalProcessor

import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo

op_list = [
    sp.IFFT(dim_indx=0),
    apo.Exponential(Lambda=200, dim_indx=0, dep_var_indx=0),
    sp.FFT(dim_indx=0),
    sp.Scale(factor=1),
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
# Next, we will need a list of parameters that will be used in the fit. the *LMFIT* library allows us to create
# a list of parameters rather easily using the `Parameters <https://lmfit.github.io/lmfit-py/parameters.html>`_ class.
# We have created a function to parse the
# ``simulator`` object for available parameters and construct an *LMFIT* ``Parameter`` object which is shown in the next two
# examples on fitting. Here, however, we will construct the parameter list explicitly to demonstrate how the parameters
# are created.

#%%


from lmfit import Minimizer, Parameters, fit_report

params = Parameters()

params.add(name="iso", value=-80)
params.add(name="eta", value=0.6, min=0, max=1)
params.add(name="zeta", value=-60)
params.add(name="Lambda", value=200)
params.add(name="factor", value=1)
params

#%%
# We will next set up an error function that will update the simulation throughout the minimization.
# We will construct a simple function here to demonstrate the *LMFIT* library, however, the next examples
# will showcase a fitting function provided in the *mrsimulator* library which automates the process.


def test_function(params, sim, post_sim):
    values = params.valuesdict()

    intensity = sim.methods[0].experiment.dependent_variables[0].components[0].real

    # Here, we update simulation parameters iso, eta, and zeta for the site object
    site = sim.spin_systems[0].sites[0]
    site.isotropic_chemical_shift = values["iso"]
    site.shielding_symmetric.eta = values["eta"]
    site.shielding_symmetric.zeta = values["zeta"]

    # here we run the simulation
    sim.run()

    post_sim.operations[3].factor = values["factor"]
    post_sim.operations[1].Lambda = values["Lambda"]
    # here we apodize the signal to simulate line broadening
    processed_data = post_sim.apply_operations()

    return intensity - processed_data.dependent_variables[0].components[0].real


#%%
# With the synthetic data, simulation, and the parameters we are ready to perform the fit. To fit, we use
# the *LMFIT* `Minimizer <https://lmfit.github.io/lmfit-py/fitting.html>`_ class.
# One consideration for the case of magic angle spinning fitting is we must use a discrete minimization method
# such as 'powell' as the chemical shift varies discretely

# %% **Step 8** Perform minimization.

minner = Minimizer(test_function, params, fcn_args=(sim, post_sim))
result = minner.minimize(method="powell")
print(fit_report(result))

#%%
# **Step 9** Plot the fitted spectrum.


plt.figsize = (4, 3)
residual = synthetic_experiment.copy()
residual[:] = result.residual
plt.plot(*synthetic_experiment.to_list(), label="Spectrum")
plt.plot(*(synthetic_experiment - residual).to_list(), "r", alpha=0.5, label="Fit")
plt.plot(*residual.to_list(), alpha=0.5, label="Residual")

plt.xlabel("Frequency / Hz")
plt.gca().invert_xaxis()
plt.grid(which="major", axis="both", linestyle="--")
plt.legend()

plt.tight_layout()
plt.show()

#%%
#
# .. [#f1] Hansen, M. R., Jakobsen, H. J., Skibsted, J., :math:`^{29}\text{Si}`
#       Chemical Shift Anisotropies in Calcium Silicates from High-Field
#       :math:`^{29}\text{Si}` MAS NMR Spectroscopy, Inorg. Chem. 2003,
#       **42**, *7*, 2368-2377.
#       `DOI: 10.1021/ic020647f <https://doi.org/10.1021/ic020647f>`_


# %%
