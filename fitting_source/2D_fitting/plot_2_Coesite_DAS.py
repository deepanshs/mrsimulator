#!/usr/bin/env python
"""
¹⁷O 2D DAS NMR of Coesite
^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# Coesite is a high-pressure (2-3 GPa) and high-temperature (700°C) polymorph of silicon
# dioxide :math:`\text{SiO}_2`. Coesite has five crystallographic :math:`^{17}\text{O}`
# sites. The experimental dataset used in this example is published in
# Grandinetti `et al.` [#f1]_
import numpy as np
import csdmpy as cp
import matplotlib.pyplot as plt
from lmfit import Minimizer

from mrsimulator import Simulator
from mrsimulator import signal_processor as sp
from mrsimulator.utils import spectral_fitting as sf
from mrsimulator.utils import get_spectral_dimensions
from mrsimulator.utils.collection import single_site_system_generator
from mrsimulator.method import Method, SpectralDimension, SpectralEvent, MixingEvent

# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
filename = "https://ssnmr.org/sites/default/files/mrsimulator/DASCoesite.csdf"
experiment = cp.load(filename)

# standard deviation of noise from the dataset
sigma = 921.6698

# For spectral fitting, we only focus on the real part of the complex dataset
experiment = experiment.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in experiment.dimensions]

# plot of the dataset.
max_amp = experiment.max()
levels = (np.arange(14) + 1) * max_amp / 15  # contours are drawn at these levels.
options = dict(levels=levels, alpha=0.75, linewidths=0.5)  # plot options

plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", **options)
ax.invert_xaxis()
ax.set_ylim(30, -30)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# **Guess model**
#
# Create a guess list of spin systems.

shifts = [29, 39, 54.8, 51, 56]  # in ppm
Cq = [6.1e6, 5.4e6, 5.5e6, 5.5e6, 5.1e6]  # in  Hz
eta = [0.1, 0.2, 0.15, 0.15, 0.3]
abundance_ratio = [1, 1, 2, 2, 2]
abundance = np.asarray(abundance_ratio) / 8 * 100  # in %

spin_systems = single_site_system_generator(
    isotope="17O",
    isotropic_chemical_shift=shifts,
    quadrupolar={"Cq": Cq, "eta": eta},
    abundance=abundance,
)

# %%
# **Method**
#
# Create the DAS method.

# Get the spectral dimension parameters from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

DAS = Method(
    channels=["17O"],
    magnetic_flux_density=11.744,  # in T
    rotor_frequency=np.inf,
    spectral_dimensions=[
        SpectralDimension(
            **spectral_dims[0],
            events=[
                SpectralEvent(
                    fraction=0.5,
                    rotor_angle=37.38 * np.pi / 180,  # in rads
                    transition_queries=[{"ch1": {"P": [-1], "D": [0]}}],
                ),
                MixingEvent(query="NoMixing"),
                SpectralEvent(
                    fraction=0.5,
                    rotor_angle=79.19 * np.pi / 180,  # in rads
                    transition_queries=[{"ch1": {"P": [-1], "D": [0]}}],
                ),
                MixingEvent(query="NoMixing"),
            ],
        ),
        # The last spectral dimension block is the direct-dimension
        SpectralDimension(
            **spectral_dims[1],
            events=[
                SpectralEvent(
                    rotor_angle=54.735 * np.pi / 180,  # in rads
                    transition_queries=[{"ch1": {"P": [-1], "D": [0]}}],
                )
            ],
        ),
    ],
    experiment=experiment,  # also add the measurement to the method.
)

# %%
# **Guess Spectrum**

# Simulation
# ----------
sim = Simulator(spin_systems=spin_systems, methods=[DAS])
sim.config.number_of_sidebands = 1  # no sidebands are required for this dataset.
sim.run()

# Post Simulation Processing
# --------------------------
processor = sp.SignalProcessor(
    operations=[
        # Gaussian convolution along both dimensions.
        sp.IFFT(dim_index=(0, 1)),
        sp.apodization.Gaussian(FWHM="0.15 kHz", dim_index=0),
        sp.apodization.Gaussian(FWHM="0.1 kHz", dim_index=1),
        sp.FFT(dim_index=(0, 1)),
        sp.Scale(factor=4e7),
    ]
)
processed_dataset = processor.apply_operations(dataset=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", **options)
ax.contour(processed_dataset, colors="r", linestyles="--", **options)
ax.invert_xaxis()
ax.set_ylim(30, -30)
plt.grid()
plt.tight_layout()
plt.show()


# %%
# Least-squares minimization with LMFIT
# -------------------------------------
# Use the :func:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params` for a quick
# setup of the fitting parameters.
params = sf.make_LMFIT_params(sim, processor)
print(params.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

# %%
# **Solve the minimizer using LMFIT**
opt = sim.optimize()
minner = Minimizer(sf.LMFIT_min_function, params, fcn_args=(sim, processor, sigma, opt))
result = minner.minimize(method="powell")
result

# %%
# The best fit solution
# ---------------------
best_fit = sf.bestfit(sim, processor)[0].real

# Plot the spectrum
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", **options)
ax.contour(best_fit, colors="r", linestyles="--", **options)
ax.invert_xaxis()
ax.set_ylim(30, -30)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# The best fit solution
# ---------------------
residuals = sf.residuals(sim, processor)[0].real

fig, ax = plt.subplots(
    1, 3, sharey=True, figsize=(10, 3.0), subplot_kw={"projection": "csdm"}
)
vmax, vmin = experiment.max(), experiment.min()
for i, dat in enumerate([experiment, best_fit, residuals]):
    ax[i].imshow(dat, aspect="auto", vmax=vmax, vmin=vmin)
    ax[i].invert_xaxis()
ax[0].set_ylim(30, -30)
plt.tight_layout()
plt.show()
# %%
# .. [#f1] Grandinetti, P. J., Baltisberger, J. H., Farnan, I., Stebbins, J. F.,
#       Werner, U. and Pines, A.
#       Solid-State :math:`^{17}\text{O}` Magic-Angle and Dynamic-Angle Spinning NMR
#       Study of the :math:`\text{SiO}_2` Polymorph Coesite, J. Phys. Chem. 1995,
#       **99**, *32*, 12341-12348.
#       `DOI: 10.1021/j100032a045 <https://doi.org/10.1021/j100032a045>`_
