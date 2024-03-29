#!/usr/bin/env python
"""
Extended Czjzek fitting of ¹³⁹La MAS NMR of La₀.₂Y₁.₈Si₂2O₇
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# The following is a demonstration on how to fit tensor distributions to an experimental
# spectrum using a structure-forward approach. The dataset was acquired from a
# sample of :math:`\text{La}_{0.2}\text{Y}_{1.8}\text{Si}_2\text{O}_7` silicate glass
# and shared by Fernańdez-Carrioń `et al.` [#f1]_.
from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.method.lib import BlochDecayCTSpectrum
from mrsimulator.models import ExtCzjzekDistribution
from mrsimulator.simulator import Isotope
from mrsimulator.utils import get_spectral_dimensions
from mrsimulator import signal_processor as sp

import numpy as np
import csdmpy as cp
import lmfit
from lmfit import Minimizer
import matplotlib.pyplot as plt


# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
host = "http://ssnmr.org/sites/default/files/"
filename = "mrsimulator/La bc - LaY90 cpmg.csdf"
experiment = cp.load(host + filename).real
experiment.x[0].to("ppm", "nmr_frequency_ratio")
experiment /= experiment.max()

plt.figure(figsize=(4, 3))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", alpha=0.5)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# Calculate Noise Standard Deviation
# ----------------------------------
loc = np.where(experiment.dimensions[0].coordinates < -5000e-6)
sigma = experiment[loc].std()
print(sigma)

# %%
# Setup Method and Processor
# --------------------------
#
# Next setup a method which simulates the experiment, again using the
# :py:meth:`~mrsimulator.utils.get_spectral_dimensions` function to get the spectral
# dimension parameters from the experimental dataset.
# We also create a signal processor object which will scale the simulated spectrum.
spec_dims = get_spectral_dimensions(experiment)
method = BlochDecayCTSpectrum(
    channels=["139La"], magnetic_flux_density=17.6, spectral_dimensions=spec_dims
)

processor = sp.SignalProcessor(operations=[sp.Scale(factor=1000)])


# %%
# Create a Lineshape Kernel
# -------------------------
#
# An approach for fitting model values to the experimental dataset could be computing
# the spectrum from a list of spin system objects with tensor parameters pulled from
# an Extended Czjzek model during each round of fitting, but this would be inefficient!
# The grid of :math:`\text{Cq}`
# and :math:`\eta_q` points will not change during the fit; we are interested in the
# probability distribution across these points.
#
# We therefore define a function which pre-computes a lineshape kernel given a range
# of :math:`\text{Cq}` and :math:`\eta_q` values, an NMR method, and an isotope which,
# for this example, is ``"139La"``.
def make_kernel(pos, method, isotope):
    """Pre-computes the kernel to use when fitting the experimental spectrum

    Arguments:
        (np.array) pos: list of Numpy array of cq and eta points
        (Method) method: mrsimulator Method object used to simulate the spectra
        (str) isotope: Isotope of the site

    Returns:
        Lineshape kernel as a numpy array
    """
    Cq, eta = np.meshgrid(pos[0], pos[1], indexing="xy")

    spin_systems = [
        SpinSystem(sites=[Site(isotope=isotope, quadrupolar=dict(Cq=cq * 1e6, eta=e))])
        for cq, e in zip(Cq.ravel(), eta.ravel())
    ]
    sim = Simulator(spin_systems=spin_systems, methods=[method])
    sim.config.number_of_sidebands = 4
    sim.config.decompose_spectrum = "spin_system"
    sim.run(pack_as_csdm=False)  # Will return spectrum as numpy array, not CSDM object

    amp = sim.methods[0].simulation.real
    return amp.T


# Create ranges to construct cq and eta grid points
cq_range = np.linspace(0, 100, num=100) * 0.8 + 25
eta_range = np.arange(21) / 20
pos = [cq_range, eta_range]
kernel = make_kernel(pos, method, "139La")


# %%
# Make spectrum from kernel and model values
# ------------------------------------------
#
# Next define a function which takes a LMFIT Parameters object, the pre-computed kernel
# and a signal processor object and returns a CSDM object holding the guess spectrum.
# Here the Parameters object holds attributes from the
# :py:class:`~mrsimulator.models.ExtCzjzekDistribution` class which is used to
# create the probability density function, :math:`{\bf f}`. Then :math:`{\bf f}` is
# multiplied with the kernel, :math:`{\bf K}`, to produce a spectrum, :math:`{\bf s}`.
#
# .. math::
#
#     {\bf s} = {\bf K \cdot f},
#
# A constant isotropic chemical shift, whose value is also in the Parameters object,
# is applied to :math:`{\bf s}` using the `Fourier shift relation
# <https://en.wikipedia.org/wiki/Fourier_transform#Translation_/_time_shifting>`__.
# Finally, the spectrum is scaled and returned.
def make_spectrum_from_parameters(params, kernel, processor, pos, distribution):
    """Makes a spectrum with values given in a parameters object and a pre-computed
    kernel.

    Arguments:
        (Parameters) params: LMFIT Parameters object with values
        (np.array) kernel: Pre-computed spectrum kernel
        (sp.SignalProcessor) processor: SignalProcessor to apply to spectrum

    Returns:
        CSDM object of spectrum
    """
    # Extract values from parameters object
    values = params.valuesdict()
    Cq = values["dist_Cq"]
    eta = values["dist_eta"]
    eps = values["dist_eps"]
    iso_shift = values["dist_iso_shift"]

    # Setup model object and get the amplitude
    distribution.symmetric_tensor.Cq = Cq
    distribution.symmetric_tensor.eta = eta
    distribution.eps = eps
    _, _, amp = distribution.pdf(pos=pos)

    # Create spectra by dotting the amplitude distribution with the kernel
    dist = np.dot(kernel, amp.ravel())

    # Pack numpy array as csdm object and apply signal processing
    guess_dataset = cp.CSDM(
        dimensions=experiment.x,
        dependent_variables=[cp.as_dependent_variable(dist)],
    )

    # Calculate isotropic shift in Hz
    larmor_freq = Isotope(symbol="139La").gyromagnetic_ratio * 17.6
    iso_shift_in_hz = larmor_freq * iso_shift

    # Apply isotropic shift using FFT shift theorem
    guess_dataset = guess_dataset.fft()
    time_coords = guess_dataset.x[0].coordinates.value
    guess_dataset.y[0].components[0] *= np.exp(
        -np.pi * 2j * iso_shift_in_hz * time_coords
    )
    guess_dataset = guess_dataset.fft()

    # Apply signal processor and return
    processor.operations[0].factor = values["sp_scale_factor"]
    guess_dataset = processor.apply_operations(guess_dataset)
    return guess_dataset.real


def residuals(exp_spectra, simulated_spectra):
    """Returns the difference between exp_spectra and simulated_spectra"""
    return exp_spectra - simulated_spectra


# %%
# Make initial guess
# ------------------
#
# Next pack the attributes to be fit into a Parameters object and plot the initial guess
# spectrum.
#
# .. note::
#
#     If you adapt this example to your own dataset, make sure the initial guess is
#     decently good, otherwise LMFIT is likely to fall into a local minima.
params = lmfit.Parameters()
params.add("dist_Cq", value=49.5)
params.add("dist_eta", value=0.55, min=0, max=1)
params.add("dist_eps", value=0.24, min=0)
params.add("dist_iso_shift", value=350)
params.add("sp_scale_factor", value=3.8e3, min=0)

# Plot the initial guess spectrum along with the experimental data
distribution = ExtCzjzekDistribution(
    symmetric_tensor={"Cq": params["dist_Cq"], "eta": params["dist_eta"]},
    eps=params["dist_eps"],
)
initial_guess = make_spectrum_from_parameters(
    params, kernel, processor, pos, distribution
)
residual_spectrum = residuals(experiment, initial_guess)

plt.figure(figsize=(4, 3))
ax = plt.subplot(projection="csdm")
ax.plot(experiment.real, "k", alpha=0.5, label="Experiment")
ax.plot(initial_guess.real, "r", alpha=0.3, label="Guess")
ax.plot(residual_spectrum.real, "b", alpha=0.3, label="Residuals")
plt.legend()
plt.grid()
plt.title("Initial Guess")
plt.tight_layout()
plt.show()


# %%
# Least-squares minimization with LMFIT
# -------------------------------------
#
# Now define minimization function per the LMFIT specifications which will return
# the difference between the experimental and guess spectrum scaled by the noise
# standard deviation which was calculated in a previous cell.
def minimization_function(
    params, experiment, processor, kernel, pos, distribution, sigma=sigma
):
    guess_spectrum = make_spectrum_from_parameters(
        params, kernel, processor, pos, distribution
    )
    residual_spectrum = residuals(experiment, guess_spectrum)
    return residual_spectrum.y[0].components[0].real / sigma


# %%
scipy_minimization_kwargs = dict(
    diff_step=1e-3,  # Increase step size
    gtol=1e-15,  # Increase global convergence requirement (default 1e-8)
    xtol=1e-15,  # Increase variable convergence requirement (default 1e-8)
    verbose=2,  # Print minimization info during each step
    loss="soft_l1",
)

minner = Minimizer(
    minimization_function,
    params,
    fcn_args=(experiment, processor, kernel, pos, distribution),
    **scipy_minimization_kwargs,
)
result = minner.minimize(method="least_squares")
best_fit_params = result.params  # Grab the Parameters object from the best fit
result


# %%
# Plot the best-fit solution
# --------------------------
#
# Finally, plot the best fit spectrum and the residuals.
final_fit = make_spectrum_from_parameters(
    best_fit_params, kernel, processor, pos, distribution
)
residual_spectrum = residuals(experiment, final_fit)

plt.figure(figsize=(4, 3))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", alpha=0.5, label="Experiment")
ax.plot(final_fit, "r", alpha=0.3, label="Fit")
ax.plot(residual_spectrum, "b", alpha=0.3, label="Residuals")
plt.legend()
ax.set_xlim(-11000, 9000)
plt.grid()
plt.title("Best Fit")
plt.tight_layout()
plt.show()

# %%
#
# .. [#f1] A. J. Fernández-Carrión, M. Allix, P. Florian, M. R. Suchomel, and A. I.
#       Becerro.
#       Revealing Structural Detail in the High Temperature La2Si2O7–Y2Si2O7 Phase
#       Diagram by Synchrotron Powder Diffraction and Nuclear Magnetic Resonance
#       Spectroscopy.
#       The Journal of Physical Chemistry C 2012 **116** (40), 21523-21535
#       `DOI: 10.1021/jp305777m <https://doi.org/10.1021/jp305777m>`_
