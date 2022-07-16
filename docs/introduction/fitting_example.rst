.. _fitting_example:

Least-Squares Fitting Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Mrsimulator** can interact with various Python data science
packages.  One such package,
`LMFIT <https://lmfit.github.io/lmfit-py/>`_, can be used to perform non-linear
least-squares analysis of experimental NMR spectra.

Here, we illustrate the use of the **mrsimulator** objects to

- import and prepare an experimental dataset for the least-squares analysis,
- create a fitting model using Simulator and SignalProcessor objects,
- use the fitting model to perform a least-squares analysis,
- extract the model parameters with uncertainties, and
- plot the experimental spectrum along with the best-fit simulation and residuals.

Import Experimental Dataset
---------------------------

In this example, you will apply the least-squares fitting procedure to a
:math:`^{27}\text{Al}` magic-angle spinning spectrum of :math:`\text{Al(acac)$_3$}`
measured with whole echo acquisition.

You will begin by importing an experimental dataset measured on a 9.4 T Bruker
AVANCE III HD NMR spectrometer into the script. Bruker datasets are saved in
folders with a number as the folder name. In this case, that folder has been
transferred from the spectrometer and renamed "Al_acac".



Download experimental dataset
'''''''''''''''''''''''''''''

For our purposes, the folder was also compressed into a zip archive and uploaded
to an internet-accessible server. You can use the code block below to download
the zip archive from the server and unzip it into the originally named folder.

.. plot::
    :context: reset

    import requests
    import zipfile
    from io import BytesIO

    file_ = "https://ssnmr.org/sites/default/files/mrsimulator/Al_acac3_0.zip"
    request = requests.get(file_)
    z = zipfile.ZipFile(BytesIO(request.content))
    z.extractall("Al_acac")


Convert experimental dataset to CSDM
''''''''''''''''''''''''''''''''''''

Now that the Bruker dataset folder is accessible to your Python code, you can
use the Python package `nmrglue <https://github.com/jjhelmus/nmrglue>`_ to
convert the dataset into a `CSDM <https://csdmpy.readthedocs.io/en/stable/>`_
object.


.. plot::
    :context: close-figs

    import nmrglue as ng

    # initialize nmrglue converter object
    converter = ng.convert.converter()

    # read in the bruker dataset file
    dic, data = ng.bruker.read("Al_acac")
    converter.from_bruker(dic, data)

    # convert to CSDM format
    csdm_ds = converter.to_csdm()

With the dataset converted into a CSDM object, plot the dataset to make sure
that you imported it correctly.

.. skip: next

.. plot::
    :context: close-figs

    import matplotlib.pyplot as plt

    plt.figure(figsize=(5, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(csdm_ds.real, label="real")
    ax.plot(csdm_ds.imag, label="imag")
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.show()

This is the raw time-domain dataset, acquired using whole-echo acquisition. The
blue and orange lines are the real and imaginary parts of the complex
time-domain signal. If you're working with your own experimental dataset and
have already processed it into the frequency domain, then you can skip the next
few steps and proceed to the `Measure Noise`_ section.

Process Experimental Dataset
----------------------------

Proceeding from here, you'll need to transform this dataset into the frequency
domain for the least-squares analysis. Before applying the Fourier transform,
however, two things need to be adjusted.

First, you need to adjust the ``coordinates_offset`` to place the time origin at
the top of the echo. You can find this time offset among the pulse sequence
parameters. If you acquired the signal with a simple Hahn-echo sequence,
i.e., :math:`\pi/2-\tau-\pi-t`, then the ``coordinates_offset`` should be the
time between the centers of the two pulses. However, there are often some
additional receiver delays before the signal acquisition begins, and those
times need to be subtracted from the interpulse spacing. In this measurement,
we determined the echo top position to be 0.00816 s. The
``coordinates_offset``, the time associated with the first point in the signal,
will need to be set to –0.00816 s. When correctly set, the time origin should
coincide with the maximum magnitude of the complex signal.

Second, you need to phase correct the time domain so that the maximum echo
amplitude is in the real part of the signal. For this operation, you can use
numpy `abs
() <https://numpy.org/doc/stable/reference/generated/numpy.absolute.html>`_ to
take the absolute value of each complex signal amplitude, and numpy `argmax
() <https://numpy.org/doc/stable/reference/generated/numpy.argmax.html>`_ to
find the time index where the absolute value of the signal is at a maximum.
Then use the signal phase at that time index to place the maximum amplitude
into the real part of the time domain signal.

Both these steps are performed by the code below.

.. skip: next

.. plot::
    :context: close-figs

    import numpy as np

    # set time origin to echo top
    csdm_ds.dimensions[0].coordinates_offset = "-0.00816 s"

    # Phase echo top, putting maximum amplitude into real part
    index = np.argmax(np.abs(csdm_ds.dependent_variables[0].components[0]))
    angle = np.angle(csdm_ds.dependent_variables[0].components[0][index])
    phased_ds = csdm_ds * np.exp(-1j * angle)

    plt.figure(figsize=(5, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(phased_ds.real, label="real")
    ax.plot(phased_ds.imag, label="imag")
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.show()

Here, you see that the echo top has been phased so that the maximum amplitude is
in the real (blue) part and that the echo top occurs at the time origin. Notice
that the echo has a slight asymmetry about the time origin after it has been
phased. The first half of the echo has a slightly stronger amplitude than the
last half. This asymmetry is due to an additional dephasing caused by
homonuclear dipolar couplings among the :math:`^{27}\text{Al}` nuclei. It may
have been possible to remove or minimize the effects of these dipolar couplings
using a higher MAS rate. Nonetheless, you can still proceed in this analysis
and, as you will see later, can model this additional decay with an ad-hoc
Gaussian convolution of the spectrum.

Next, create a SignalProcessor object to apply the Fourier transform operation
to the CSDM object ``exp_spectrum``. Note that with a correctly set time
origin, the :py:meth:`~mrsimulator.signal_processor.FFT` operation
automatically applies the appropriate first-order phase correction to the
spectrum after performing the fast Fourier transform. After performing the
Fourier transform, convert the coordinate units of the CSDM dimension from
frequency to a frequency ratio using the
`to()
<https://csdmpy.readthedocs.io/en/stable/api/Dimensions.html#csdmpy.Dimension.to>`_
method of the
`Dimension <https://csdmpy.readthedocs.io/en/stable/api/Dimensions.html>`_ object.

.. skip: next

.. plot::
    :context: close-figs

    from mrsimulator import signal_processor as sp

    ft = sp.SignalProcessor(operations=[sp.FFT()])
    exp_spectrum = ft.apply_operations(dataset=phased_ds)
    exp_spectrum.dimensions[0].to("ppm", "nmr_frequency_ratio")

    fig, ax = plt.subplots(1, 2, figsize=(9, 3.5), subplot_kw={"projection": "csdm"})
    ax[0].plot(exp_spectrum.real)
    ax[0].plot(exp_spectrum.imag)
    ax[0].set_title("Full Spectrum")
    ax[0].grid()
    ax[1].plot(exp_spectrum.real, label="real")
    ax[1].plot(exp_spectrum.imag, label="imag")
    ax[1].set_title("Zoomed Spectrum")
    ax[1].set_xlim(-15, 15)
    ax[1].grid()
    plt.tight_layout()
    plt.legend()
    plt.show()

.. Again, the blue and orange lines are the real and imaginary parts of the complex
.. frequency-domain spectrum.

.. _Measure Noise:

Measure Noise
-------------

Now that you have an adequately phased frequency domain dataset, you'll need to
take the real part of the spectrum for the rest of the analysis, i.e., remove
the imaginary part.

The least-squares analysis also needs the standard deviation of the noise in the
spectrum. We can obtain that from the spectrum regions below -20 ppm or above
20 ppm, where there is no signal amplitude. To accomplish this, you can use
numpy
`where() <https://numpy.org/doc/stable/reference/generated/numpy.where.html>`_. It
evaluates a condition for each item in the list and returns the indexes for
those items where the condition is true. With the indexes returned by
`where() <https://numpy.org/doc/stable/reference/generated/numpy.where.html>`_, you
can calculate the standard deviation of the noise region with numpy
`std() <https://numpy.org/doc/stable/reference/generated/numpy.std.html>`_.

.. skip: next

.. plot::
    :context: close-figs

    # Use only the real part of the spectrum
    exp_spectrum = exp_spectrum.real

    # Use region below -20 ppm to calculate the noise standard deviation
    loc = np.where(exp_spectrum.dimensions[0].coordinates < -20e-6)
    sigma = exp_spectrum[loc].std()

You can now move to the next step and create the fitting model.

Create Fitting Model
--------------------

To create a proper fitting model, you'll need more information about the nuclei
being observed, the material's phase, and some idea about the local structure
around the atoms holding the observed nuclei. In this example, you know that
you are working with :math:`^{27}\text{Al}`, a quadrupolar nucleus with a
half-integer spin of 5/2, and that the material, :math:`\text{Al(acac)$_3$}`,
is a solid polycrystalline sample. The symmetry of the
first-coordination sphere around aluminum is likely low enough to generate a
large electric field gradient, and hence a sizeable quadrupolar coupling
constant for :math:`^{27}\text{Al}`. These details are usually sorted out
before the NMR measurement and used to choose the appropriate NMR methods for
the sample. In this example, the measurement was performed under magic-angle
spinning at a rotation rate of 12.5 kHz. Due to the expected large quadrupolar
coupling, relatively low power rf pulses were used to excite only the
central :math:`m = \tfrac{1}{2}\rightarrow-\tfrac{1}{2}` transition of
:math:`^{27}\text{Al}`. The central transition is much narrower and more easily
detected than the other transitions.  Armed with this understanding of the
sample and method, you can proceed to create the fitting model.

Start by creating the Method object to model the experimental method used to
acquire the spectrum. Choose the
:py:meth:`~mrsimulator.method.lib.base.BlochDecayCTSpectrum()` method since the
measurement is designed to excite only the central transition of the
:math:`^{27}\text{Al}` nuclei. From the CSDM object holding the experimental
spectrum, i.e., ``exp_spectrum``, you can extract the relevant parameters for
the ``spectral_dimensions`` attribute of the
BlochDecayCTSpectrum method using the
fitting utility function
:py:meth:`~mrsimulator.utils.get_spectral_dimensions`. The experimental
measurement parameters associated with the method attributes
``magnetic_flux_density`` and ``rotor_frequency`` are also used in creating
this BlochDecayCTSpectrum method.
Finally, every Method object has the ``experiment`` attribute used to hold the
experimental spectrum that is to be modeled with the Method object.

.. skip: next

.. plot::
    :context: close-figs

    from mrsimulator.method.lib import BlochDecayCTSpectrum
    from mrsimulator.utils import get_spectral_dimensions

    spectral_dims = get_spectral_dimensions(exp_spectrum)
    MAS = BlochDecayCTSpectrum(
        channels=["27Al"],
        magnetic_flux_density=9.4,  # in T
        rotor_frequency=12500,  # in Hz
        spectral_dimensions=spectral_dims,
        experiment=exp_spectrum,  # add the measurement to the method.
    )


To build a spin system, you need to know how many magnetically inequivalent
nuclei are in the sample and if there are couplings between them. Inspection of
the spectrum reveals an anisotropic lineshape that appears to be characteristic
of the second-order MAS lineshape of a single site. Knowing this requires that
you are already familiar with such lineshapes (**mrsimulator** can help with
that!). One might also hypothesize that there may be other sites with lower
intensity present in the spectrum, or perhaps the spectrum, as noted earlier,
is from a distribution of :math:`^{27}\text{Al}` sites with very similar efg
tensor parameters and dipolar couplings among them. These are all valid
hypotheses and could be used to create more elaborate and perhaps even more
realistic spin system models. For now, you can choose the simplest spin system
model with a single
:math:`^{27}\text{Al}` site,  as shown in the code below.

.. skip: next

.. plot::
    :context: close-figs

    from mrsimulator import Site, SpinSystem, Simulator

    site = Site(
        isotope="27Al",
        isotropic_chemical_shift=5,
        quadrupolar={"Cq": 3e6, "eta": 0.0},
    )
    sys = SpinSystem(sites=[site])

The tensor parameters above are an educated guess for the tensor parameters,
which can be iteratively refined using the code that follows.


Create the simulator object initialized with the SpinSystem and Method objects
and run.

.. skip: next

.. plot::
    :context: close-figs

    sim = Simulator(spin_systems=[sys], methods=[MAS])
    sim.run()

Before comparing the simulation to the experimental spectrum, you need to add
the Gaussian line broadening to the simulation. Setup a SignalProcessor object
to do a Gaussian lineshape convolution with an FWHM of 50 Hz.

Additionally, you must scale the simulation in intensity to match the
experimental spectrum. You may have noticed in earlier plots that the vertical
axis of the experimental spectrum plot was on the order of 1e6. Use numpy
`max() <https://numpy.org/doc/stable/reference/generated/numpy.maximum.html>`_ to
get the highest amplitude, set that as the factor as a Scale operation in the
SignalProcessor.

.. skip: next

.. plot::
    :context: close-figs

    # Post Simulation Processing
    # --------------------------
    relative_intensity_factor = exp_spectrum.max() / sim.methods[0].simulation.max()
    processor = sp.SignalProcessor(operations=[
            sp.IFFT(),
            sp.apodization.Gaussian(FWHM="50 Hz"),
            sp.FFT(),
            sp.Scale(factor=relative_intensity_factor)
        ]
    )
    processed_dataset=processor.apply_operations(dataset=sim.methods[0].simulation).real


You now have set up and simulated the first guess in modeling the experimental
spectrum. Plot it and see how it compares to the experimental spectrum.

.. skip: next

.. plot::
    :context: close-figs

    # Plot of the guess spectrum
    # --------------------------
    plt.figure(figsize=(6, 3.0))
    ax = plt.subplot(projection="csdm")
    ax.plot(exp_spectrum, label="Experiment")
    ax.plot(processed_dataset, label="guess spectrum")
    ax.set_xlim(-15, 15)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()


The fit parameters are the spin system tensor and signal processor parameters.
If your initial guess is not so good, you could iteratively change the fit
parameters until your simulation is closer to the experimental spectrum. This
will ensure faster convergence to the best-fit parameters and could prevent the
least-squares analysis from falling into false minima on the chi-squared
surface. For this example, however, the above initial guess should be good enough


Perform Least-Squares Analysis
------------------------------

Up to this point in the discussion, you've done little more than what you've
learned earlier in setting up a simulation with mrsimulator. Except now,
you're ready to leverage the power of `LMFIT
<https://lmfit.github.io/lmfit-py/>`_ to obtain the best-fit parameters.

Define the fit parameters
'''''''''''''''''''''''''


Begin by using an **mrsimulator** utility function
:py:meth:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params` to extract a
list of LMFIT parameters from the Simulator and SignalProcessor objects.

.. skip: next

.. plot::
    :context: close-figs

    from mrsimulator.utils import spectral_fitting as sf
    fit_parameters = sf.make_LMFIT_params(sim, processor)
    print(fit_parameters.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

.. parsed-literal::

    Name                                      Value      Min      Max     Vary     Expr
    SP_0_operation_1_Gaussian_FWHM               50     -inf      inf     True     None
    SP_0_operation_3_Scale_factor          4.322e+06    -inf      inf     True     None
    sys_0_abundance                             100        0      100    False      100
    sys_0_site_0_isotropic_chemical_shift         5     -inf      inf     True     None
    sys_0_site_0_quadrupolar_Cq               3e+06     -inf      inf     True     None
    sys_0_site_0_quadrupolar_eta                  0        0        1     True     None
    None

The output of the ``print()`` statement, shown above, gives the table of the
LMFIT parameters created
by :py:meth:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params`. The
returned ``fit_parameters`` is a dictionary with each fit parameter object
identified by a string.  LMFIT does not allow  special characters
such as ``[``, ``]`` or ``.`` in the parameter string identifiers.  Therefore, when
the
:py:meth:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params` function
creates the LMFIT parameters dictionary, it flattens the variable namespace
into a string with these special characters replaced by a ``_``. For example,

**"sim.spin_systems[0].sites[1].quadrupolar.Cq"** :math:`\rightarrow`
**"sys_0_site_1_quadrupolar_Cq"**

or

**"sp[0].operation[3].factor"** :math:`\rightarrow` **"SP_0_operation_3_Scale_factor"**.

Using these parameter string names, you can access and change any of its LMFIT
parameter attributes, i.e.,
``value``, ``min``, ``max``, ``vary``, ``expr``. For example, using the code below, you
can set the quadrupolar asymmetry parameter value to be zero and request that
it be held constant during the fit.

.. skip: next

.. plot::
    :context: close-figs

    fit_parameters["sys_0_site_0_quadrupolar_eta"].value = 0
    fit_parameters["sys_0_site_0_quadrupolar_eta"].vary = False


.. warning::

    First-principles DFT calculations based on structural hypotheses can sometimes
    help determine the initial guess for some parameters; however, they are rarely
    accurate enough, even when using the correct structure, to be used as fixed
    parameters in a least-squares analysis of an experimental spectrum.



Define and minimize the chi-squared function
''''''''''''''''''''''''''''''''''''''''''''

To perform a least-squares analysis, `LMFIT
<https://lmfit.github.io/lmfit-py/>`_ needs a chi-squared function. LMFIT
expects this function to return a list of residuals (difference between model
and data) divided by the experimental noise standard deviation. Mrsimulator
comes with a pre-built chi-squared
function :py:meth:`~mrsimulator.utils.spectral_fitting.LMFIT_min_function`
which takes the Simulator, SignalProcessor, and the experimental noise standard
deviation as function arguments.


Perform the chi-squared minimization
''''''''''''''''''''''''''''''''''''

The least-squares analysis is performed by creating an `LMFIT
<https://lmfit.github.io/lmfit-py/>`_ `Minimizer
<https://lmfit-py.readthedocs.io/en/latest/fitting.html#lmfit.minimizer.Minimizer>`_
object initialized with a chi-squared function and the fit parameters
(``fit_parameters``). Any additional objects needed to evaluate the chi-squared
function are placed in ``fcn_args``.
For :py:meth:`~mrsimulator.utils.spectral_fitting.LMFIT_min_function`,
``fcn_args`` needs to hold the Simulator, SignalProcessor, and the experimental
noise standard deviation.

After the
`minimize() <https://lmfit-py.readthedocs.io/en/latest/fitting.html#lmfit.minimizer.minimize>`_
function of the
`Minimizer <https://lmfit-py.readthedocs.io/en/latest/fitting.html#lmfit.minimizer.Minimizer>`_
object exits, the parameters in the Simulator and SignalProcessor are updated
with the best-fit parameters and the results of the least-squares analysis is
returned as an
`MinimizerResult <https://lmfit-py.readthedocs.io/en/latest/fitting.html#lmfit.minimizer.MinimizerResult>`_
object containing the optimized parameters and several goodness-of-fit
statistics.

Use the code below to create and initialize the ``Minimizer`` object, run the
minimization, and print the
`MinimizerResult <https://lmfit-py.readthedocs.io/en/latest/fitting.html#lmfit.minimizer.MinimizerResult>`_.

.. skip: next

.. plot::
    :context: close-figs

    from lmfit import Minimizer
    minner = Minimizer(sf.LMFIT_min_function, fit_parameters, fcn_args=(sim, processor, sigma))
    result = minner.minimize()
    result


.. figure:: ../_static/FitStatistics1.*
    :width: 1200
    :alt: figure
    :align: center

From the printout of the
`MinimizerResult <https://lmfit-py.readthedocs.io/en/latest/fitting.html#lmfit.minimizer.MinimizerResult>`_
above, you can find the best-fit parameters and their associated uncertainties
from least-squares analysis.

.. warning::

    A word of caution about best-fit parameter uncertainties: If the model is
    accurate, then you expect the residuals to be pure noise, i.e., a histogram
    of the residuals should arise from a Gaussian parent distribution with a
    mean of zero. Therefore, at the very least, you should inspect a plot of
    the residuals and, even better, check that a histogram of the residuals is
    consistent with a Gaussian parent distribution.

    If this is not true, then
    the parameter uncertainties from the least-squares analysis will be
    underestimated. Such discrepancies between the experimental and simulated
    spectra can often arise from measurement artifacts, e.g., receiver
    deadtimes, non-uniform excitation, etc. They can also arise from an
    inadequate model (spin systems and method) for the spectrum.



Compare experimental and best-fit spectra with residuals
''''''''''''''''''''''''''''''''''''''''''''''''''''''''

You can now plot the experimental and best-fit simulated spectra along with the
residuals.  Use the **mrsimulator** utility
function :py:meth:`~mrsimulator.utils.spectral_fitting.bestfit`
and :py:meth:`~mrsimulator.utils.spectral_fitting.residuals` to extract the
best-fit simulation and the residuals as CSDM objects.

.. skip: next

.. plot::
    :context: close-figs

    best_fit = sf.bestfit(sim, processor)[0].real
    residuals = sf.residuals(sim, processor)[0].real

    # Plot the spectrum
    plt.figure(figsize=(6, 3.0))
    ax = plt.subplot(projection="csdm")
    ax.plot(exp_spectrum, label="Experiment")
    ax.plot(best_fit, alpha=0.75, label="Best Fit")
    ax.plot(residuals, alpha=0.75, label="Residuals")
    ax.set_xlim(-15, 15)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()


The Minimizer will improve the fit parameters even if the initial parameters are
far from the best-fit values. However, if the initial parameters are too far
away, the Minimizer may not reach the best-fit parameters in a single run. If
you think that may be the case, you can re-extract a new initial guess from the
Simulator and SignalProcessor objects using
:py:meth:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params`, create and
initialize a new Minimizer object as before, and run again, i.e., rerun the
code starting at the beginning of this section. You may see that the fit
improves and gives a lower chi-squared value.

In the least-square analysis above, you had locked the quadrupolar asymmetry
parameter to a value of zero, which is reasonably close to the true value. At
such low values, the quadrupolar asymmetry parameter is correlated to the
Gaussian line broadening FWHM in the fit. Set the quadrupolar asymmetry
parameter to be a fit parameter, and rerun the analysis.

.. skip: next

.. plot::
    :context: close-figs

    fit_parameters["sys_0_site_0_quadrupolar_eta"].value = 0
    fit_parameters["sys_0_site_0_quadrupolar_eta"].vary = True
    minner = Minimizer(sf.LMFIT_min_function, fit_parameters, fcn_args=(sim, processor, sigma))
    result = minner.minimize()
    best_fit = sf.bestfit(sim, processor)[0].real
    residuals = sf.residuals(sim, processor)[0].real

    # Plot the spectrum
    plt.figure(figsize=(6, 3.0))
    ax = plt.subplot(projection="csdm")
    ax.plot(exp_spectrum, label="Experiment")
    ax.plot(best_fit, alpha=0.75, label="Best Fit")
    ax.plot(residuals, alpha=0.75, label="Residuals")
    ax.set_xlim(-15, 15)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

.. figure:: ../_static/FitStatistics2.*
    :width: 1200
    :alt: figure
    :align: center

You see a slight improvement in the fit, with the asymmetry parameter increasing
from 0 to 0.136, and the Gaussian FWHM decreased from 106 to 63 Hz. The
MinimizerResult printout also shows a correlation of -0.74 between these two
parameters.

We close this section by noting that a compelling feature of mrsimulator& LMFit
is that you can perform a simultaneous spectra fit from different methods for a
single set of spin system parameters. Check out all the examples in
the :ref:`fitting_examples`, notably the
:ref:`sphx_glr_fitting_1D_fitting_plot_2_13C_glycine_multi_spectra_fit.py` example.


.. plot::
    :include-source: False

    import shutil

    shutil.rmtree("Al_acac")
