.. _fitting_example:

Least-Squares Fitting Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``mrsimulator`` can interact with various Python data science 
packages.  One such package, 
`LMFIT <https://lmfit.github.io/lmfit-py/>`_, can be used to perform non-linear 
least-squares analysis of experimental NMR spectra. 

Here, we illustrate the use of the mrsimulator objects to

- import and prepare an experimental dataset for the least-squares analysis,
- create a fitting model using Simulator and SignalProcessor objects,
- use the fitting model to perform a least-squares analysis,
- extract the model parameters with uncertainties, and
- plot the experimental spectrum along with the best-fit simulation and residuals.

Import Experimental Dataset
---------------------------

In this example, you will apply the least-squares fitting procedure to a 
:math:`^{27}\text{Al}` magic-angle spinning spectrum of :math:`\text{Al
(acac)$_3$}` measured with whole echo acquisition.

You will begin by importing an experimental dataset, measured on a 9.4 T Bruker
AVANCE III HD NMR spectrometer into the script. Bruker datasets are saved in
folders with a number as the folder name. In this case, that folder has been
transferred from the spectrometer and renamed "Al_acac". 

For our purposes, the folder was also compressed into a zip archive and uploaded
to an internet-accessible server.  You can used the code block below to download
the zip archive from the server and unzip it into the originally named folder.

.. plot::
    :context: close-figs

    import requests
    import zipfile
    from io import BytesIO
    file_ = "https://ssnmr.org/sites/default/files/mrsimulator/Al_acac3_0.zip"
    request = requests.get(file_)
    z = zipfile.ZipFile(BytesIO(request.content))
    z.extractall("Al_acac")


Now that the Bruker dataset folder is accessible to your Python code, you can use the
Python package `nmrglue <https://github.com/jjhelmus/nmrglue>`_ to convert the 
dataset into a `CSDM <https://csdmpy.readthedocs.io/en/stable/>`_ object.


.. plot::
    :context: close-figs

    import nmrglue as ng

    # initialize nmrglue converter object
    converter = ng.convert.converter()

    # read in the bruker dataset file
    dic, data = ng.bruker.read("Al_acac") 

    converter.from_bruker(dic,data)

    # convert to CSDM format
    csdm_ds = converter.to_csdm()

Now that the dataset is converted into a CSDM object plot the dataset to make
sure that it was imported correctly.

.. plot::
    :context: close-figs

    import matplotlib.pyplot as plt
    plt.figure(figsize = (5, 3))  # set the figure size
    ax = plt.subplot(projection = "csdm")
    ax.plot(csdm_ds.real)
    ax.plot(csdm_ds.imag)
    plt.tight_layout()
    plt.grid()
    plt.show()

This is the raw time-domain dataset, acquired using whole-echo acquisition. The
blue and orange lines are the real and imaginary parts of the complex
time-domain signal. If you're working with your own experimental dataset and have already 
processed it into the frequency domain, then you can skip the next few steps and proceed 
to the `Measure Noise`_ section.

Process Experimental Dataset
----------------------------

Proceeding from here, you'll need to transform this dataset into the frequency
domain for the least-squares analysis. Before applying the Fourier transform,
however, there are two things that need to be adjusted. 

First, you need to adjust the ``coordinates_offset`` to place the time origin at
the top of the echo. This time offset can be found among the pulse sequence
parameters. If the signal was acquired with a simple Hahn-echo sequence,
i.e., :math:`\pi/2-\tau-\pi-t`, then the ``coordinates_offset`` should be the
time between the centers of the two pulses. However, there are often some
additional receiver delays before the signal acquisition begins, and those
times need to be subtracted from the interpulse spacing. In this measurement,
we determine the echo top position to be 0.00816 s. The ``coordinates_offset``,
which is the time associated with the first point in the signal, will need to
be set to –0.00816 s. When correctly set, the time origin should coincide with 
the maximum magnitude of the complex signal.

Second, you need to phase correct the time domain so that the maximum echo
amplitude is in the real part of the signal. For this operation, you can use numpy
`abs() <https://numpy.org/doc/stable/reference/generated/numpy.absolute.html>`_ to
take the absolute value of each complex signal amplitude, and 
numpy `argmax() <https://numpy.org/doc/stable/reference/generated/numpy.argmax.html>`_ 
to find the time index where the absolute value of the signal is at a maximum. Then use 
the signal phase at that time index to place the maximum amplitude into the real part
of the time domain signal.

Both these steps are performed by the code below.

.. plot::
    :context: close-figs

    import numpy as np

    # set time origin to echo top
    csdm_ds.dimensions[0].coordinates_offset = "-0.00816 s"

    # Phase echo top, putting maximum amplitude into real part
    index = np.argmax(np.abs(csdm_ds.dependent_variables[0].components[0]))
    angle = np.angle(csdm_ds.dependent_variables[0].components[0][index])
    phased_ds = csdm_ds * np.exp(-1j*angle)

    plt.figure(figsize = (5, 3))  # set the figure size
    ax = plt.subplot(projection = "csdm")
    ax.plot(phased_ds.real)
    ax.plot(phased_ds.imag)
    plt.tight_layout()
    plt.grid()
    plt.show()

Here, you see that the echo top has been phased so that the maximum amplitude is
in the real (blue) part and that the echo top occurs at the time origin. 

Next, create a SignalProcessor object to apply the Fourier transform operation to the 
CSDM object ``exp_spectrum``.  Note that with a correctly set time origin, the ``FFT`` 
operation automatically applies the appropriate first-order phase correction 
to the spectrum after performing the fast Fourier transform.  After performing the Fourier
transform, convert the coordinate units of the CSDM dimension from frequency to a
frequency ratio using the 
`to() <https://csdmpy.readthedocs.io/en/stable/api/Dimensions.html#csdmpy.Dimension.to>`_ 
method of the `Dimension <https://csdmpy.readthedocs.io/en/stable/api/Dimensions.html>`_ object.

.. plot::
    :context: close-figs

    from mrsimulator import signal_processor as sp

    ft = sp.SignalProcessor(operations = [sp.FFT()])
    exp_spectrum = ft.apply_operations(dataset = phased_ds)
    exp_spectrum.dimensions[0].to("ppm", "nmr_frequency_ratio")

    fig, ax = plt.subplots(1, 2, figsize = (9, 3.5), subplot_kw = {"projection": "csdm"})
    ax[0].plot(exp_spectrum.real)
    ax[0].plot(exp_spectrum.imag)
    ax[0].set_title("Full Spectrum")
    ax[0].grid()
    ax[1].plot(exp_spectrum.real)
    ax[1].plot(exp_spectrum.imag)
    ax[1].set_title("Zoomed Spectrum")
    ax[1].set_xlim(-15,15)
    ax[1].grid()
    plt.tight_layout()
    plt.show()

Again, the blue and orange lines are the real and imaginary parts of the complex frequency-domain 
spectrum.

.. _Measure Noise:

Measure Noise
-------------

Now that you have an adequately phased frequency domain dataset, you'll need to take the
real part of the spectrum for the rest of the analysis, i.e., remove the imaginary 
part. 

The least-squares analysis also needs the standard deviation of the noise in the
spectrum. We can obtain that from the spectrum regions below -20 ppm or
above 20 ppm, where there is no signal amplitude.  To accomplish this, you can
use numpy `where() <https://numpy.org/doc/stable/reference/generated/numpy.where.html>`_.
It evaluates a condition for each item in the list, and return the indexes for those 
items where the condition is true.   With the indexes returned by 
`where() <https://numpy.org/doc/stable/reference/generated/numpy.where.html>`_,
the standard deviation of the noise region can be caculated with numpy 
`std() <https://numpy.org/doc/stable/reference/generated/numpy.std.html>`_.

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
around the atoms holding the observed nuclei. In this example, you know that you
are working with :math:`^{27}\text{Al}`, a quadrupolar nucleus with a
half-integer spin of 5/2, and that the material, :math:`\text{Al(acac)$_3$}`, is a solid
polycrystalline sample. The symmetry of the first-coordination sphere around
aluminum is likely low enough to generate a large electric field gradient, and
hence a sizeable quadrupolar coupling constant for :math:`^{27}\text{Al}`. These
details are usually sorted out before the NMR measurement and used to choose
the appropriate NMR methods for the sample. In this example, the measurement
was performed under magic-angle spinning at a rotation rate of 12.5 kHz. Due to
the expected large quadrupolar coupling, relatively low power rf pulses were
used to excite only the central :math:`m = \tfrac{1}{2}\rightarrow-\tfrac{1}
{2}` transition of :math:`^{27}\text{Al}`. The central transition is much narrower and
more easily detected than the other transitions.

Armed with this understanding of the sample and method, you can proceed to create
the fitting model. You might begin by setting up the spin system. But here again, you are
faced with needing more information about the nuclei being observed, i.e., you
need to know how many magnetically inequivalent nuclei are in the sample, and
if there are any couplings between nuclei. Inspection of the spectrum reveals
an anisotropic lineshape that appears to be characteristic of the second-order
MAS lineshape of a single site. Knowing this requires that you are already
familiar with such lineshapes (``mrsimulator`` can help with that!). One
might also hypothesize that there may be other sites with lower intensity
present in the spectrum, or perhaps the spectrum is from a distribution of
sites with very similar NMR tensor parameters as well as having dipole couplings
among them. These are all valid hypotheses and could be used to create more 
elaborate spin system models. For now, you can invoke Occam's razor and choose 
the simplest spin system model with a single :math:`^{27}\text{Al}` site,  as 
shown in the code below.

.. plot::
    :context: close-figs

    from mrsimulator import Site, SpinSystem, Simulator

    site = Site(
        isotope = "27Al",
        isotropic_chemical_shift = 5, 
        quadrupolar = {"Cq" : 3e6, "eta" : 0.0},
    )
    sys = SpinSystem(sites = [site]) 

The tensor parameters above are an educated guess for the tensor parameters, 
which can be iteratively refined using the code that follows.

Next, create the Method object to model the experimental method used to
acquire the spectrum. It is a straightforward procedure in this case. Choose
the ``BlochDecayCTSpectrum`` method since the measurement is designed to
excite only the central transition of the :math:`^{27}\text{Al}` nuclei. From
the CSDM object holding the experimental spectrum, i.e., ``exp_spectrum``, you
can extract the relevant parameters for the ``spectral_dimension`` attribute of
the ``BlochDecayCTSpectrum`` method using the fitting utility function
``get_spectral_dimensions()``. The experimental measurement parameters
associated with the method attributes ``magnetic_flux_density`` and
``rotor_frequency`` are also used in creating this ``BlochDecayCTSpectrum``
method. Finally, every Method object has ``experiment`` attribute used to hold
the experimental spectrum that is to be modeled with the Method object.

.. plot::
    :context: close-figs

    from mrsimulator.method.lib import BlochDecayCTSpectrum
    from mrsimulator.utils import get_spectral_dimensions

    spectral_dims = get_spectral_dimensions(exp_spectrum)
    MAS = BlochDecayCTSpectrum(
        channels = ["27Al"],
        magnetic_flux_density = 9.4,  # in T
        rotor_frequency = 12500,  # in Hz
        spectral_dimensions = spectral_dims,
        experiment = exp_spectrum,  # add the measurement to the method.
    )

Create the simulator object initialized with the SpinSystem and
Method objects and run.

.. plot::
    :context: close-figs

    sim = Simulator(spin_systems = [sys], methods = [MAS])
    sim.run()

Before comparing the simulation to the experimental spectrum, you need to add
some line broadening to the simulation.  Setup a SignalProcessor object to do
a Gaussian lineshape convolution with a FWHM of 50 Hz.

Additionally, the simulation needs to be scaled in intensity to
match the experimental spectrum. You may have noticed in earlier plots that the
vertical axis of the experimental spectrum plot was on the order of 1e6.  Use 
numpy `max() <https://numpy.org/doc/stable/reference/generated/numpy.maximum.html>`_  
to get the highest amplitude and set that as the factor as a Scale
operation in the SignalProcessor.

.. plot::
    :context: close-figs

    # Post Simulation Processing
    # --------------------------
    processor = sp.SignalProcessor(operations=[
            sp.IFFT(),
            sp.apodization.Gaussian(FWHM = "50 Hz"),
            sp.FFT(),
            sp.Scale(factor = exp_spectrum.max())
        ]
    )
    processed_dataset = processor.apply_operations(dataset = sim.methods[0].simulation).real


You now have set up and run a simulation of the first guess in modeling the experimental spectrum.
Plot it and see how it compares to the experimental spectrum.

.. plot::
    :context: close-figs

    # Plot of the guess spectrum
    # --------------------------
    plt.figure(figsize = (6, 3.0))
    ax = plt.subplot(projection="csdm")
    ax.plot(exp_spectrum.real, label = "Experiment")
    ax.plot(processed_dataset.real, label = "guess spectrum")
    ax.set_xlim(-15, 15)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()


The fit parameters are the spin system tensor and signal processor
parameters. If your initial guess was not so good, you could iteratively change
the fit parameters until your simulation is closer to the experimental spectrum.
This will ensure faster convergence to the best-fit parameters and could
prevent the least-squares analysis from falling into false minima on the
chi-squared surface.


Perform Least-Squares Analysis
------------------------------

Up to this point in the discussion, you've done little more than what you've
learned earlier in setting up a simulation with ``mrsimulator``. Except now,
you're ready to leverage the power of LMFIT to obtain the best-fit parameters.
Begin by using an ``mrsimulator`` utility function :py:meth:`~make_LMFIT_params` 
to extract a list of LMFIT parameters from the Simulator and SignalProcessor objects.

.. plot::
    :context: close-figs

    from mrsimulator.utils import spectral_fitting as sf
    params = sf.make_LMFIT_params(sim, processor)
    print(params.pretty_print(columns = ["value", "min", "max", "vary", "expr"]))

.. parsed-literal::

    Name                                      Value      Min      Max     Vary     Expr
    SP_0_operation_1_Gaussian_FWHM               50     -inf      inf     True     None
    SP_0_operation_3_Scale_factor           2.5e+06     -inf      inf     True     None
    sys_0_abundance                             100        0      100    False      100
    sys_0_site_0_isotropic_chemical_shift         5     -inf      inf     True     None
    sys_0_site_0_quadrupolar_Cq             2.9e+06     -inf      inf     True     None
    sys_0_site_0_quadrupolar_eta                0.2        0        1     True     None
    None

The output of the ``print()`` statement, shown above, gives the table of the
LMFIT parameters.  Here, you can determine which parameters are fit and which
are fixed.  

How can you change params attributes, e.g., Vary from True to False?


.. note::

    First-principles DFT calculations based on structural hypotheses can sometimes help determine 
    the initial guess for some parameters, however, they are rarely accurate enough–even when 
    using the correct structure–to be used as "ground-truth" fixed parameters in a least-squares 
    analysis of an experimental spectrum. 


The least-squares analysis is performed by creating a `LMFIT
<https://lmfit.github.io/lmfit-py/>`_ `Minimizer
<https://lmfit-py.readthedocs.io/en/latest/fitting.html#lmfit.minimizer.Minimizer>`_
object initialized with chi-squared function, the fit parameters (``params``), and 
any additional objects needed to evaluate the chi-squared function.  Here, you will 
set the chi-squared function to the ``mrsimulator`` utility function 
``sf.LMFIT_min_function`` and ``fcn_args`` to hold the Simulator, SignalProcessor, 
and the noise standard deviation of the experimental spectrum.

After minimize() exits the parameters in the Simulator and SignalProcessor will be updated, and the
results of the least-squares analysis are returned as an object containing the optimized parameters 
and several goodness-of-fit statistics in ``result``.  

.. plot::
    :context: close-figs

    from lmfit import Minimizer
    minner = Minimizer(sf.LMFIT_min_function, params, fcn_args=(sim, processor, sigma))
    result = minner.minimize()
    result


.. figure:: ../_static/FitStatistics.*
    :width: 800
    :alt: figure
    :align: center

You should also plot the experimental and simulated spectra along with the residuals.  Use
the ``mrsimulator`` utility function ``sf.bestfit()`` and ``sf.residuals()`` to extract
the best-fit simulation and the residuals as CSDM objects.

.. plot::
    :context: close-figs

    best_fit = sf.bestfit(sim, processor)[0]
    residuals = sf.residuals(sim, processor)[0]

    # Plot the spectrum
    plt.figure(figsize = (6, 3.0))
    ax = plt.subplot(projection = "csdm")
    ax.plot(exp_spectrum, label = "Experiment")
    ax.plot(best_fit, alpha=0.75, label = "Best Fit")
    ax.plot(residuals, alpha=0.75, label = "Residuals")
    ax.set_xlim(-15, 15)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()


The Minimizer will improve the fit parameters even if the initial parameters guess
is far from the best-fit values.  However, if the initial guess is too far away, the
Minimizer may not reach the best-fit parameters in a single run.  If you think that 
may be the case, you can re-extract a new initial guess from the Simulator and 
SignalProcessor objects using ``make_LMFIT_params()``, create and initialize a new 
Minimizer object as before, and run again, i.e., restart at the beginning of this
section.  You may see that the fit improves and obtain a lower chi-squared value.





We close this section by noting that a compelling feature of mrsimulator+LMFit
is that you can perform a simultaneous fit of spectra from different methods
for a single set of spin system parameters. Check out all the examples in
the :ref:`fitting_examples`.

.. plot::
    :include-source: False

    import shutil

    shutil.rmtree("Al_acac")

