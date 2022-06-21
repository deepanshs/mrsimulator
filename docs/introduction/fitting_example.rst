.. _fitting_example:

Least-Squares Fitting Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``mrsimulator`` can interact with a variety of other Python data science 
packages.  One such package, 
`LMFIT <https://lmfit.github.io/lmfit-py/>`_, can be used to perform non-linear 
least-squares analysis of experimental NMR spectra. 

Here, we illustrate the use of the mrsimulator objects to

- import and prepare the experimental dataset for the least-squares analysis,
- create a fitting model using Simulator and SignalProcessor objects,
- use the fitting model to perform a least-squares analysis,
- extract the model parameters with uncertainties, and
- plot the experimental spectrum along with the best-fit simulation and residuals.

Importing and processing the experimental dataset
-------------------------------------------------

In this example, we apply the least-squares fitting procedure to the 
:math:`^{27}\text{Al}` magic-angle spinning spectrum of :math:`\text{Al(acac)$_2$}`
measured with whole echo acquisition.

We begin by importing the experimental datatset, measured on a 9.4 T
Bruker AVANCE IIIHD NMR spectrometer, into the script.  Bruker datasets are 
saved in folders with a number as the folder name.  In this case, that folder 
has been transferred from the spectrometer and renamed to "Al_acac".  

For our purposes, the folder was also compressed into zip archive and uploaded to an
internet accessible server.  If you already have the folder available on your
local machine you can skip this step.

.. plot::
    :context: close-figs

    import requests
    import zipfile
    from io import BytesIO
    file_ = "https://ssnmr.org/sites/default/files/mrsimulator/Al_acac3_0.zip"
    request = requests.get(file_)
    z = zipfile.ZipFile(BytesIO(request.content))
    z.extractall("Al_acac")

Once we have Bruker dataset folder is accessible to our Python code, we use the
Python package `nmrglue <https://github.com/jjhelmus/nmrglue>`_ to convert the 
dataset into `CSDM <https://csdmpy.readthedocs.io/en/stable/>`__ format.

.. plot::
    :context: close-figs

    import nmrglue as ng

    # initialize nmrglue converter object
    converter = ng.convert.converter()

    # read in the bruker data file
    dic, data = ng.bruker.read("Al_acac") 

    converter.from_bruker(dic,data, remove_digital_filter=False)

    # convert to CSDM format
    csdm_ds = converter.to_csdm()

Now that the dataset is converted into a CSDM object, let's plot the
dataset to make sure that it was imported correctly.

.. plot::
    :context: close-figs

    import matplotlib.pyplot as plt
    plt.figure(figsize=(5, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(csdm_ds.real)
    ax.plot(csdm_ds.imag)
    plt.tight_layout()
    plt.grid()
    plt.show()

This is the raw time-domain dataset, acquired using whole-echo acquisition.
The blue and orange lines are the real and imaginary parts, respectively,
of the complex time domain signal.  If you've already processed your dataset
into frequency domain, then you can skip the next few steps and proceed to 
Creating the Fitting Model.

Proceeding from here, we'll need to transform this dataset into 
the frequency domain for the least-squares analysis.  Before applying the Fourier 
transform, however, we'll need to adjust the  ``coordinates_offset`` to place the 
time origin at the top of the echo. This time offset can be found among the pulse 
sequence parameters.  If the signal was acquired with a simple Hahn-echo sequence, 
i.e., :math:`\pi/2-\tau-\pi-t`, then the ``coordinates_offset`` should be 
the time between the centers of the two pulses.  However, there is often some 
additional receiver delays before the signal acquisition begins and those 
times need to be subtracted from the interpulse spacing.   In this measurement,
we determine the echo top position to be 0.00816 s.  The ``coordinates_offset``, 
which will be the time associated with the first point in the signal, will be –0.00816 s.  
When properly set, the time origin should coincide with the maximum magnitude of the complex
signal.

Additionally, we need to phase correct the time domain so that the maximum echo amplitude 
is in the real part of the signal.  For this operation, we'll use numpy ``max()`` to find 
the time where the magnitude of the signal is at a maximum, and then use that signal phase 
to place the maximum amplitude into the real part of the time domain signal.

.. plot::
    :context: close-figs

    import numpy as np

    # set time origin to echo top
    csdm_ds.dimensions[0].coordinates_offset = "-0.00816 s" 
    phased_ds = csdm_ds*np.exp(-1j*(np.pi+np.angle(csdm_ds.max()).value))
    plt.figure(figsize=(5, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(phased_ds.real)
    ax.plot(phased_ds.imag)
    plt.tight_layout()
    plt.grid()
    plt.show()

Here, you see that the echo top has been phased so that the maximum amplitude is in the real part,
and that the echo top occurs at the time origin.   With a properly set time origin, the Fourier 
transform operation can apply the appropriate first-order phase correction to the spectrum after
performing the fast Fourier transform, as shown in the code below.

.. plot::
    :context: close-figs

    from mrsimulator import signal_processing as sp

    ft = sp.SignalProcessor(operations=[sp.FFT()])
    exp_spectrum = ft.apply_operations(data=phased_ds)
    exp_spectrum.x[0].to("ppm", "nmr_frequency_ratio")

    fig, ax = plt.subplots(1, 2, figsize=(9, 3.5), subplot_kw={"projection": "csdm"})
    ax[0].plot(exp_spectrum.real)
    ax[0].plot(exp_spectrum.imag)
    ax[0].set_title("Full Spectrum")
    ax[0].grid()
    ax[1].plot(exp_spectrum.real)
    ax[1].plot(exp_spectrum.imag)
    ax[1].set_title("Zoomed Spectrum")
    ax[1].set_xlim(-20,20)
    ax[1].grid()
    plt.tight_layout()
    plt.show()

Now that we have a properly phased frequency domain dataset, we use only the real part of the spectrum
in the analysis, i.e., remove the imaginary part.  Additionally, the least-squares analysis also 
needs the standard deviation of the noise in the spectrum.  We can obtain that from the regions
of the spectrum from -40 to -10 ppm and from 10 to 40 ppm, where there is no signal amplitude.

.. plot::
    :context: close-figs

    exp_spectrum = exp_spectrum.real
    sigma = 0.03 #need code here to determine sigma

We can now move to the next step and create the fitting model.

Creating the Fitting Model
--------------------------

NMR spectra are like dog breeds; each can appear and behave quite differently. To
create a proper fitting model, we need more information about the nuclei being observed,
the material's phase, and some idea about the local structure around the atoms
holding the observed nuclei. In this example, we know that we are working with :math:`^{27}\text{Al}`, 
a quadrupolar nucleus with a half-integer spin of 5/2. The material, :math:`\text{Al(acac)$_2$}`, 
is a solid polycrystalline sample. The symmetry of the first-coordination sphere around aluminum
is likely low enough to generate a large electric field gradient, and hence sizeable quadrupolar
coupling constant for :math:`^{27}\text{Al}`. These details are usually sorted out before the
NMR measurement and used to choose the appropriate NMR methods for the sample. In
this example, the measurement was performed under magic-angle spinning at a rotation rate of 12.5 kHz.
Due to the expected large quadrupolar coupling, relatively low power rf pulses were used to excite
only the central :math:`m = \tfrac{1}{2}\rightarrow-\tfrac{1}{2}` transition of :math:`^{27}\text{Al}`.
This transition is much narrower and more easily detected than the other single-quantum transitions.

Armed with this understanding of the sample and method, we can proceed to create the fitting model.
We begin by setting up the spin system. Here again, we are faced with needing more information about
the nuclei being observed, i.e., we need to know how many magnetically inequivalent nuclei are 
in the sample.  Inspection of the spectrum reveals an anisotropic lineshape that appears to
be characteristic of the second-order MAS lineshape of a single site.  Knowing this requires that you
are already familiar with such lineshapes (something that ``mrsimulator`` can help with!).  One might
also hypothesize that there may be other sites with lower intensity present in the spectrum, or perhaps 
that the spectrum is from a distribution of sites with very similar NMR tensor parameters. These 
are all valid hypotheses and could be used to create more elaborate spin system models.  For now, we
invoke Occam's razor and choose the simplest spin system model with a single :math:`^{27}\text{Al}` site, 
as shown in the code below.

.. plot::
    :context: close-figs

    from mrsimulator import Site, SpinSystem, Simulator

    site = Site(
        isotope="27Al",
        isotropic_chemical_shift=5, 
        quadrupolar = {"Cq":2.9e6, "eta":0.2},
    )
    sys = SpinSystem(sites = [site]) 

We used an educated guess for the tensor parameters, which can be iteratively refined using the code that
follows.

Next, we create the Method object to model the experimental method used to acquire the spectrum. It is a
straightforward procedure in this case. We choose the ``BlochDecayCTSpectrum`` method since the measurement was
designed to excite only the central transition of the :math:`^{27}\text{Al}` nuclei. From the CSDM object holding
the experimental spectrum, i.e., ``exp_spectrum``, we can extract the relevant parameters for the ``spectral_dimension``
attribute of the ``BlochDecayCTSpectrum`` method using the fitting utility function ``get_spectral_dimensions()``.
The experimental measurement parameters associated with the method attributes ``magnetic_flux_density`` 
and ``rotor_frequency`` are also used in creating this ``BlochDecayCTSpectrum`` method. Finally, every Method object
has ``experiment`` attribute used to hold the experimental spectrum that is to be modeled with the Method object.

Next, the simulator object is created and initialized with the SpinSystem and Method objects, and run.

.. plot::
    :context: close-figs

    from mrsimulator.method.lib import BlochDecayCTSpectrum
    from mrsimulator.utils import get_spectral_dimensions

    spectral_dims = get_spectral_dimensions(exp_spectrum)
    MAS = BlochDecayCTSpectrum(
        channels=["27Al"],
        magnetic_flux_density=9.4,  # in T
        rotor_frequency=12500,  # in Hz
        spectral_dimensions= spectral_dims,
        experiment=exp_spectrum,  # add the measurement to the method.
    )
    sim = Simulator(spin_systems=[sys], methods=[MAS])
    sim.run()

Before comparing the simulation to the experimental spectrum, we need to add some line broadening to the 
simulation in the form of a Gaussian lineshape convolution.  Additionally, the simulation needs to be
scaled in intensity to match that of the experimental spectrum.  These two operations are performed using
the SignalProcessor object created in the code below.   The final spectrum, intended to model the 
experimental spectrum, is plotted after the SignalProcessor object has operated on the simulated spectrum.

.. plot::
    :context: close-figs

    # Post Simulation Processing
    # --------------------------
    processor = sp.SignalProcessor(operations=[
            sp.IFFT(),
            sp.apodization.Gaussian(FWHM="50 Hz"),
            sp.FFT(),
            sp.Scale(factor=2.5e6)
        ]
    )
    processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

    # Plot of the guess spectrum
    # --------------------------
    plt.figure(figsize=(6, 3.0))
    ax = plt.subplot(projection="csdm")
    ax.plot(exp_spectrum.real, "k", linewidth=1, label="Experiment")
    ax.plot(processed_data.real, "b",  linewidth=1, label="guess spectrum") #alpha=0.75,
    ax.set_xlim(-20, 20)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

Above is the experimental spectrum along with the simulation using our initial guesses for the fit parameters, 
i.e., the spin system tensor and signal processor parameters.  If our initial guess was not so good, then we 
would iteratively change the fit parameters until our simulation is reasonably close to the experimental 
spectrum.  This is done to ensure faster convergence to the best-fit parameters and could prevent the
least-squares analysis from falling into false minima on the chi-squared surface.

Up to this point in the discussion, we've done little more than what we've learned earlier in setting up a 
simulation with ``mrsimulator``.  Except now, we're ready to leverage the power of LMFIT to obtain the 
best-fit parameters.  We begin by using an ``mrsimulator`` utility function ``make_LMFIT_params()`` for 
extracting a list of LMFIT parameters from extracting the Simulator and SignalProcessor objects.

.. plot::
    :context: close-figs

    from mrsimulator.utils import spectral_fitting as sf
    params = sf.make_LMFIT_params(sim, processor)
    print(params.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

.. plot::
    :context: close-figs

    from lmfit import Minimizer
    minner = Minimizer(sf.LMFIT_min_function, params, fcn_args=(sim, processor, sigma))
    result = minner.minimize()
    result

We close this section by noting that a particularly powerful feature of mrsimulator+LMFit is that you can perform a simultaneous fit of spectra 
from different  methods for a single set of spin system parameters. Check out all the examples in the :ref:`fitting_examples`.

.. plot::
    :include-source: False

    import shutil

    shutil.rmtree("Al_acac")


