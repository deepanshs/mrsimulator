
Post Simulation Signal Processing
=================================

Introduction
------------

After running a simulation, you may need to apply some post-simulation signal processing. For example, you may want to scale the intensities to
match the experiment, add line-broadening, or simulate signal artifact such as sinc
wiggles. There are many signal-processing libraries, such as Numpy and Scipy, that you
may use to accomplish this. Although, in NMR, certain operations like applying
line-broadening, is so regularly used that it soon becomes inconvenient to having to
write your own set of code. For this reason, the ``mrsimulator`` package offers some
frequently used signal-processing operations.

The following section will demonstrate the use of the
:class:`~mrsimulator.signal_processing.SignalProcessor` class in applying various
operations to the simulation data.

**Setup for the figures**

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> import matplotlib as mpl
    >>> import matplotlib.pyplot as plt
    ...
    >>> # global plot configuration
    >>> font = {"size": 11}
    >>> mpl.rc("font", **font)
    >>> mpl.rcParams["figure.figsize"] = [6, 3.5]


Simulating spectrum
-------------------
Please refer to the :ref:`using_objects` for a detailed description
of how to set up a simulation. Here, we will create a hypothetical simulation from two
single-site spin systems to illustrate the use of the post-simulation signal processing
module.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> from mrsimulator import Simulator, SpinSystem, Site
    >>> from mrsimulator.methods import BlochDecaySpectrum
    ...
    >>> # **Step 1** Create the site and spin system objects
    >>> site1 = Site(
    ...     isotope="29Si",
    ...     isotropic_chemical_shift=-75.0,  # in ppm
    ...     shielding_symmetric={"zeta": -60, "eta": 0.6},  # zeta in ppm
    ... )
    >>> site2 = Site(
    ...     isotope="29Si",
    ...     isotropic_chemical_shift=-80.0,  # in ppm
    ...     shielding_symmetric={"zeta": -70, "eta": 0.5},  # zeta in ppm
    ... )
    ...
    >>> sites = [site1, site2]
    >>> labels = ["Sys-1", "Sys-2"]
    >>> spin_systems = [SpinSystem(name=l, sites=[s]) for l, s in zip(labels, sites)]
    ...
    >>> # **Step 2** Create a Bloch decay spectrum method.
    >>> method_object = BlochDecaySpectrum(
    ...     channels=["29Si"],
    ...     magnetic_flux_density=7.1,  # in T
    ...     rotor_angle=54.735 * 3.1415 / 180,  # in rads
    ...     rotor_frequency=780,  # in Hz
    ...     spectral_dimensions=[
    ...         {
    ...             "count": 2048,
    ...             "spectral_width": 25000,  # in Hz
    ...             "reference_offset": -5070,  # in Hz
    ...             "label": "$^{29}$Si resonances",
    ...         }
    ...     ],
    ... )
    ...
    >>> # **Step 3** Create the simulation and add the spin system and method objects.
    >>> sim = Simulator()
    >>> sim.spin_systems = spin_systems
    >>> sim.methods = [method_object]
    ...
    >>> # **Step 4** Simulate the spectra.
    >>> sim.run()

The plot the spectrum is shown below.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> ax = plt.subplot(projection="csdm") # doctest: +SKIP
    >>> ax.plot(sim.methods[0].simulation, color="black", linewidth=1) # doctest: +SKIP
    >>> ax.set_xlim(-200, 50) # doctest: +SKIP
    >>> ax.invert_xaxis() # doctest: +SKIP
    >>> plt.tight_layout() # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

.. _fig1_signal_process:
.. figure:: _static/null.*

    1D :math:`^{29}\text{Si}` MAS simulation of two single-site spin system.

Post-simulating processing
--------------------------

Signal processing is a series of operations that are applied to the dataset. In this
workflow, the result from the previous operation becomes the input for the next
operation. In the ``mrsimulator`` library, we define this series as a list of operations.

Setting a list of operations
''''''''''''''''''''''''''''

All signal processing operations are located in the `signal_processing` module of the
``mrsimulator`` library. Within the module is the `apodization` sub-module. An
apodization is a point-wise multiplication operation of the input signal with the
apodizing vector. Please read our :ref:`operations_api` documentation for a complete
list of operations.

Import the module and sub-module as

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> import mrsimulator.signal_processing as sp
    >>> import mrsimulator.signal_processing.apodization as apo

In the following example, we show the application of a single operation—-convoluting
the frequency spectrum with a Gaussian lineshape, that is, simulating a Gaussian
line-broadening--using the :class:`~mrsimulator.signal_processing.SignalProcessor`
class.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> # list of processing operations
    >>> post_sim = sp.SignalProcessor(
    ...     operations=[
    ...         sp.IFFT(), apo.Gaussian(sigma=100), sp.FFT()
    ...     ]
    ... )

The required attribute of the ``SignalProcessor`` class, `operations`, holds the list of
operations that gets applied to the input dataset. The above set of operations is for a
frequency domain input signal undergoing a Gaussian convolution of 100 Hz. In this scheme,
the operations list will first perform an inverse Fourier Transform to convert
the frequency domain signal to the time domain. Next, the time domain signal is apodized
by a Gaussian function with a broadening factor of 100 Hz, followed by a forward Fourier
transformation transforming the signal back to the frequency domain.

.. note::
    For almost all NMR spectrum, the post-simulation processing is a convolution, including
    the line-broadening. The convolution theorem states that under suitable conditions, the
    Fourier transform of a convolution of two signals is the pointwise product of their
    Fourier transforms.


Applying operation to the spectrum
''''''''''''''''''''''''''''''''''

To apply the above list of operations to the simulation/input data, use the
:meth:`~mrsimulator.signal_processing.SignalProcessor.apply_operations` method of the
``SignalProcessor`` instance as follows,

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> processed_data = post_sim.apply_operations(data=sim.methods[0].simulation)

The `data` is the required argument of the `apply_operations` method, whose value is a
CSDM object holding the dataset. The variable `processed_data` holds the output, that is,
the processed data. The plot of the processed signal is shown below.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> ax = plt.gca(projection="csdm") # doctest: +SKIP
    >>> ax.plot(processed_data, color="black", linewidth=1) # doctest: +SKIP
    >>> ax.set_xlim(-200, 50) # doctest: +SKIP
    >>> ax.invert_xaxis() # doctest: +SKIP
    >>> plt.tight_layout() # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

.. _fig2_signal_process:
.. figure:: _static/null.*

    1D :math:`^{29}\text{Si}` MAS simulation of two single-site spin system with a
    100 Hz Gaussian convolution.

Applying operation to the sub-spectra
'''''''''''''''''''''''''''''''''''''

.. spectrum and follow up by decomposing the spectrum and processing each signal
.. independently.
.. The above code resulted in the same processing to be applied
.. to both signals because in the simulation the signals were not
.. seperated.

It is not uncommon for the NMR spectrum to compose of sub-spectrum, from different
sites/systems, exhibiting differential relaxations, and therefore, have different
extents of line-broadening. The reason for this differential relaxation behavior is
not the focus of this sub-section. Here, we show how one can simulate such spectra
using the operations list.

Before we can move forward, you will first need to identify these sub-systems and
simulate individual spectra for these systems. In this example, we will treat the two
spin systems as the two different spin environments exhibiting different
relaxations/line-broadening. To simulate the sub-spectrum from the individual
spin systems, modify the value of the :attr:`~mrsimulator.Simulator.config` attribute
as follows, and re-run the simulation.
Refer to the :ref:`config_simulator` section for further details.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> sim.config.decompose_spectrum = "spin_system"
    >>> sim.run()

.. Note, in the previous example, both sites/spin systems got the same extent of Gaussian
.. line-broadening. The following example illustrates how you can apply you might want to apply a different set of
.. In order to apply different processes to each signal,
.. we must set the simulation config to decompose the spectrum.
.. Steps 1-3 will be the same and we will start at step 4.
.. #
.. **Step 4** Decompose spectrum and run simulation.
.. sim.config.decompose_spectrum = "spin_system"
.. sim.run()
..  plt.xlabel("$^{29}$Si frequency / ppm")
..  plt.xlim(x.value.max(), x.value.min())
..  plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)

The above code generates two spectra, each corresponding to a spin system.
The plot of the spectra is shown below.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> ax = plt.gca(projection="csdm") # doctest: +SKIP
    >>> ax.plot(sim.methods[0].simulation) # doctest: +SKIP
    >>> ax.set_xlim(-200, 50) # doctest: +SKIP
    >>> ax.invert_xaxis() # doctest: +SKIP
    >>> plt.tight_layout() # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

.. _fig3_signal_process:
.. figure:: _static/null.*

    Two 1D :math:`^{29}\text{Si}` MAS simulations, shown in blue and organe, for the two
    single-site spin systems.

Because the simulation is stored as a CSDM [#f1]_ object, each sub-spectrum is a
dependent-variable of the CSDM object, sharing the same frequency dimension.
When using the list of the operations, you may selectively apply a given operation to a
specific dependent-variable by specifying the index of the corresponding
dependent-variable as an argument to the operation class. Note, the order of the
dependent-variables is the same as the order of the spin systems. Use the `dep_var_indx`
argument of the operation to specify the index. Consider the following list of
operations.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> post_sim = sp.SignalProcessor(
    ...     operations=[
    ...         sp.IFFT(), # convert to time-domain
    ...         apo.Gaussian(sigma=50, dep_var_indx=0),
    ...         apo.Exponential(FWHM=200, dep_var_indx=1),
    ...         sp.FFT(), # convert to frequency-domain
    ...     ]
    ... )

The above operations list first applies an inverse Fourier transformation,
followed by a Gaussian apodization on the dependent variable at index 0 (spin system
labeled as `sys1`), followed by an Exponential apodization on the dependent
variable at index 1 (spin system labeled as `sys2`), and finally a forward Fourier
transform. Note, the FFT and IFFT operations apply on all dependent-variables.

As before, apply the operations with the
:meth:`~mrsimulator.signal_processing.SignalProcessor.apply_operations` method.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> processed_data = post_sim.apply_operations(data=sim.methods[0].simulation)

The plot of the processed spectrum is shown below.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> ax = plt.gca(projection="csdm") # doctest: +SKIP
    >>> ax.plot(processed_data, alpha=0.9)  # doctest: +SKIP
    >>> ax.set_xlim(-200, 50) # doctest: +SKIP
    >>> ax.invert_xaxis() # doctest: +SKIP
    >>> plt.tight_layout()  # doctest: +SKIP
    >>> plt.show()  # doctest: +SKIP

.. _fig4_signal_process:
.. figure:: _static/null.*

    Two 1D :math:`^{29}\text{Si}` MAS simulations, shown in blue and organe, for the two
    single-site spin systems with a 50 Hz Gaussian and 200 Hz Lorentzian convolution,
    respectively.

Serializing the operations list
-------------------------------

You may also serialize the operations list using the
:meth:`~mrsimulator.signal_processing.SignalProcessor.to_dict_with_units`
method, as follows

.. doctest::

    >>> from pprint import pprint
    >>> pprint(post_sim.to_dict_with_units())
    {'operations': [{'dim_indx': 0, 'function': 'IFFT'},
                    {'dep_var_indx': 0,
                     'dim_indx': 0,
                     'function': 'apodization',
                     'sigma': '50.0 Hz',
                     'type': 'Gaussian'},
                    {'FWHM': '200.0 Hz',
                     'dep_var_indx': 1,
                     'dim_indx': 0,
                     'function': 'apodization',
                     'type': 'Exponential'},
                    {'dim_indx': 0, 'function': 'FFT'}]}

.. [#f1] Srivastava, D. J., Vosegaard, T., Massiot, D., Grandinetti, P. J.,
            Core Scientific Dataset Model: A lightweight and portable model and
            file format for multi-dimensional scientific data, PLOS ONE,
            **15**, 1-38, (2020).
            `DOI:10.1371/journal.pone.0225953 <https://doi.org/10.1371/journal.pone.0225953>`_
