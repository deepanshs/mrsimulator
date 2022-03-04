.. _getting_started:

===============
Getting Started
===============

We have put together some introductory examples which outline the basic use of ``mrsimulator``.
For more detailed documentation on the usage of ``mrsimulator`` classes, see the
User Documentation section. Also, check out our :ref:`example_gallery` and
:ref:`fitting_examples`.

Spin System
-----------

An NMR spin system is an isolated system of sites (spins) and couplings. Spin systems
can include as many sites and couplings necessary to model a sample. For this
introductory example, we will create a coupled :math:`^1\text{H}` - :math:`^{13}\text{C}`
spin system.
First we will construct two :ref:`site_documentation` objects for the :math:`^1\text{H}` and
:math:`^{13}\text{C}` sites.

.. plot::
    :context: close-figs

    # Import the Site and SymmetricTensor classes
    from mrsimulator import Site
    from mrsimulator.spin_system.tensors import SymmetricTensor

    # Create the Site objects
    H_site = Site(isotope="1H")
    C_site = Site(
        isotope="13C",
        isotropic_chemical_shift=100.0,  # in ppm
        shielding_symmetric=SymmetricTensor(
            zeta=70.0,  # in ppm
            eta=0.5,
        ),
    )

We now have two variables, ``H_site`` and ``C_site``, which are :ref:`site_api` objects. ``H_site``
represents a proton site with no chemical shift. ``C_site`` represents a carbon-13 site with
a chemical shift of 100 ppm as well as a shielding component represented by :ref:`sy_api`
object. We parametrize tensors using the Haeberlen convention.

Next we will define a dipolar coupling by creating a :ref:`coupling_documentation` object.

.. plot::
    :context: close-figs

    # Import the Coupling class
    from mrsimulator import Coupling

    # Create the Coupling object
    coupling = Coupling(
        site_index=[0, 1],
        dipolar=SymmetricTensor(D=-2e4),  # in Hz
    )

Now we have all the pieces needed to create the spin system.
If you need to create an uncoupled spin system, simply omit the ``couplings`` attribute.

.. plot::
    :context: close-figs

    # Import the SpinSystem class
    from mrsimulator import SpinSystem

    # Create the SpinSystem object
    spin_system = SpinSystem(
        sites=[H_site, C_site],
        couplings=[coupling],
    )

Thats it! We have created a spin system whose spectrum is ready to be simulated.

Methods
-------

A :ref`method_documentation` object describes an NMR method. For this introduction, we will use
the :py:class:`~mrsimulator.methods.BlochDecaySpectrum` which is one of the pre-defined methods.
Some attributes of the method still need to be provided as seen below.

.. plot::
    :context: close-figs

    # Import the BlochDecaySpectrum class
    from mrsimulator.methods import BlochDecaySpectrum

    # Create a BlochDecaySpectrum object
    method = BlochDecaySpectrum(
        channels=["13C"],
        magnetic_flux_density=9.4,  # in T
        rotor_angle=54.735 * 3.14159 / 180,  # in rad (magic angle)
        rotor_frequency=3000,  # in Hz
        spectral_dimensions=[
            dict(
                count=2048,
                spectral_width=80e3,  # in Hz
                reference_offset=6e3,  # in Hz
                label=r"$^{13}$C resonances",
            )
        ],
    )

The variable ``method`` defines a Bloch decay MAS method for the :math:`^{13}\text{C}` channel.
A Bloch decay method only has one spectral dimension and this specific spectral dimension has
2048 points spanning 80 kHz with a reference offset of 6 kHz.

.. ((The method is looking at)) a the :math:`^{13}\text{C}` channel in a 9.4 tesla environment while the
.. sample spins at 3 kHz at the magic angle. We also have a single spectral dimension  which
.. defines a frequency dimension with 2048 points, spanning 80 kHz with a reference offset of
.. 6 kHz. :ref:`spec_dim_documentation`

Now all we need is to put our :ref:`spin_sys` and :ref:`method_api` objects together and simulate
the spectrum.

Simulator
---------

At the heart of ``mrsimulator`` is the :ref:`simulator_documentation` object which performs
the calculation of the NMR spectrum. Lets create the :ref:`simulator_api` object:

.. plot::
    :context: close-figs

    # Import the Simulator class
    from mrsimulator import Simulator

    # Create a Simulator object
    sim = Simulator()

Each :ref:`simulator_api` object holds a list of :ref:`spin_sys` objects and a list of :ref:`method_api`
objects. Below we add the spin system and method objects we previously defined:

.. plot::
    :context: close-figs

    # Add the SpinSystem and Method objects
    sim.spin_systems = [spin_system]
    sim.methods = [method]

Now to simulate the spectrum we need to call :py:meth:`~mrsimulator.Simulator.run`
on our :ref:`simulator_api` object.

.. plot::
    :context: close-figs

    sim.run()

The simulated spectrum is calculated and stored in the method object. Next we process and
plot the data

.. note:: In ``mrsimulator``, all resonance frequencies are calculated assuming the
    weakly-coupled (Zeeman) basis for the spin system.

Signal Processing
-----------------

``mrsimulator`` performs all calculations in the frequency domain, so plotting the dataset now
would show only delta functions. For this reason, we have the :ref:`signal_processing_documentation`
object which applies post-processing to the data after simulation.

Here we apply 200 Hz of exponential line broadening.

.. plot::
    :context: close-figs

    from mrsimulator import signal_processing as sp

    # Create the SignalProcessor object
    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(),
            sp.apodization.Exponential(FWHM="200 Hz"),
            sp.FFT(),
        ]
    )

    # Apply the processor to the simulation data
    processed_data = processor.apply_operations(data=sim.methods[0].simulation)

Each :ref:`signal_processing_api` object has a list of operations which are applied sequentially to
a dataset. For a comprehensive list of operations and how to use the signal processing object,
see the :ref:`signal_processing_documentation` documentation page.

Plotting the Data
-----------------

We end this example by using the `matplotlib <https://matplotlib.org/stable/>`_ Python library
to plot the simulated dataset.
:numref:`fig1-getting-started` depicts the plot of the simulated spectrum.

Below is the code used to generate the image:

.. skip: next

.. plot::
    :context: close-figs

    import matplotlib.pyplot as plt
    plt.figure(figsize=(6.6, 4))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(processed_data.real)
    ax.invert_xaxis()  # reverse x-axis
    plt.tight_layout(pad=0.1)
    plt.show()

.. .. _fig1-getting-started:
.. .. figure:: ../_static/getting_started.png
..     :figwidth: 75%

    A simulated MAS spectrum of :math:`^{13}\text{C}`.
