.. _getting_started:

===============
Getting Started
===============

In ``mrsimulator``, the user initializes objects from three mrsimulator classes: :ref:`spin_system_documentation`,
:ref:`method_documentation`, and :ref:`simulator_documentation`. :ref:`spin_system_documentation` defines the
spin system tensor parameters used to generate a particular subspectrum, and :ref:`method_documentation` defines the
parameters for the particular NMR measurement to be simulated. A list of :ref:`method_documentation` and
:ref:`spin_system_documentation` objects are used to initialize a Simulator object, which is then used to generate
the corresponding NMR spectra--returned as a CSDM object in each Method object. For more information on the CSDM
(Core Scientific Dataset Model), see the `csdmpy documentation <https://csdmpy.readthedocs.io/en/stable/>`__. There
is an additional class, :ref:`signal_processing_documentation`, for applying various post-simulation signal processing
operations to CSDM dataset objects. All objects can be serialized. We adopt the
`Javascript Object Notation (JSON) <https://www.json.org>`__ as the file-serialization format for the model because it
is human-readable if properly organized and easily integrable with numerous programming languages and related software
packages. It is also the preferred serialization for data exchange in web-based applications.

Here, we have put together some introductory examples which outline the basic use of ``mrsimulator``.
See the User Documentation section for more detailed documentation on the usage of ``mrsimulator`` classes.
Also, check out our :ref:`example_gallery` and :ref:`fitting_examples`.

Spin System
-----------

An NMR spin system is an isolated system of sites (spins) and couplings. Spin systems can include as many sites and couplings
as necessary to model a sample. For this introductory example, we will create a coupled :math:`^1\text{H}` - :math:`^{13}\text{C}`
spin system.  First we will construct two :ref:`site_documentation` objects for the :math:`^1\text{H}` and
:math:`^{13}\text{C}` sites.

.. plot::
    :context: reset

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
    my_sites=[H_site, C_site]

We now have two variables, ``H_site`` and ``C_site``, which are :ref:`site_api` objects. ``H_site``
represents a proton site with zero (default) chemical shift. ``C_site`` represents a carbon-13 site with
a chemical shift of 100 ppm as well as a shielding component represented by :ref:`sy_api`
object. We parametrize tensors using the Haeberlen convention. A Site object has default values for unspecified attributes.
All spin interaction parameters, e.g., isotropic chemical shift and other coupling parameters, are initialized to zero.
Additionally, the default isotope is ``1H``. For example, the code above could have used ``H_site = Site()``.

Next, we will define a dipolar coupling by creating a :ref:`coupling_documentation` object.

.. plot::
    :context: close-figs

    # Import the Coupling class
    from mrsimulator import Coupling

    # Create the Coupling object
    coupling = Coupling(
        site_index=[0, 1],
        dipolar=SymmetricTensor(D=-2e4),  # in Hz
    )

Couplings between Sites are specified using the indexes of the Sites in the list variable ``my_sites``.  
Now we have all the pieces needed to create the spin system.  If you need to create an uncoupled spin system, omit the ``couplings`` attribute.

.. plot::
    :context: close-figs

    # Import the SpinSystem class
    from mrsimulator import SpinSystem

    # Create the SpinSystem object
    spin_system = SpinSystem(
        sites=my_sites,
        couplings=[coupling],
    )

That's it! We have created a spin system whose spectrum is ready to be simulated.

Methods
-------

A :ref:`method_documentation` object describes an NMR method. For this introduction, we will use
the :py:class:`~mrsimulator.method.lib.BlochDecaySpectrum`, which is one of the pre-defined methods.
Some attributes of the Method need to be provided, as shown below.

.. plot::
    :context: close-figs

    # Import the BlochDecaySpectrum class
    from mrsimulator.method.lib import BlochDecaySpectrum
    from mrsimulator.method import SpectralDimension

    # Create a BlochDecaySpectrum object
    method = BlochDecaySpectrum(
        channels=["13C"],
        magnetic_flux_density=9.4,  # in T
        rotor_angle=54.735 * 3.14159 / 180,  # in rad (magic angle)
        rotor_frequency=3000,  # in Hz
        spectral_dimensions=[
            SpectralDimension(
                count=2048,
                spectral_width=80e3,  # in Hz
                reference_offset=6e3,  # in Hz
                label=r"$^{13}$C resonances",
            )
        ],
    )

The variable ``method`` defines a Bloch decay MAS method for the :math:`^{13}\text{C}` channel.
A Bloch decay method only has one spectral dimension, and this specific spectral dimension has
2048 points, spanning 80 kHz with a reference offset of 6 kHz.

.. ((The method is looking at)) a the :math:`^{13}\text{C}` channel in a 9.4 tesla environment while the
.. sample spins at 3 kHz at the magic angle. We also have a single spectral dimension which
.. defines a frequency dimension with 2048 points, spanning 80 kHz with a reference offset of
.. 6 kHz. :ref:`spec_dim_documentation`

Next, we put the SpinSystem and Method objects together to simulate the spectrum.

Simulator
---------

At the heart of ``mrsimulator`` is the :ref:`simulator_documentation` object, which calculates the NMR
spectrum. Let us first create a :ref:`simulator_api` object, initialized with our previously defined spin 
system and method, and then call :py:meth:`~mrsimulator.Simulator.run` on our :ref:`simulator_api` object.

.. plot::
    :context: close-figs

    # Import the Simulator class
    from mrsimulator import Simulator

    # Create a Simulator object
    sim = Simulator(spin_systems = [spin_system], methods = [method])
    sim.run()

The simulated spectrum is calculated and stored in the Method object.

.. note::
    In ``mrsimulator``, all resonance frequencies are calculated assuming the
    weakly-coupled (Zeeman) basis for the spin system.

Signal Processing
-----------------

``mrsimulator`` performs all calculations in the frequency domain.  Plotting the spectrum in this example would
show only delta functions. For this reason, we use the :ref:`signal_processing_documentation` object to add line
broadening to the simulated spectrum.  Below, we create a SignalProcessing object to do a convolution of the simulated
spectrum with a Lorentzian distribution with a full-width-half-maximum of 200 Hz.  This is performed in the time
domain by first applying an inverse fast Fourier transform, an apodization with an exponential decay, followed by
a fast Fourier transform back into the frequency domain.

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

A :ref:`signal_processing_api` object holds a list of operations applied sequentially to a dataset.
For a comprehensive list of operations and further details on using the :ref:`signal_processing_api` object,
see the :ref:`signal_processing_documentation` documentation page.

Plotting the Simulation
-----------------------

We end this example by using the Python package `matplotlib <https://matplotlib.org/stable/>`_
to plot the simulated dataset.  Below is code that can be used to generate an image and a pdf
file of the simulated spectrum:

.. _fig1-getting-started:
.. skip: next

.. plot::
    :context: close-figs
    :caption: A simulated :math:`^{13}\text{C}` MAS spectrum.

    import matplotlib.pyplot as plt
    plt.figure(figsize=(5, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(processed_data.real)
    ax.invert_xaxis()  # reverse x-axis
    plt.tight_layout(pad=0.1)
    plt.savefig("spectrum.pdf")
    plt.show()

The ``plt.savefig("spectrum.pdf")`` line creates a pdf file that can be edited in a vector graphics editor such as
Adobe Illustrator.

Saving the Simulation dataset
-----------------------------
``mrsimulator`` is designed to be part of a larger data workflow involving other software packages. 
For this larger context, ``mrsimulator`` uses the Core Scientific Dataset Model (CSDM) for importing 
and exporting your datasets. CSDM is a lightweight, portable, human-readable, and versatile standard 
for intra- and interdisciplinary exchange of scientific datasets. The model supports multi-dimensional 
datasets with a multi-component dependent variable discretely sampled at unique points in a 
multi-dimensional independent variable space. It can also hold correlated datasets assuming the different 
physical quantities (dependent variables) are sampled on the same orthogonal grid of independent variables. 
The CSDM can also serve as a re-usable building block in the development of more sophisticated portable 
scientific dataset file standards.

``mrsimulator`` also uses CSDM as its object model for simulated and experimental datasets. Any CSDM object in 
``mrsimulator`` can be serialized as a JavaScript Object Notation (JSON) file using its ``save()`` method. 
For example, the simulation after the signal processing step above is saved as a csdf file as shown below.

.. plot::
    :context: close-figs

    processed_data.save("processed_simulation.csdf")

For more information on the CSDM format, see the `csdmpy documentation <https://csdmpy.readthedocs.io/en/stable/>`__.

.. plot::
    :include-source: False

    import os
    from os.path import isfile

    if isfile("spectrum.pdf"): os.remove("spectrum.pdf")
    if isfile("processed_simulation.csdf"): os.remove("processed_simulation.csdf")
