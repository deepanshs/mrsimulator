.. _getting_started_coupled_spin_system_etoh:

==================================
Coupled Spin System: Using objects
==================================

In this example, we will simulate the :math:`^1\text{H}` NMR spectrum of
ethanol using the core ``mrsimulator`` objects. Let’s start by importing
all the necessary packages.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> import matplotlib.pyplot as plt
    >>> from mrsimulator import Simulator, SpinSystem, Site, Coupling
    >>> from mrsimulator.methods import BlochDecaySpectrum

Setting up coupled SpinSystem objects
-------------------------------------

Sites
^^^^^

An NMR spin system is defined by an isolated set of sites (spins) and couplings. You can
make a spin system as large and complex as needed, but for this example, we will build
the most abundant isotopomer of the ethanol molecule (all carbons are :math:`^{12}\text{C}`,
and the oxygen is :math:`^{16}\text{O}`). We start by defining the three distinct proton
sites and then build a list to hold all the sites (site indices correspond to atoms as
shown in the structure below).

.. figure:: _static/iso1.*
    :width: 200
    :alt: image
    :align: center

    A representation of Ethanol molecule. The proton subscripts correspond to the site
    indexes used in the spin system.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> H_CH3 = Site(isotope='1H', isotropic_chemical_shift=1.226)  # methyl proton
    >>> H_CH2 = Site(isotope='1H', isotropic_chemical_shift=2.61)  # methylene proton
    >>> H_OH = Site(isotope='1H', isotropic_chemical_shift=3.687)  # hydroxyl proton
    ...
    >>> etoh_sites = [H_CH3, H_CH3, H_CH3, H_CH2, H_CH2, H_OH]

Couplings
^^^^^^^^^

Now, we need to define the :math:`^3J_{HH}` couplings that cause the splittings
we're used to seeing in the spectrum of ethanol. In ``mrsimulator``, all Couplings
are defined using the :class:`~mrsimulator.Coupling` class. Let's start by defining
a coupling between a methyl and a methylene proton.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> HH_coupling_1 = Coupling(site_index=[0, 3], isotropic_j=7)
    >>> HH_coupling_1.property_units
    {'isotropic_j': 'Hz'}

The attribute *site_index* holds a pair of integers, where each integer is the index
of the coupled site object. The attribute *isotropic_j* is the isotropic *J*-coupling
between the coupled sites in units of *Hz*. Like every other object, the information on
the default unit is held with the ``property_units``  attribute.
In the above example, we define a coupling between site 0 (methyl) and site 3 (methylene).
The indexes 0 and 3 are relative to the list of site objects in ``etoh_sites``. The
isotropic *J*-coupling is 7 Hz.
Now, we define the rest of the methyl-methylene couplings and make a list to hold them all.

.. note::
    ``mrsimulator`` library does not support strong couplings, so we will be neglecting them.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> HH_coupling_2 = Coupling(site_index=[0, 4], isotropic_j=7)
    >>> HH_coupling_3 = Coupling(site_index=[1, 3], isotropic_j=7)
    >>> HH_coupling_4 = Coupling(site_index=[1, 4], isotropic_j=7)
    >>> HH_coupling_5 = Coupling(site_index=[2, 3], isotropic_j=7)
    >>> HH_coupling_6 = Coupling(site_index=[2, 4], isotropic_j=7)
    >>>
    >>> etoh_couplings = [
    ...     HH_coupling_1,
    ...     HH_coupling_2,
    ...     HH_coupling_3,
    ...     HH_coupling_4,
    ...     HH_coupling_5,
    ...     HH_coupling_6,
    ... ]


Spin system
^^^^^^^^^^^

Now, we add the sites and couplings to the spin system object.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> etoh = SpinSystem(sites=etoh_sites, couplings=etoh_couplings)

We have successfully built our ethanol spin system! If you need to
create more spin systems, repeat these instructions, but for this
example, we will stick with a single spin system.

Setting up the Method objects
-----------------------------
Next, we create a method to simulate a simple 1D pulse-acquire
:math:`^1H` spectrum.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> method_H = BlochDecaySpectrum(
    ...     channels=['1H'],
    ...     magnetic_flux_density=9.4,  # T
    ...     spectral_dimensions=[{
    ...         "count": 3000,
    ...         "spectral_width": 1.5e3,  # in Hz
    ...         "reference_offset": 940,  # in Hz
    ...         "label": "$^{1}$H frequency",
    ...     }],
    ... )


In the above code, *channels* is a list of isotope symbols that a method
will use. The Bloch Decay method only uses one channel, and in this case
we are simulating a :math:`^1\text{H}` spectrum. *magnetic_flux_density*
describes the environment under which the resonance frequency is
evaluated. *spectral_dimensions* contains a list of spectral dimensions
(only one for the Bloch Decay method). In this case, we define a frequency
dimension with 3000 points, spanning 1.5 kHz with a reference offset of 940 Hz.

You can create as many methods as you need, but in this case we will
stick with the one method.

Running simulation
------------------
Next, we need to create an instance of the simulator object and then
add our spin system and method to it. Then, we run the simulator with
the :meth:`~mrsimulator.Simulator.run` method.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> sim = Simulator()
    >>> sim.spin_systems = [etoh]
    >>> sim.methods = [method_H]
    >>> sim.run()

The simulator object has now processed the method with our spin system
and has stored the result in the simulation attribute of that method.
Let’s get the data from the method so we can plot it.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> H_data = sim.methods[0].simulation

Visualizing the dataset
-----------------------
Now that we have our data, let’s plot the spectrum using matplotlib!

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> plt.figure(figsize=(10, 4)) # set the figure size  # doctest: +SKIP
    >>> ax = plt.subplot(projection='csdm')  # doctest: +SKIP
    >>> ax.plot(H_data.real, color="black", linewidth=0.5)  # doctest: +SKIP
    >>> ax.set_xlim(4, 0.75)  # doctest: +SKIP
    >>> plt.tight_layout()  # doctest: +SKIP
    >>> plt.show()  # doctest: +SKIP

.. _fig-etoh-getting-started-coupled:
.. figure:: _static/null.*
    :alt: _images/null.png

    An example :math:`^{1}\text{H}` NMR spectrum simulation of Ethanol.
