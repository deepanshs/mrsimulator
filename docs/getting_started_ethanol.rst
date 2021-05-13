.. _getting_started_coupled_spin_system_etoh:

===============================
The Basics: Coupled Spin System
===============================

In this example, we will simulate the :math:`^1H` NMR spectrum of
ethanol. Let’s start by importing all the necessary packages.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from mrsimulator import Simulator, SpinSystem, Site, Coupling
    >>> from mrsimulator.methods import BlochDecaySpectrum
    >>> from mrsimulator import signal_processing as sp

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

Next, we define :math:`^3J_{HH}` couplings and make a list to hold them
all. mrsimulator does not support strong couplings, so we will be
neglecting them.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> # all methyl-methylene coupling pairs
    >>> HH_coupling_1 = Coupling(site_index=[0, 3], isotropic_j=7)
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
    ...         "count": 16000,
    ...         "spectral_width": 1.5e3,  # in Hz
    ...         "reference_offset": 950,  # in Hz
    ...         "label": "$^{1}$H frequency",
    ...     }],
    ... )


In the above code, *channels* is a list of isotope symbols that a method
will use. The Bloch Decay method only uses one channel, and in this case
we are simulating a :math:`^1H` spectrum. *magnetic_flux_density*
describes the environment under which the resonance frequency is
evaluated. *spectral_dimensions* contains a list of spectral dimensions
(only one for the Bloch Decay method). In this case, we define a
frequency dimension with 16,000 points, spanning 1.5 kHz with a
reference offset of 950 Hz.

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
Now that we have our data, let’s add some post-simulation processing. We
define a SignalProcessor object that adds an exponential apodization of
1 Hz and then apply this processor on our data from the simulation.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> processor = sp.SignalProcessor(
    ...     operations=[
    ...         sp.IFFT(),
    ...         sp.apodization.Exponential(FWHM="1 Hz"),
    ...         sp.FFT(),
    ...     ]
    ... )
    >>> processed_H_data = processor.apply_operations(data=H_data)


**Plot**

Now, let’s plot the spectrum using matplotlib!

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> plt.figure(figsize=(6, 4)) # set the figure size  # doctest: +SKIP
    >>> ax = plt.subplot(projection='csdm')  # doctest: +SKIP
    >>> ax.plot(processed_H_data.real, color="black", linewidth=0.5)  # doctest: +SKIP
    >>> ax.invert_xaxis()  # doctest: +SKIP
    >>> plt.tight_layout()  # doctest: +SKIP
    >>> plt.show()  # doctest: +SKIP

.. _fig-etoh-getting-started-coupled:
.. figure:: _static/null.*
    :alt: _images/null.png

    An example :math:`^{1}\text{H}` NMR spectrum simulation of Ethanol.
