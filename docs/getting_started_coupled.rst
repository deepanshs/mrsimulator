
.. _getting_started_coupled:

===============================
The Basics: Coupled Spin System
===============================

In this next example, we demonstrate the simulation of a coupled spin system. We assume that
you are familiar with the setup of an uncoupled spin system. If not, please first read the
:ref:`getting_started` section before continuing.

Import the :ref:`simulator_api` class and create an instance as follows,

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> from mrsimulator import Simulator
    >>> sim = Simulator()

Here, the variable ``sim`` is an instance of the :ref:`simulator_api` class.

Setting up coupled SpinSystem objects
-------------------------------------
An NMR spin system is an isolated system of sites (spins) and couplings. You may
construct a spin system with as many sites and couplings as necessary; for this
example, we stick to a two-site coupled spin system. Let’s start by first building
the sites.

Let's start with spin-1/2, :math:`^{29}\text{Si}`, and spin-5/2, :math:`^{27}\text{Al}`,
isotopes, and create sites.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> sites = [
    ...   {
    ...     "isotope": "29Si",
    ...     "isotropic_chemical_shift": "-101.1 ppm",
    ...     "shielding_symmetric": {"zeta": "70.5 ppm", "eta": 0.5},
    ...   },
    ...   {
    ...     "isotope": "27Al",
    ...     "isotropic_chemical_shift": "50.1 ppm",
    ...     "quadrupolar": {"Cq": "5.1 MHz", "eta": 0.4},
    ...   }
    ... ]

In the above code, the variable ``sites`` is a list of two sites. The site at index 0
describes a :math:`^{29}\text{Si}` isotope with a -101.1 ppm isotropic chemical shift along
with the symmetric part of the nuclear shielding tensor, described here with the parameters
`zeta` and `eta` using the Haeberlen convention. The site at index 1 describes a
:math:`^{27}\text{Al}` isotope with a 50.1 ppm isotropic chemical shift along with the
symmetric quadrupolar tensor, described here with the parameters `Cq` and `eta`.

That's it for the sites! Now let's create the coupling between the two sites.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> coupling_AX = {
    ...     "site_index": [0, 1],
    ...     "isotropic_j": "200 Hz",
    ... }

In the above code, ``coupling_AX`` is a simplified python dictionary representation of the
:ref:`coupling_api` object. Here, the coupling is between the sites at index 0 and 1 from the
previously defined ``sites`` list. In addition, we also provide a 200 Hz isotropic `J`-coupling
between the sites with the key `isotropic_j`.

Now that we have the sites and coupling, we can create a two-site coupled spin system as follows,

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> coupled_spin_system = {
    ...     "name": "AX system",
    ...     "description": "A test weakly coupled AX spin system",
    ...     "sites": sites,              # from the above code
    ...     "couplings": [coupling_AX],   # from the above code
    ...     "abundance": "100%",
    ... }

As mentioned before, a spin system is a collection of sites and couplings. In the above
example, we have created a coupled spin system with two sites and one coupling. Here,
the attribute `sites` hold a list of two sites. The attributes `name`, `description`, and
`abundance` are optional.


Similar to the previous example, import the SpinSystem class and use its
:meth:`~mrsimulator.SpinSystem.parse_dict_with_units` method to parse the python
dictionary and create an instance of the spin system class, as follows,

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> from mrsimulator import SpinSystem
    >>> system_object_1 = SpinSystem.parse_dict_with_units(coupled_spin_system)

.. note:: We provide the :meth:`~mrsimulator.SpinSystem.parse_dict_with_units` method
    because it allows the user to create spin systems, where the attribute value is a
    physical quantity, represented as a string with a value and a unit.
    Physical quantities remove the ambiguity in the units, which is otherwise
    a source of general confusion within many scientific applications. With this said,
    parsing physical quantities can add significant overhead when used in an iterative
    algorithm, such as the least-squares minimization. In such cases, we recommend
    defining objects directly. See the :ref:`using_objects` for details.

We have successfully created a coupled spin system object. To create more spin system objects,
repeat the above set of instructions. In this example, we stick with a single
spin system object. Once all spin system objects are ready, add these objects to the
instance of the Simulator class, as follows

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> sim.spin_systems += [system_object_1] # add all spin system objects.


Setting up the Method objects
-----------------------------

Let's start with the :func:`~mrsimulator.methods.BlochDecaySpectrum` method.
The following is a python dictionary representation of the BlochDecaySpectrum method.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> method_dict_Si = {
    ...     "channels": ["29Si"],
    ...     "magnetic_flux_density": "9.4 T",
    ...     "rotor_angle": "54.735 deg",
    ...     "rotor_frequency": "3 kHz",
    ...     "spectral_dimensions": [{
    ...         "count": 16384,
    ...         "spectral_width": "25 kHz",
    ...         "reference_offset": "-8 kHz",
    ...         "label": r"$^{29}$Si resonances",
    ...     }]
    ... }

Here, the key `channels` is a list of isotope symbols over which the method is applied.
A Bloch Decay method only has a single channel. In this example, it is given a value
of ``29Si``, which implies that the simulated spectrum from this method will comprise
frequency components arising from the :math:`^{29}\text{Si}` resonances.
The keys `magnetic_flux_density`, `rotor_angle`, and `rotor_frequency` collectively
describe the spin environment under which the resonance frequency is evaluated.
The key `spectral_dimensions` is a list of spectral dimensions. A Bloch Decay method
only has one spectral dimension. In this example, the spectral dimension defines a
frequency dimension with 2048 points, spanning 25 kHz with a reference offset of
-8 kHz.

Let's create another method, :func:`~mrsimulator.methods.BlochDecayCTSpectrum` method,
which is a central transition selective Bloch decay spectrum method. The method is defined similarly
to the Bloch decay spectrum method. The following is a python dictionary representation of the method.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> method_dict_Al = {
    ...     "channels": ["27Al"],
    ...     "magnetic_flux_density": "9.4 T",
    ...     "rotor_angle": "54.735 deg",
    ...     "rotor_frequency": "8 kHz",
    ...     "spectral_dimensions": [{
    ...         "count": 2048,
    ...         "spectral_width": "25 kHz",
    ...         "reference_offset": "2 kHz",
    ...         "label": r"$^{27}$Al resonances",
    ...     }]
    ... }

Like before, you may parse the above two method dictionaries using the
:meth:`~mrsimulator.methods.BlochDecaySpectrum.parse_dict_with_units` function of the respective
methods. Import the two method class and create an instance of the methods as follows

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> from mrsimulator.methods import BlochDecaySpectrum, BlochDecayCTSpectrum
    >>> method_object_Si = BlochDecaySpectrum.parse_dict_with_units(method_dict_Si)
    >>> method_object_Al = BlochDecayCTSpectrum.parse_dict_with_units(method_dict_Al)

Here, ``method_object_Si`` and ``method_object_Al`` are the instances of the
:class:`~mrsimulator.Method` class.

Likewise, you may create additional method objects. In this example, we
stick with the two methods. Finally, add all the method objects to the instance of the
Simulator class, ``sim``, as follows,

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> sim.methods += [method_object_Si, method_object_Al] # add all methods.

Running simulation
------------------

To simulate the spectrum, run the simulator with the :meth:`~mrsimulator.Simulator.run`
method, as follows,

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> sim.run()

.. note:: In ``mrsimulator``, all resonant frequencies are calculated assuming the
    weakly-coupled (Zeeman) basis for the spin system.

The simulator object, ``sim``, will process every method over all the spin systems and
store the result in the :attr:`~mrsimulator.Method.simulation` attribute of the
respective Method object. In this example, we have two methods. You may access
the simulation data for these methods as,

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> data_0 = sim.methods[0].simulation
    >>> data_1 = sim.methods[1].simulation
    >>> # data_n = sim.method[n].simulation # when there are multiple methods.

Here, ``data_0`` is a CSDM object holding the simulation data from the method
at index 0 of the :attr:`~mrsimulator.Simulator.methods` attribute from the ``sim``
object, that is, the :math:`^{29}\text{Si}` spectrum. The ``data_1`` is a CSDM object
corresponding to the method at index 1, that is, the :math:`^{27}\text{Al}` spectrum.

.. seealso::
    **CSDM:** The core scientific dataset model (CSDM) [#f1]_ is a lightweight and portable
    file format model for multi-dimensional scientific datasets and is supported by numerous
    NMR software---DMFIT, SIMPSON, jsNMR, and RMN. We also provide a python package
    `csdmpy <https://csdmpy.readthedocs.io/en/stable/>`_.

Visualizing the dataset
-----------------------

At this point, you may continue with additional post-simulation processing.
We end this example with a plot of the data from the simulation. :numref:`fig1-getting-started-coupled`
and :numref:`fig2-getting-started-coupled` depicts the plot of :math:`^{29}\text{Si}` and
:math:`^{27}\text{Al}` NMR simulated spectrum.

For a quick plot of the csdm data, you may use the `csdmpy <https://csdmpy.readthedocs.io/en/stable/>`_
library. The `csdmpy` package uses the matplotlib library to produce basic plots.
You may optionally customize the plot using matplotlib methods.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> import matplotlib.pyplot as plt
    >>> plt.figure(figsize=(6, 3.5)) # set the figure size # doctest: +SKIP
    >>> ax = plt.subplot(projection='csdm') # doctest: +SKIP
    >>> ax.plot(data_0, linewidth=1.5) # doctest: +SKIP
    >>> ax.invert_xaxis() # reverse x-axis # doctest: +SKIP
    >>> plt.tight_layout(pad=0.1) # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

.. _fig1-getting-started-coupled:
.. figure:: _static/null.*
    :alt: _images/null.png

    An example :math:`^{29}\text{Si}` NMR simulation from a coupled Si-Al spin system.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> plt.figure(figsize=(6, 3.5)) # set the figure size # doctest: +SKIP
    >>> ax = plt.subplot(projection='csdm') # doctest: +SKIP
    >>> ax.plot(data_1, linewidth=1.5) # doctest: +SKIP
    >>> ax.invert_xaxis() # reverse x-axis # doctest: +SKIP
    >>> plt.tight_layout(pad=0.1) # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

.. _fig2-getting-started-coupled:
.. figure:: _static/null.*
    :alt: _images/null.png

    An example :math:`^{27}\text{Al}` NMR simulation from a coupled Si-Al spin system.


.. [#f1] Srivastava, D. J., Vosegaard, T., Massiot, D., Grandinetti, P. J.
        Core Scientific Dataset Model: A lightweight and portable model and file format
        for multi-dimensional scientific data. PLOS ONE, 2020, **15**, 1.
        `DOI 10.1371/e0225953 <https://doi.org/10.1371/journal.pone.0225953>`_
