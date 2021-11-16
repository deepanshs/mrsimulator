
.. _getting_started:

==========
The Basics
==========

We have put together a set of guidelines for using the ``mrsimulator`` package. We
encourage our users to follow these guidelines for consistency. In
``mrsimulator``, the solid-state nuclear magnetic resonance (ssNMR) spectrum is
calculated through an instance of the :ref:`simulator_api` class.

Import the :ref:`simulator_api` class using

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> from mrsimulator import Simulator

and create an instance as follows,

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> sim = Simulator()

Here, the variable ``sim`` is an instance of the :ref:`simulator_api` class. The two
attributes of this class that you will frequently use are the
:attr:`~mrsimulator.Simulator.spin_systems` and
:attr:`~mrsimulator.Simulator.methods`, whose values are a list of
:ref:`spin_sys_api` and :ref:`method_api` objects,
respectively. The default value of these attributes is an empty list.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> sim.spin_systems
    []
    >>> sim.methods
    []


Before you can start simulating the NMR spectrum, you need to understand the role of
the SpinSystem and Method objects. The following provides a brief description of the
respective objects.

.. For more information, we recommend reading :ref:`dictionary_objects`
.. and :ref:`dimension`.


Setting up the SpinSystem objects
---------------------------------
An NMR spin system is an isolated system of sites (spins) and couplings. You may
construct a spin system with as many sites and couplings as necessary; for this
example, we stick to a single-site spin system. Let’s start by first building
a site.

A site object is a collection of attributes that describe site-specific interactions.
In NMR, these spin interactions are described by a second-rank tensor.
Site-specific interactions include the interaction between the magnetic dipole moment
of the nucleus and the surrounding magnetic field and the interaction between the
electric quadrupole moment of the nucleus with the surrounding electric field gradient.
The latter is zero for sites with the spin quantum number, :math:`I=1/2`.

Let's start with a spin-1/2 isotope, :math:`^{29}\text{Si}`, and create a site.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> from mrsimulator import Site
    >>> from mrsimulator.spin_system.tensors import SymmetricTensor
    >>> the_site = Site(
    ...     isotope="29Si",
    ...     isotropic_chemical_shift=-101.1,  # in ppm
    ...     shielding_symmetric=SymmetricTensor(
    ...         zeta=70.5,  # in ppm
    ...         eta=0.5,
    ...     )
    ... )

This :ref:`site_api` describes a :math:`^{29}\text{Si}` isotope with a
-101.1 ppm isotropic chemical shift along with the symmetric part of the nuclear
shielding anisotropy tensor, described here with the parameters *zeta* and *eta* using
the Haeberlen convention.

.. note::

    We strongly recommend directly instantiating objects in your code since it improves
    readability. However, all objects in ``mrsimulator`` can be represented by and
    constructed from a dictionary.

    Looking at the :ref:`site_api` object we can call the :py:meth:`~mrsimulator.Site.json`
    method to get the JSON-style dictionary. Note that values are serialized with units.

    .. plot::
        :format: doctest
        :context: close-figs
        :include-source:

        >>> import pprint
        >>> py_dict = the_site.json()
        >>> pprint.pprint(py_dict)
        {'isotope': '29Si',
         'isotropic_chemical_shift': '-101.1 ppm',
         'shielding_symmetric': {'eta': 0.5, 'zeta': '70.5 ppm'}}

    Likewise a site object can be constructed from a dictionary by calling the
    :py:meth:`mrsimulator.Site.parse_dict_with_units` method.

    .. plot::
        :format: doctest
        :context: close-figs
        :include-source:

        >>> same_site = Site.parse_dict_with_units(py_dict)
        >>> the_site == same_site
        True

    We provide the ``parse_dict_with_units`` method
    because it allows the user to create an object where each attribute value is a
    physical quantity represented by a string with a value and a unit.
    Physical quantities remove the ambiguity in the units, which is otherwise
    a source of general confusion within many scientific applications. With this said,
    parsing physical quantities will add significant overhead when used in an iterative
    algorithm, such as the least-squares minimization.

That's it! Now that we have a site, we can create a single-site spin system following,

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> from mrsimulator import SpinSystem
    >>> system_object_1 = SpinSystem(
    ...     name="site A",
    ...     description="A test 29Si site",
    ...     sites=[the_site],  # from the above code
    ...     abundance=80,  # percentage
    ... }

As mentioned before, a spin system is a collection of sites and couplings. In the above
example, we have created a spin system with a single site and no coupling. Here, the
attribute *sites* hold a list of sites. The *abundance* attribute is an optional percentage
describing the abundance of the spin system and defaults to 100%. The attributes *name* and
*description* are optional and have no effect on the resulting spectrum.

..  .. seealso:: :ref:`dictionary_objects`, :ref:`spin_system` and :ref:`site`.

..  Until now, we have only created a python dictionary representation of a spin system. To
..  run the simulation, you need to create an instance of the
..  :py:class:`~mrsimulator.SpinSystem` class. Import the SpinSystem class and use it's
..  :py:meth:`~mrsimulator.SpinSystem.parse_dict_with_units` method to parse the python
..  dictionary and create an instance of the spin system class, as follows,
..
..  .. plot::
..      :format: doctest
..      :context: close-figs
..      :include-source:
..
..      >>> from mrsimulator import SpinSystem
..      >>> system_object_1 = SpinSystem.parse_dict_with_units(the_spin_system)

We have successfully created a spin system object. To create more spin system objects,
repeat the above set of instructions. In this example, we stick with a single
spin system object. Once all spin system objects are ready, add these objects to the
instance of the Simulator class, as follows

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> sim.spin_systems = [system_object_1]


Setting up the Method objects
-----------------------------

A :ref:`method_api` object is a collection of attributes that describe an NMR method.
In ``mrsimulator``, all methods are described through five keywords -

.. cssclass:: table-bordered

.. list-table::
  :widths: 25 75
  :header-rows: 1

  * - Keywords
    - Description
  * - channels
    - A list of isotope symbols over which the given method applies.
  * - magnetic_flux_density
    - The macroscopic magnetic flux density of the applied external magnetic field.
  * - rotor_angle
    - The angle between the sample rotation axis and the applied external magnetic field.
  * - rotor_frequency
    - The sample rotation frequency.
  * - spectral_dimensions
    - A list of spectral dimensions. The coordinates along each spectral dimension is
      described with the keywords, *count* (:math:`N`), *spectral_width*
      (:math:`\nu_\text{sw}`), and *reference_offset* (:math:`\nu_0`). The
      coordinates are evaluated as,

      .. math::
        \left([0, 1, 2, ... N-1] - \frac{T}{2}\right) \frac{\nu_\text{sw}}{N} + \nu_0

      where :math:`T=N` when :math:`N` is even else :math:`T=N-1`.

Let's start with the simplest method, the :func:`~mrsimulator.methods.BlochDecaySpectrum`.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> from mrsimulator.methods import BlochDecaySpectrum
    >>> from mrsimulator.method.spectral_dimensions import SpectralDimension

    >>> the_method = BlochDecaySpectrum(
    ...     channels=["29Si"],
    ...     magnetic_flux_density=9.4,  # in T
    ...     rotor_angle=0.9553166,  # in rad (magic angle)
    ...     rotor_frequency=0,
    ...     spectral_dimensions=[SpectralDimension(
    ...         count=2048,
    ...         spectral_width=25e3,    # in Hz
    ...         reference_offset=-8e3,  # in Hz
    ...         label=r"$^{29}$Si resonances",
    ...     )]
    ... )

Here, the key *channels* is a list of isotope symbols over which the method is applied.
A Bloch Decay method only has a single channel. In this example, it is given a value
of ``29Si``, which implies that the simulated spectrum from this method will comprise
frequency components arising from the :math:`^{29}\text{Si}` resonances.
The keys *magnetic_flux_density*, *rotor_angle*, and *rotor_frequency* collectively
describe the spin environment under which the resonance frequency is evaluated.
The key *spectral_dimensions* is a list of spectral dimensions. A Bloch Decay method
only has one spectral dimension. In this example, the spectral dimension defines a
frequency dimension with 2048 points, spanning 25 kHz with a reference offset of
-8 kHz.

..  Like before, you may parse the above ``method_dict`` using the
..  :py:meth:`~mrsimulator.methods.BlochDecaySpectrum.parse_dict_with_units` function of the
..  method. Import the BlochDecaySpectrum class and create an instance of the method,
..  following,
..
..  .. plot::
..      :format: doctest
..      :context: close-figs
..      :include-source:
..
..      >>> from mrsimulator.methods import BlochDecaySpectrum
..      >>> method_object = BlochDecaySpectrum.parse_dict_with_units(method_dict)
..
..  Here, ``the_method`` is an instance of the :py:class:`~mrsimulator.Method` class.

A simulator can hold multiple method objects like in our link(¹³C MAS NMR of Glycine (CSA) multi-spectra fit).
In this example, we stick with a single method. Finally, add all the method objects,
in this case, ``the_method``, to the instance of the Simulator class, ``sim``, as follows,

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> sim.methods = [the_method]

Running simulation
------------------

To simulate the spectrum, run the simulator with the :py:meth:`~mrsimulator.Simulator.run`
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
respective Method object. In this example, we have a single method. You may access
the simulation data for this method as,

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> data_0 = sim.methods[0].simulation
    >>> # data_n = sim.method[n].simulation # when there are multiple methods.

Here, ``data_0`` is a CSDM object holding the simulation data from the method
at index 0 of the :attr:`~mrsimulator.Simulator.methods` attribute from the ``sim``
object.

.. seealso::
    **CSDM:** The core scientific dataset model (CSDM) [#f1]_ is a lightweight and portable
    file format model for multi-dimensional scientific datasets and is supported by numerous
    NMR software---DMFIT, SIMPSON, jsNMR, and RMN. We also provide a python package
    `csdmpy <https://csdmpy.readthedocs.io/en/stable/>`_.

Visualizing the dataset
-----------------------

At this point, you may continue with additional post-simulation processing.
We end this example with a plot of the data from the simulation.
:numref:`fig1-getting-started` depicts the plot of the simulated spectrum.

For a quick plot of the csdm data, you may use the `csdmpy <https://csdmpy.readthedocs.io/en/stable/>`_
library. The *csdmpy* package uses the matplotlib library to produce basic plots.
You may optionally customize the plot using matplotlib methods.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> import matplotlib.pyplot as plt
    >>> plt.figure(figsize=(6, 3.5)) # set the figure size  # doctest: +SKIP
    >>> ax = plt.subplot(projection='csdm')  # doctest: +SKIP
    >>> ax.plot(data_0, linewidth=1.5)  # doctest: +SKIP
    >>> ax.invert_xaxis() # reverse x-axis  # doctest: +SKIP
    >>> plt.tight_layout(pad=0.1)  # doctest: +SKIP
    >>> plt.show()  # doctest: +SKIP

.. _fig1-getting-started:
.. figure:: _static/null.*
    :alt: _images/null.png

    An example of solid-state static NMR spectrum simulation.


.. **For advanced users**

.. Advanced uses may prefer to apply some more sophisticated processing or use some other
.. plotting libraries. For those users, you may extract the data from the csdm object
.. as a list of arrays using the `to_list() <https://csdmpy.readthedocs.io/en/stable/api/CSDM.html#csdmpy.CSDM.to_list>`_
.. method of the csdm object, following,

.. .. plot::
..     :format: doctest
..     :context: close-figs
..     :include-source:

..     >>> x, y = data_0.to_list()

.. Here, ``x`` is a quantity array and contains the coordinates of the spectral dimension
.. in units of ppm, and ``y`` is the response array.

.. The following is a matplotlib script
.. which uses the above ``x``, and ``y`` variables to generate a similar plot shown in
.. :numref:`fig1-getting-started`.

.. .. doctest::

..     >>> import matplotlib.pyplot as plt
..     >>> def plot(x, y):
..     ...     plt.figure(figsize=(4,3))
..     ...     plt.plot(x,y)
..     ...     plt.xlim([x.value.max(), x.value.min()]) # for reverse axis
..     ...     plt.xlabel(f'frequency ratio / {str(x.unit)}')
..     ...     plt.tight_layout()
..     ...     plt.show()

..     >>> plot(x, y)  # doctest:+SKIP

.. .. testsetup::
..    >>> plot_save(freq, amp, "example")  # doctest: +SKIP

.. [#f1] Srivastava, D. J., Vosegaard, T., Massiot, D., Grandinetti, P. J.
        Core Scientific Dataset Model: A lightweight and portable model and file format
        for multi-dimensional scientific data. PLOS ONE, 2020, **15**, 1.
        `DOI 10.1371/e0225953 <https://doi.org/10.1371/journal.pone.0225953>`_
