
.. _getting_started:

==============================================
Getting started with `Mrsimulator`: The basics
==============================================

We have put together a set of guidelines for using the methods and
objects of the `Mrsimulator` package. We encourage our users
to follow these guidelines to promote consistency.
In `Mrsimulator`, the solid-state nuclear magnetic resonance (ssNMR)
lineshape is calculated through an instance of the :ref:`simulator_api`
class.

Import the :ref:`simulator_api` class using

.. doctest::

    >>> from mrsimulator import Simulator

and create an instance as follows,

.. doctest::

    >>> sim = Simulator()

Here, the variable ``sim`` is an instance of the
:ref:`simulator_api` class. The two attributes of this class that you will
frequently use are the :attr:`~mrsimulator.Simulator.isotopomers` and
:attr:`~mrsimulator.Simulator.methods`, whose values are a list of
:ref:`isotopomer_api` and :ref:`method_api` objects,
respectively. The default value of these attributes is an empty list.

.. doctest::

    >>> sim.isotopomers
    []
    >>> sim.methods
    []


Before you can start simulating
NMR lineshapes, you need to understand the role of Isotopomer and
Method objects. The following provides a brief description of the respective
objects.

.. For more information, we recommend reading :ref:`dictionary_objects`
.. and :ref:`dimension`.


Setting up the Isotopomer objects
---------------------------------
We consider an isotopomer as an isolated spin-system, which may contain multiple sites
and couplings between these sites. You may construct an isotopomer (spin-system) with
as many sites and couplings as necessary, for this example, we stick to a single site
spin-system, that is, an isotopomer with a single site. Let’s start by first building
a site.

A site object is a collection of parameters that describe the single site interactions.
In NMR, these spin interactions are described by a second-rank tensor. The single-site
interactions include the interaction between the magnetic dipole moment of the nucleus
and the surrounding magnetic field and the interaction between the electric quadrupole
moment of the nucleus with the surrounding electric field gradient. The later is zero
for sites with the spin quantum number, :math:`I=1/2`.

Let's start with a spin-1/2 isotope, :math:`^{29}\text{Si}`, and create
a site.

.. code-block:: python

    the_site = {
        "isotope": "29Si",
        "isotropic_chemical_shift": "-101.1 ppm",
        "shielding_symmetric": {"zeta": "70.5 ppm", "eta": 0.5},
    }

.. testsetup::
    >>> the_site = {
    ...     "isotope": "29Si",
    ...     "isotropic_chemical_shift": "-101.1 ppm",
    ...     "shielding_symmetric": {
    ...         "zeta": "70.5 ppm",
    ...         "eta": 0.5
    ...     }
    ... }

In the above code, ``the_site`` is a simplified python dictionary
representation of a :ref:`site_api` object. This site describes a
:math:`^{29}\text{Si}` isotope with a -101.1 ppm isotropic chemical shift
along with the symmetric part of the nuclear shielding anisotropy,
described here with parameters `zeta` and `eta` using Haeberlen convention.

That's it! Now that we have the site, we can create an isotopomer as,

.. code-block:: python

    the_isotopomer = {
        "name": "site A",
        "sites": [the_site],  # from the above code
        "abundance": "80%",
    }

.. testsetup::
    >>> the_isotopomer = {"name": "site A", "sites": [ the_site ],
    ...     "abundance": "80%"}

As mentioned before, an isotopomer is a collection of sites and couplings. In the above
example, we have an isotopomer with a single site and no couplings. Here, the attribute
`sites` hold a list of sites. The attributes `name` and `abundance` are optional.

..  .. seealso:: :ref:`dictionary_objects`, :ref:`isotopomer` and :ref:`site`.

Until now, we have only created a python dictionary representation of an isotopomer. To
run the simulation, you need to create an instance of the Isotopomer class. Use the
:meth:`~mrsimulator.Isotopomer.parse_dict_with_units` method of the Isotopomer class to
create the Isotopomer object, as follows,

    >>> from mrsimulator import Isotopomer
    >>> isotopomer_object = Isotopomer.parse_dict_with_units(the_isotopomer)

.. note:: We provide the :meth:`~mrsimulator.Isotopomer.parse_dict_with_units` method
    because it allows the user to add isotopomers, where the attribute value is a
    physical quantity, represented as a string with a value and a unit.
    Physical quantities remove the ambiguity in the units, which is otherwise
    a source of common confusion within many scientific applications. With this said,
    parsing physical quantities can add significant overhead when used in an iterative
    algorithm, such as the least-squares minimization. In such cases, we recommend
    defining objects directly. See the next topic for details.

We have successfully created one isotopomer. To create more isotopomers, repeat the
above set of instructions. In this example, we stick with a single isotopomer.
Once all the isotopomers ready, add the isotopomer objects to the instance of the
Simulator class, as follows

    >>> sim.isotopomers += [isotopomer_object] # add all isotopomers.


Setting up the Method objects
-----------------------------

A :ref:`method_api` object is a collection of parameters that describe an NMR method.
In Mrsimulator, all methods are described through five keywords -

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
      described with the keywords, `count` (:math:`N`), `spectral_width`
      (:math:`\nu_\text{sw}`), and `reference_offset` (:math:`\nu_0`). The
      coordinates are given as,

      .. math::
        \left([0 ... N-1] - \frac{T}{2}\right) \frac{\nu_\text{sw}}{N} + \nu_0

      where :math:`T=N` when :math:`N` is even else :math:`T=N-1`.

Let's start with the simplest method, the :func:`~mrsimulator.methods.BlochDecaySpectrum`.
The following is a python dictionary representation of the BlochDecaySpectrum method.

.. code-block::  python

    method_dict = {
        "channels": ["29Si"],
        "magnetic_flux_density": "9.4 T",
        "rotor_angle": "54.735 deg",
        "rotor_frequency": "0 Hz",
        "spectral_dimensions": [{
            "count": 2048,
            "spectral_width": "25 kHz",
            "reference_offset": "-8 kHz",
        }]
    }

.. testsetup::
    >>> method_dict = {
    ...     "channels": ["29Si"],
    ...     "magnetic_flux_density": "9.4 T",
    ...     "rotor_angle": "54.735 deg",
    ...     "rotor_frequency": "0 Hz",
    ...     "spectral_dimensions": [{
    ...         "count": 2048,
    ...         "spectral_width": "25 kHz",
    ...         "reference_offset": "-8 kHz",
    ...     }]
    ... }

Here, the key `channels` is a list of isotope symbols over which the method is applied. A
Bloch Decay method only has a single channel. In this example, it is given a value
of ``29Si``, which implies that the simulated lineshape from this method will comprise
frequency components arising from :math:`^{29}\text{Si}` resonances.
The keys `magnetic_flux_density`, `rotor_angle`, and `rotor_frequency` collectively
describe the spin-environment under which the resonance frequency is evaluated.
The key `spectral_dimensions` is a list of spectral dimensions. A Bloch Decay method only
has one spectral dimension.

Like before, you may parse the above ``method_dict`` to a create an instance using the
:meth:`~mrsimulator.methods.BlochDecaySpectrum.parse_dict_with_units` method as,

.. doctest::

    >>> from mrsimulator.methods import BlochDecaySpectrum
    >>> method_object = BlochDecaySpectrum.parse_dict_with_units(method_dict)

In the above example, the variable ``method_object`` is an instance of the
:class:`~mrsimulator.Method` class.

Likewise, you may create multiple method objects. In this example, we
stick with a single method. Finally, add all the method objects, in this case,
``method_object``, to the instance of the Simulator class, ``sim``, as follows,

.. doctest::

    >>> sim.methods += [method_object] # add all methods.

Running simulation
------------------

To simulate the line-shape, run the simulator with the
:meth:`~mrsimulator.Simulator.run` method, as follows,

.. note:: In Mrsimulator, the resonant frequencies are calculated assuming the
    weakly-coupled (Zeeman) basis for the spin-system.

.. doctest::

    >>> sim.run()

The simulator object, ``sim``, will process every method within the object over all
isotopomers and store the result in the :attr:`~mrsimulator.Method.simulation`
attribute of the respective Method object. In this example, we have a single method.
We can access the simulation data for this method as,

.. doctest::

    >>> data_1 = sim.methods[0].simulation
    >>> # data_n = sim.method[n].simulation # when there are multiple methods.

Here, ``data_1`` is a CSDM object holding the simulation data from the method
at index 0 of the :attr:`~mrsimulator.Simulator.methods` attribute from the ``sim``
object.

Visualizing the dataset
-----------------------

At this point, you may continue with additional post-simulation processing.
We end this example with a plot of the data from the simulation.
The following is the simulated lineshape plotted using the matplotlib
library.

.. doctest::

    >>> import matplotlib.pyplot as plt
    >>> def plot(x, y):
    ...     plt.figure(figsize=(4,3))
    ...     plt.plot(x,y)
    ...     plt.xlim([x.value.max(), x.value.min()]) # for reverse axis
    ...     plt.xlabel(f'frequency ratio / {str(x.unit)}')
    ...     plt.tight_layout()
    ...     plt.show()

    >>> plot(*data_1.to_list())  # doctest:+SKIP

.. .. testsetup::
..    >>> plot_save(freq, amp, "example")  # doctest: +SKIP

.. _fig2:
.. figure:: _images/example.*
    :figclass: figure

    An example static solid state NMR lineshape simulation.
