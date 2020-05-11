
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
:attr:`~mrsimulator.Simulator.methods`, whose value is a list of
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
We consider an isotopomer as an isolated spin-system with multiple sites and couplings
between the sites.
You may construct an isotopomer (spin-system) with as many sites and coupling as
you may desire, for this example, we will stick to a single site spin-system, that
is, an isotopomer with a single site. Let's start by first building a site.

A site object is a collection of parameters that describe the single site interactions.
In NMR, these spin interactions are described by a second-rank tensor. The interactions
that involve a single site are the nuclear shielding interaction between the magnetic
dipole moment of the nucleus and the surrounding magnetic field, and the interaction
between the electric quadrupole moment of the nucleus with the surrounding electric
field gradient. The later is zero for sites with the spin quantum number, :math:`I=1/2`.

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
        "sites": [the_site],  # from previous code
        "abundance": "80%",
    }

.. testsetup::
    >>> the_isotopomer = {"name": "site A", "sites": [ the_site ],
    ...     "abundance": "80%"}

As mentioned before, an isotopomer is a collection of sites and couplings. In the above
example, we created an isotopomer with a single site and no coupling. Here, the attribute
`sites` hold a list of sites. The attributes `name` and `abundance` are optional.

..  .. seealso:: :ref:`dictionary_objects`, :ref:`isotopomer` and :ref:`site`.

Until now, we have only created a python dictionary representation of an isotopomer. To
run the simulation, you first need to create an instance of the Isotopomer class. Use the :meth:`~mrsimulator.Isotopomer.parse_dict_with_units` method of the Isotopomer class to
create the Isotopomer object, as follows,

    >>> from mrsimulator import Isotopomer
    >>> isotopomer_object = Isotopomer.parse_dict_with_units(the_isotopomer)

.. note:: We provide the :meth:`~mrsimulator.Isotopomer.parse_dict_with_units` method
    because it allows the user to add isotopomer, where the attribute value is a physical
    quantity. Physical quantities remove the ambiguity in the units, which is otherwise
    a source of common confusion within many scientific applications. With this said, parsing
    physical quantities can add significant overhead when used in an iterative algorithm,
    such as the least-squares minimization. In such cases, we recommend defining objects
    directly. See the next section.

You may create as many isotopomer objects as necessary, although in this
example, we stick with a single isotopomer.
Finally, add the isotopomer objects, in this case, the variable ``isotopomer_object``,
to the instance of the Simulator class, ``sim``, as follows,

    >>> sim.isotopomers += [isotopomer_object]


Setting up the Method objects
-----------------------------

A :ref:`method_api` object describes the NMR method.
Let's start with the simplest NMR method, :func:`~mrsimulator.methods.BlochDecaySpectrum`.

A Bloch Decay Spectrum method is a single event method. The following a python
dictionary representation of the BlochDecaySpectrum method.

.. code-block::  python

    method_dict = {
        "channels": ["29Si"],
        "magnetic_flux_density": "9.4 T",
        "rotor_angle": "54.735 deg",
        "rotor_frequency": "0 Hz",
        "dimensions": [{
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
    ...     "dimensions": [{
    ...         "count": 2048,
    ...         "spectral_width": "25 kHz",
    ...         "reference_offset": "-8 kHz",
    ...     }]
    ... }

Here, the key `channels` is a list of isotope symbols over which the method is applied. A
Bloch Decay method only has a single channel, which in this example, is given a value
of ``29Si``, which implies that the simulated lineshape from this method will comprise of
frequency components arising from :math:`^{29}\text{Si}` resonances.
The keys `magnetic_flux_density`, `rotor_angle`, and `rotor_frequency` collectively
describe the spin-environment under which the resonance frequency is evaluated.
The key `dimensions` is a list of spectral dimensions. A Bloch Decay method has
a single spectral dimension. Each spectral dimension is described by the keys,
`count`, `spectral_width`, and `reference_offset`, which collectively define the
grid coordinates along the spectral dimension,

.. math::

    \text{coordinates} = ([0 ... \text{count}] - N/2) \text{spectral_width}/\text{count} + \text{reference_offset}

where `N=count` when count is even else `N=count-1`

As before, you may parse the above ``method_dict`` to a create an instance the
method as,

..  .. seealso:: :ref:`dimension`.

.. doctest::

    >>> from mrsimulator.methods import BlochDecaySpectrum
    >>> method_object = BlochDecaySpectrum.parse_dict_with_units(method_dict)

In the above example, the variable ``method`` is an instance of the BlochDecaySpectrum
object.

You may create multiple method objects if required. In this example, we
stick with a single method. Finally, add the method, in this case, ``method_object``,
to the instance of the Simulator class, ``sim``, as follows,

.. doctest::

    >>> sim.methods += [method_object]

Running simulation
------------------

To simulate the line-shape, run the simulator with the
:meth:`~mrsimulator.Simulator.run` method, as follows,

.. note:: In Mrsimulator, the resonant frequencies are calculated assuming the
    weakly-coupled (Zeeman) basis for the spin-system.

.. doctest::

    >>> sim.run()

The simulator object, ``sim``, will process all the methods within the object
and store the result in the :attr:`~mrsimulator.Method.simulation`
attribute of the Method object. You may access the simulation as follows,

.. doctest::

    >>> data1 = sim.methods[0].simulation

Here, ``data1`` is the CSDM object holding the simulation data from the method
at index 0 of the :attr:`~mrsimulator.Simulator.methods` attribute. We use CSDM
format because of its versatility in handling multi-dimensional datasets.

At this point, you may continue with additional post-simulation processing.
We end this example with the plot of the data from the simulation. Because NMR data
is visualized on a dimensionless frequency ratio axis, we first convert the frequency
to a `ppm` scale, following,

.. doctest::

    >>> data1.dimensions[0].to('ppm', 'nmr_frequency_ratio')

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

    >>> plot(*data1.to_list())

.. .. testsetup::
..    >>> plot_save(freq, amp, "example")  # doctest: +SKIP

.. _fig2:
.. figure:: _images/example.*
    :figclass: figure-polaroid

    An example static solid state NMR lineshape simulation
