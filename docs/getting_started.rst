
.. _getting_started:

.. testsetup::

    >>> import matplotlib
    >>> font = {'family': 'Helvetica', 'weight': 'light', 'size': 9}
    >>> matplotlib.rc('font', **font)
    >>> from os import path

    >>> import matplotlib.pyplot as plt
    >>> def plot_save(x, y, filename):
    ...     plt.figure(figsize=(4, 3))
    ...     plt.plot(x, y, linewidth=1)
    ...     plt.xlim([x.value.max(), x.value.min()])
    ...     plt.xlabel(f"frequency ratio / {str(x.unit)}", **font)
    ...     plt.grid(color='gray', linestyle='--', linewidth=1.0, alpha=0.25)
    ...     plt.tight_layout(h_pad=0, w_pad=0, pad=0)
    ...
    ...     filename = path.split(filename)[1]
    ...     filepath = './docs/_images'
    ...     pth = path.join(filepath, filename)
    ...     plt.savefig(pth+'.pdf')
    ...     plt.savefig(pth+'.png', dpi=100)
    ...     plt.close()


==============================================
Getting started with `Mrsimulator`: The basics
==============================================

We have put together a set of guidelines for using various methods and
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
frequently use are :attr:`~mrsimulator.simulator.Simulator.isotopomers` and
:attr:`~mrsimulator.simulator.Simulator.dimensions`, whose value is a list of
:ref:`isotopomer_api` and :ref:`dimension_api` objects,
respectively. The default value of these attributes is an empty list.

.. doctest::

    >>> sim.isotopomers
    []
    >>> sim.dimensions
    []


Before you can start simulating
NMR lineshapes, you need to understand the role of Isotopomer and
Dimension objects. The following provides a brief description of the respective
objects.

.. For more information, we recommend reading :ref:`dictionary_objects`
.. and :ref:`dimension`.


Setting up Isotopomer objects
-----------------------------
For all practical purposes, an isotopomer may be described as an isolated
spin-system containing multiple sites and couplings between them. In the
current version, we focus on a single site spin-system, that is,
an isotopomer with a single site. Let's start by building a site.

In NMR, an active site may be described by a second-rank nuclear shielding
interaction tensor and additionally a second-rank electric quadrupole
interaction tensor for isotopes with the spin quantum number :math:`I>1/2`.
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
along with nuclear shielding anisotropy, described here with parameters `zeta`
and `eta` using Haeberlen convention.

Now that we have our site, we can create an isotopomer with this site, as
follows,

.. code-block:: python

    the_isotopomer = {
        "name": "site A",
        "sites": [the_site],  # from previous code
        "abundance": "80%",
    }

.. testsetup::
    >>> the_isotopomer = {"name": "site A", "sites": [ the_site ],
    ...     "abundance": "80%"}

The above isotopomer contains ``the_site`` as the value of the `sites`
attribute of the isotopomer. In addition to the site, we have also provided
an optional `name` and `abundance` to the isotopomer.

..  .. seealso:: :ref:`dictionary_objects`, :ref:`isotopomer` and :ref:`site`.

An instance of an Isotopomer class may be created from the above python
dictionary, using the :meth:`~mrsimulator.Isotopomer.parse_dict_with_units`
method of the Isotopomer class as follows,

    >>> from mrsimulator import Isotopomer
    >>> isotopomer_object = Isotopomer.parse_dict_with_units(the_isotopomer)

Here, the ``isotopomer_object`` is an instance of the Isotopomer class.
You may create as many isotopomer objects as necessary, although in this
example, we stick with a single isotopomer. Finally, add the isotopomer
objects, in this case, the variable ``isotopomer_object``, to the instance of
the Simulator class, ``sim``, as follows,

    >>> sim.isotopomers += [isotopomer_object]


Setting up Dimension objects
----------------------------

The :ref:`dimension_api` object describes a spectroscopic dimension of the
NMR spectrum. The number of dimension objects required in the simulation
depends on the dimensionality of the problem. For example, a one-pulse acquired
1-D spectrum requires a single dimension object, while two-dimensional spectrum
requires two dimension objects. In `Mrsimulator`, the Dimension object is
designed to include attributes required for evaluating the spectrum/line-shape
along that dimension.

.. code-block:: python

    dimension = {
        "isotope": "29Si",
        "magnetic_flux_density": "9.4 T",
        "rotor_angle": "54.735 deg",
        "rotor_frequency": "0 kHz",
        "number_of_points": 2048,
        "spectral_width": "25 kHz",
        "reference_offset": "-8 kHz",
    }

.. testsetup::
    >>> dimension = {
    ...     "isotope": "29Si",
    ...     "magnetic_flux_density": "9.4 T",
    ...     "rotor_angle": "54.735 deg",
    ...     "rotor_frequency": "0 kHz",
    ...     "number_of_points": 2048,
    ...     "spectral_width": "25 kHz",
    ...     "reference_offset": "-8 kHz"
    ... }

In the above example, the variable ``dimension`` holds a python dictionary
representation of a Dimension object. Here, the value of the
`isotope` key is the isotope symbol of the observed nucleus. A value, ``29Si``,
implies that the simulated lineshape will comprise of frequency components
arising from :math:`^{29}\text{Si}` resonances.
The keys `magnetic_flux_density`, `rotor_angle`, and `rotor_frequency`
collectively describe the spin-environment, while the keys `number_of_points`,
`spectral_width`, and `reference_offset` describes the grid coordinates
along the dimension at which the spectrum is evaluated.

..  .. seealso:: :ref:`dimension`.


An instance of a Dimension object may be created from a python dictionary,
such as the one shown above, using the
:meth:`~mrsimulator.Dimension.parse_dict_with_units` method from the
:ref:`dimension_api` class as follows,

    >>> from mrsimulator import Dimension
    >>> dimension_object = Dimension.parse_dict_with_units(dimension)

You may create multiple dimension objects as required by the
experiment. In this example, we stick with a single dimension.
Finally, add the dimensions, in this case, ``dimension_object``,
to the instance of the Simulator class, ``sim``, as follows,

    >>> sim.dimensions += [dimension_object]

Setting up the NMR method
-------------------------

Besides, setting up the list of isotopomer and dimension objects, you also need
to specify an NMR method that will be used in generating the line-shape.
Note, while the list of isotopomer objects are independent of the NMR method,
the ordered list of dimension objects dependents on the specified NMR method.
In this example, we illustrate the use of a single pulse acquisition method,
referred here as `one_d_spectrum`. This method requires a single
dimension.

.. seealso:: :ref:`methods_api`

Import the method as

.. doctest::

    >>> from mrsimulator.methods import one_d_spectrum

Running the simulator
---------------------

To simulate the line-shape, run the simulator with the
:meth:`~mrsimulator.simulator.Simulator.run` method, as follows,

.. doctest::

    >>> freq, amp = sim.run(method=one_d_spectrum)

In the above code, ``freq`` and ``amp`` are the dimensionless frequency
ratio given in `ppm` and the corresponding amplitude of the spectrum. The
following is a figure of the above lineshape plotted using the matplotlib
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

    >>> plot(freq, amp) # doctest: +SKIP

.. .. testsetup::

..    >>> plot_save(freq, amp, "example.pdf")

.. figure:: _images/example.*
    :figclass: figure-polaroid
