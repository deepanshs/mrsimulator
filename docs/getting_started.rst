
.. _getting_started:

.. >>> font = {'family': 'Helvetica', 'weight': 'light', 'size': 9};
.. >>> matplotlib.rc('font', **font)

.. testsetup::

    >>> import matplotlib
    >>> from os import path

    >>> import matplotlib.pyplot as plt
    >>> def plot_save(x, y, filename):
    ...     plt.figure(figsize=(4, 3))
    ...     plt.plot(x, y, linewidth=1)
    ...     plt.xlim([x.value.max(), x.value.min()])
    ...     plt.xlabel(f"frequency ratio / {str(x.unit)}")
    ...     plt.grid(color='gray', linestyle='--', linewidth=1.0, alpha=0.25)
    ...     plt.tight_layout(h_pad=0, w_pad=0, pad=0)
    ...
    ...     filename = path.split(filename)[1]
    ...     filepath = './docs/_images'
    ...     pth = path.join(filepath, filename)
    ...     plt.savefig(pth+'.pdf')
    ...     plt.savefig(pth+'.png', dpi=100)
    ...     plt.close()


==================================
Getting started with `mrsimulator`
==================================

We have put together a set of guidelines for using various methods and
attributes of `mrsimulator` package. We encourage the users
to follow these guidelines to promote consistency amongst others.
In `mrsimulator`, the solid state nuclear magnetic resonance (ssNMR)
lineshape is calculated through an instance of the :ref:`simulator_api`
class.

Import the :ref:`simulator_api` class,

.. doctest::

    >>> from mrsimulator import Simulator

and create an instance as follows,

.. doctest::

    >>> sim = Simulator()

Here, ``sim`` is a variable with an instance of the :ref:`simulator_api`
class. The two attributes of this class that you will often use are
:attr:`~mrsimulator.Simulator.isotopomers` and
:attr:`~mrsimulator.Simulator.dimensions`, whose value is a list of
:ref:`isotopomer_api` and :ref:`dimension_api` objects,
respectively. The default value of these attributes is an empty list.

.. doctest::

    >>> sim.isotopomers
    []
    >>> sim.dimensions
    []


Before you can start simulating
NMR lineshapes, you need to understand the role of the Isotopomer and
Dimension objects. We recommend starting with
:ref:`dictionary_objects` and :ref:`dimension`.


Setting up Isotopomer objects
-----------------------------
For all practical purposes, an isotopomer is an isolated spin-system with
multiple sites and couplings between them. For simplicity, in this section,
we concern ourselves with a single site spin-system, that is, an
isotopomer with a single site. Shown below is an example of a single-site
isotopomer, expressed as a python dictionary.

.. code-block:: py
    :linenos:

    >>> isotopomer_dict = {
    ...     "sites": [
    ...         {
    ...             "isotope": "29Si",
    ...             "isotropic_chemical_shift": "-101.1 ppm",
    ...             "shielding_symmetric": {
    ...                 "zeta": "70.5 ppm",
    ...                 "eta": 0.5
    ...             }
    ...         }
    ...     ]
    ... }

The above isotopomer contains a ``sites`` keyword, at line 2, whose value is
a list of sites defined within the isotopomer. In this examples, we have
defined a single site, again as a python dictionary, (lines 3-10) containing
site specific information such as, the site isotope (line 4), the isotropic
chemical shift (line 5), and the parameters from the second-rank irreducible
symmetric nuclear shielding tensor---the shielding asymmetry (:math:`\zeta`)
at line 7, and the shielding anisotropy (:math:`\eta`) at line 8, expressed
using Haeberlen convension.
For additional information see :ref:`isotopomer` and :ref:`site`.

.. seealso:: :ref:`dictionary_objects`

An instance of an Isotopomer object may be created from a python dictionary,
such as the one shown above, using the :ref:`isotopomer_api`
class.

    >>> from mrsimulator import Isotopomer
    >>> isotopomer_object = Isotopomer.parse_dict_with_units(isotopomer_dict)

You may create as many isotopomer objects as necessary, although in this
example, we stick with a single isotopomer. Finally, add the isotopomer
objects, in this case, the variable ``isotopomer_object``, to the instance of
the Simulator class, ``sim``, as follows

    >>> sim.isotopomers.append(isotopomer_object)


Setting up Dimension objects
----------------------------

A spectroscopic dimension is a dimension of the NMR spectrum. The number of
spectroscopic dimensions depends on the dimensionality of the experiment. For
example, a one-pulse acquired 1-D spectrum consists of a single spectroscopic
dimension, while two-dimensional experiments will consist of two spectroscopic
dimensions. In `mrsimulator`, we have designed the spectroscopic dimension to
includes keywords that are required in evaluating the spectrum/line-shape along
the dimension. A spectroscopic dimension may be defined as follows,

.. code-block:: py
    :linenos:

    >>> dimension = {
    ...     "isotope": "29Si",
    ...     "magnetic_flux_density": "9.4 T",
    ...     "rotor_angle": "54.735 deg",
    ...     "rotor_frequency": "0 kHz",
    ...     "number_of_points": 2048,
    ...     "spectral_width": "25 kHz",
    ...     "reference_offset": "-8 kHz"
    ... }

In the above example, the variable ``dimension``, holds a python dictionary
representation of the spectroscopic dimension. Here, the value of the `isotope`
key is the isotope symbol of the observed nucleus. A value, ``29Si``, means
that the simulated lineshape arises from :math:`^{29}\text{Si}` resonances.
The keys `magnetic_flux_density`, `rotor_angle`, and `rotor_frequency`
collectively define the spin-environment, while the keys `number_of_points`,
`spectral_width`, and `reference_offset` describes the grid coordinates
along the spectroscopic dimension at which the spectrum is evaluated.

.. seealso:: :ref:`dimension`.


An instance of a spectroscopic dimension object may be created from
a python dictionary, such as the one shown above, using the
:ref:`dimension_api` class.

    >>> from mrsimulator import Dimension
    >>> spectrum_object = Dimension.parse_dict_with_units(dimension)

You may create multiple spectroscopic dimension objects as required by the
experiment. In this example, we stick with a single spectroscopic dimension.
Finally, add the spectroscopic dimensions, in this case, ``spectrum_object``,
to the instance of the Simulator class, ``sim``, as follows

    >>> sim.dimensions = [spectrum_object]

Setting up the NMR method
-------------------------

Beside, the list of isotopomer and spectroscopic dimension objects,
`mrsimulator` also requires an NMR method to simulate a line-shape.
Note, while the list isotopomer objects are independent of the NMR method, the
ordered list of spectroscopic dimension objects is dependent on the NMR method.
In this example, we illustrate the use of one pulse acquisition method,
referred here as, `one_d_spectrum`. This method requires a single spectroscopic
dimension.

.. seealso:: :ref:`methods_api`

Import the method as

.. doctest::

    >>> from mrsimulator.methods import one_d_spectrum

and run the simulation using

.. doctest::

    >>> freq, amp = sim.run(method=one_d_spectrum)

In the above code, the ``freq`` and ``amp`` are the dimensionless frequency
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
