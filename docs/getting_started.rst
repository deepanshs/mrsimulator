
.. _getting_started:

.. testsetup::

    >>> import matplotlib
    >>> font = {'family': 'Helvetica', 'weight': 'light', 'size': 9};
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
:attr:`~mrsimulator.Simulator.spectrum`. The default value of these
attributes is an empty list.

.. doctest::

    >>> sim.isotopomers
    []
    >>> sim.spectrum
    []

The value of the attribute :attr:`~mrsimulator.Simulator.isotopomers` is a
list of :ref:`isotopomer_api` objects while the value of the attribute
:attr:`~mrsimulator.Simulator.spectrum` is a list of
:ref:`spectroscopicDimension_api` objects. Before you can start simulating
NMR lineshapes, you need to understand the role of the Isotopomer and
SpectroscopicDimension objects. We recommend starting with
:ref:`dictionary_objects` and :ref:`spectroscopic_dimension`.


Setting up Isotopomer objects
-----------------------------
We define an isotopomer as an isolated spin-system containing
multiple sites and couplings between them. In this example, however, we
concern ourselves with a single site spin-system, that is, an isotopomer with
a single site. Shown below is an example of such isotopomer expressed as a
python dictionary.

.. code-block:: py

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

In the above example, ``isotopomer_dict``, represents an isotopomer.
This isotopomer contains a ``sites`` keyword whose value is a list of sites
defined within the isotopomer. In this examples, we have defined a single site,
again as a python dictionary, containing site specific information such as,
the site isotope, the isotropic chemical shift, and the parameters from the
irreducible second rank symmetric nuclear shielding tensor---asymmetry
(:math:`\zeta`), and anisotropy (:math:`\eta`), expressed in Haeberlen
convension. For additional information see :ref:`isotopomer` and :ref:`site`.

.. seealso:: :ref:`dictionary_objects`

An instance of an Isotopomer object may be created from a python dictionary,
such as the one shown above, using the :ref:`isotopomer_api`
class.

    >>> from mrsimulator import Isotopomer
    >>> isotopomer_object = Isotopomer.parse_dict_with_units(isotopomer_dict)

You may create as many isotopomer objects as necessary, although in this
example, we stick with a single isotopomer. Finally, add the isotopomer
objects, in this case, ``isotopomer_object``, to the instance of the Simulator
class, ``sim``, as follows

    >>> sim.isotopomers = [isotopomer_object]


Setting up SpectroscopicDimension objects
-----------------------------------------

A spectroscopic dimension is a dimension of the NMR spectrum. The number of
spectroscopic dimensions depends on the dimensionality of the NMR experiment.
A one pulse acquired spectrum will consist of a single spectroscopic dimension,
while a two-dimensional experiments will consist of two spectroscopic
dimensions. In `mrsimulator`, a spectroscopic dimension includes keywords that
are required in evaluating the spectrum/line-shape along the dimension.
A spectroscopic dimension may be defined as follows,

.. code-block:: py

    >>> dimension = {
    ...     "isotope": "29Si",
    ...     "magnetic_flux_density": "9.4 T",
    ...     "rotor_angle": "54.735 deg",
    ...     "rotor_frequency": "0 kHz",
    ...     "number_of_points": 2048,
    ...     "spectral_width": "25 kHz",
    ...     "reference_offset": "8 kHz"
    ... }

In the above example, ``dimension``, is a spectroscopic dimension represented
as a python dictionary. Here, the value of the `isotope` key is the isotope
symbol of the observed nucleus. The keys `magnetic_flux_density`,
`rotor_angle`, and `rotor_frequency` define the spin-environment, while
`number_of_points`, `spectral_width`, and `reference_offset` define the
grid points along the spectroscopic dimension at which the spectrum is
evaluated.

.. seealso:: :ref:`spectroscopic_dimension`.


An instance of a spectroscopic dimension object may be created from
a python dictionary, such as the one shown above, using the
:ref:`spectroscopicDimension_api` class.

    >>> from mrsimulator import SpectroscopicDimension
    >>> spectrum_object = SpectroscopicDimension.parse_dict_with_units(dimension)

You may create multiple spectroscopic dimension objects as required by the
experiment. In this example, we stick with a single spectroscopic dimension.
Finally, add the spectroscopic dimensions, in this case, ``spectrum_object``,
to the instance of the Simulator class, ``sim``, as follows

    >>> sim.spectrum = [spectrum_object]

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

    >>> freq, amp = sim.run(method=one_d_spectrum, verbose=1)
    `one_d_spectrum` method simulation parameters.
    ---------------------------------------------
    The magnetic flux density is 9.4 T.
    Sample rotation angle is 0.9553059660790962 rad.
    Sample rotation frequency is 0.0 Hz.
    Simulating 29Si(I=0.5, precession frequency = -79.571 MHz) isotope.
    Recording 29Si spectrum with 2048 points over 25000.0 Hz bandwidth
    and a reference offset of 8000.0 Hz.
    <BLANKLINE>
    29Si site 0 from isotopomer 0 @ 100% abundance
    ----------------------------------------------
    Isotropic chemical shift (δ) = 8044.628099999999 Hz
    Shielding anisotropy (ζ) = -5609.7555 Hz
    Shielding asymmetry (η) = 0.5
    Shielding orientation = [alpha = 0.0, beta = 0.0, gamma = 0.0]

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

    >>> plot(freq, amp)

.. testsetup::

    >>> plot_save(freq, amp, "example.pdf")

.. figure:: _images/example.*
    :figclass: figure-polaroid
