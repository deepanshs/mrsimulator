

.. _using_objects:

=================================================
Getting started with `Mrsimulator`: Using objects
=================================================

In the previous section on getting started, we showed an example where
we parse python dictionaries to create :ref:`isotopomer_api` and
:ref:`dimension_api` objects. In this section, we'll illustrate
the use of core `Mrsimulator` objects to construct isotopomer and dimension
objects.

Let's start by importing the objects.

.. doctest::

    >>> from mrsimulator import Simulator, Isotopomer, Site, Dimension
    >>> from mrsimulator import SymmetricTensor as st

.. note::
    Unlike python dictionary objects from our last example, when using
    Mrsimulator objects, the attribute value is given as a number rather than
    a string with a number and a unit. We assume default units for class
    attributes. To learn more about the default units, please refer to the
    documentation of the respective class.
    For the convenience of our users, we have added an attribute,
    ``property_units``, to every class that holds the default unit of the
    respective class attributes.

Site object
-----------
As the name suggests, a :ref:`site_api` object is used in creating sites. For
example,

.. doctest::

    >>> C13A = Site(isotope='13C')

The above code creates a site with :math:`^{13}\text{C}` isotope. Because, no
further information is delivered to the site object, other attributes such as
the isotropic chemical shift assume their default values.

.. doctest::

    >>> C13A.isotropic_chemical_shift # value is given in ppm
    0

Here, the isotropic chemical shift is given in ppm. This information is also
delivered with the ``property_units`` attribute of the instance. For example,

.. doctest::

    >>> C13A.property_units
    {'isotropic_chemical_shift': 'ppm'}

Let's create few more sites.

.. doctest::

    >>> C13B = Site(isotope='13C', isotropic_chemical_shift=-10)
    >>> H1 = Site(isotope='1H', shielding_symmetric=st(zeta=5.1, eta=0.1))
    >>> O17 = Site(isotope='17O', isotropic_chemical_shift=41.7, quadrupolar=st(Cq=5.15e6, eta=0.21))

The site, ``C13B``, is the second :math:`^{13}\text{C}` site with an isotropic
chemical shift of -10 ppm.

In creating the site, ``H1``, we use the :ref:`symmetric_tensor_api` object to
describe a traceless symmetric second-rank irreducible nuclear shielding
tensor, using the attributes `zeta` and `eta`, respectively.
The parameter `zeta` and `eta` are defined as per
Haeberlen convention and describe the anisotropy and asymmetry parameters of
the tensor, respectively.
Once again, the default unit for the attributes from the SymmetricTensor class
may be found in the ``property_units`` attribute of the respective class,
such as,

.. doctest::

    >>> H1.shielding_symmetric.property_units
    {'zeta': 'ppm', 'alpha': 'rad', 'beta': 'rad', 'gamma': 'rad'}

For site, ``O17``, we once again use the SymmetricTensor class, only this time
to describe a traceless symmetric second-rank irreducible electric quadrupole
tensor, using the attributes `Cq` and `eta`, respectively. The parameter `Cq`
is the quadrupole coupling constant, and `eta` is the asymmetry parameters of
the quadrupole tensor, respectively. Once again, the default unit for the
attributes from the SymmetricTensor object can be found in the
``property_units`` attribute of the respective class, such as,

.. doctest::

    >>> O17.quadrupolar.property_units
    {'Cq': 'Hz', 'alpha': 'rad', 'beta': 'rad', 'gamma': 'rad'}

.. .. tip::
..     If `Mrsimulator`, we use SymmetricTensor objects in about every solid-state
..     NMR line-shape simulation. We recommend reading :ref:`symmetric_tensor_api`.


Isotopomer object
-----------------

An isotopomer object contains sites and couplings along with the abundance
of the said isotopomer. In this version, we focus on isotopomers with a single
site, and therefore the couplings are irrelevant. Let's use the
sites we have already created to set up isotopomers.

.. doctest::

    >>> isotopomer_1 = Isotopomer(name='C13A', sites=[C13A], abundance=20)
    >>> isotopomer_2 = Isotopomer(name='C13B', sites=[C13B], abundance=56)
    >>> isotopomer_3 = Isotopomer(name='H1', sites=[H1], abundance=100)
    >>> isotopomer_4 = Isotopomer(name='O17', sites=[O17], abundance=1)

Here we have created four isotopomers, each with a single site.


Dimension object
----------------
Likewise, we can create a :ref:`dimension_api` object following,

.. doctest::

    >>> dimension_1 = Dimension(isotope='13C', number_of_points=2046, spectral_width=25000)

Here, we have set up a dimension, ``dimension_1``, which ready to record
:math:`^{13}\text{C}` resonances over 25 kHz with 2046 points. The
unspecified attributes, such as `rotor_frequency`, `rotor_angle`,
`magnetic_flux_density` assume their default values,

.. doctest::

    >>> dimension_1.rotor_frequency # in Hz (static)
    0
    >>> dimension_1.rotor_angle # magic angle in rad
    0.9553166
    >>> dimension_1.magnetic_flux_density # in T
    9.4

The unit associated with the attributes from the dimension object can be
similarly accessed via,

.. doctest::

    >>> dimension_1.property_units
    {'spectral_width': 'Hz',
     'reference_offset': 'Hz',
     'magnetic_flux_density': 'T',
     'rotor_frequency': 'Hz',
     'rotor_angle': 'rad'}

Simulator object
----------------

The use of simulator object is the same as described in the previous
section.

.. doctest::

    >>> sim = Simulator()
    >>> sim.isotopomers += [isotopomer_1, isotopomer_2, isotopomer_3, isotopomer_4]
    >>> sim.dimensions += [dimension_1]

A quick run
-----------

Let's import the `one_d_spectrum` method, do a quick run of the simulator,
and observe the spectrum. But before, here is the plotting method we'll
use to plot the spectrum from the following simulations.

    >>> import matplotlib.pyplot as plt
    >>> def plot(x, y):
    ...     plt.figure(figsize=(4, 3))
    ...     plt.plot(x, y, linewidth=1)
    ...     plt.xlim([x.value.max(), x.value.min()])
    ...     plt.xlabel(f"frequency ratio / {str(x.unit)}")
    ...     plt.grid(color='gray', linestyle='--', linewidth=1.0, alpha=0.25)
    ...     plt.tight_layout(h_pad=0, w_pad=0, pad=0)

And now, a quick run.

.. doctest::

    >>> from mrsimulator.methods import one_d_spectrum
    >>> x, y = sim.run(method=one_d_spectrum)
    >>> plot(x,y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_1')

.. figure:: _images/example_1.*
    :figclass: figure-polaroid

Tweak the sites and simulate
****************************

Let's add shielding tensors to sites ``C13A`` and ``C13B``.

.. doctest::

    >>> C13A.shielding_symmetric = st(zeta=80, eta=0.5)
    >>> C13B.shielding_symmetric = st(zeta=-100, eta=0.25)
    >>> x, y = sim.run(method=one_d_spectrum)
    >>> plot(x,y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_2')

.. figure:: _images/example_2.*
    :figclass: figure-polaroid

.. note::
    Because the objects in python are passed as reference, we were able to
    modify the ``C13A`` and ``C13B`` Site objects without having to reassemble
    the isotopomer or dimension objects.

Turn up the rotor frequency and simulate
****************************************

.. doctest::

    >>> dimension_1.rotor_frequency = 1000 # in Hz
    >>> x, y = sim.run(method=one_d_spectrum)
    >>> plot(x,y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_3')

.. figure:: _images/example_3.*
    :figclass: figure-polaroid


Change the rotor angle and simulate
***********************************

.. doctest::

    >>> dimension_1.rotor_angle = 90*3.1415926/180 # 90 degree in radian
    >>> x, y = sim.run(method=one_d_spectrum)
    >>> plot(x,y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_4')

.. figure:: _images/example_4.*
    :figclass: figure-polaroid

Switch to 1H and simulate
*************************

.. doctest::

    >>> dimension_1.isotope = '1H'
    >>> x, y = sim.run(method=one_d_spectrum)
    >>> plot(x,y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_5')

.. figure:: _images/example_5.*
    :figclass: figure-polaroid

Switch to 17O and simulate
**************************

.. doctest::

    >>> dimension_1.isotope = '17O'
    >>> dimension_1.rotor_angle = 0.9553166 # magic angle is rad
    >>> dimension_1.rotor_frequency = 15000 # Hz
    >>> x, y = sim.run(method=one_d_spectrum)
    >>> plot(x,y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_6')

.. figure:: _images/example_6.*
    :figclass: figure-polaroid
