

.. _using_objects:

.. .. image:: https://mybinder.org/badge_logo.svg
..  :target: https://mybinder.org/v2/gh/DeepanshS/mrsimulator/master?filepath=jupyternotebooks%2F

=================================================
Getting started with `Mrsimulator`: Using objects
=================================================

In the previous section on getting started, we showed an example where
we parse python dictionaries to create :ref:`isotopomer_api` and
:ref:`method_api` objects. In this section, we'll illustrate how we can
achieve the same result using the core `Mrsimulator` objects.

Let's start by importing the objects.

.. doctest::

    >>> from mrsimulator import Simulator, Isotopomer, Site
    >>> from mrsimulator.methods import BlochDecaySpectrum

.. note::
    Unlike python dictionary objects from our last example, when using
    Mrsimulator objects, the attribute value is given as a number rather than
    a string with a number and a unit. We assume default units for the class
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

The above code creates a site with a :math:`^{13}\text{C}` isotope. Because, no
further information is delivered to the site object, other attributes such as
the isotropic chemical shift assume their default value.

.. doctest::

    >>> C13A.isotropic_chemical_shift # value is given in ppm
    0

Here, the isotropic chemical shift is given in ppm. This information is also
present in the ``property_units`` attribute of the instance. For example,

.. doctest::

    >>> C13A.property_units
    {'isotropic_chemical_shift': 'ppm'}

Let's create few more sites.

.. doctest::

    >>> C13B = Site(isotope='13C', isotropic_chemical_shift=-10)
    >>> H1 = Site(isotope='1H', shielding_symmetric=dict(zeta=5.1, eta=0.1))
    >>> O17 = Site(isotope='17O', isotropic_chemical_shift=41.7, quadrupolar=dict(Cq=5.15e6, eta=0.21))

The site, ``C13B``, is the second :math:`^{13}\text{C}` site with an isotropic
chemical shift of -10 ppm.

In creating the site, ``H1``, we use the dict object to
describe a traceless symmetric second-rank irreducible nuclear shielding
tensor, using the attributes `zeta` and `eta`, respectively.
The parameter `zeta` and `eta` are defined as per
Haeberlen convention and describes the anisotropy and asymmetry parameter of
the tensor, respectively.
Similarly, the default unit of the attributes from the SymmetricTensor class
may be found with the ``property_units`` attribute, such as

.. doctest::

    >>> H1.shielding_symmetric.property_units
    {'zeta': 'ppm', 'alpha': 'rad', 'beta': 'rad', 'gamma': 'rad'}

For site, ``O17``, we once again make use of the dict object, only this time
to describe a traceless symmetric second-rank irreducible electric quadrupole
tensor, using the attributes `Cq` and `eta`, respectively. The parameter `Cq`
is the quadrupole coupling constant, and `eta` is the asymmetry parameters of
the quadrupole tensor, respectively.
The default unit of the attributes from the SymmetricTensor class may be found
in the ``property_units`` attribute, such as

.. doctest::

    >>> O17.quadrupolar.property_units
    {'Cq': 'Hz', 'alpha': 'rad', 'beta': 'rad', 'gamma': 'rad'}


Isotopomer object
-----------------

An isotopomer object contains sites and couplings along with the abundance
of the respective isotopomer. In this version, we focus on isotopomers with a
single site, and therefore the couplings are irrelevant.

Let's use the sites we have already created to set up isotopomers.

.. doctest::

    >>> isotopomer_1 = Isotopomer(name='C13A', sites=[C13A], abundance=20)
    >>> isotopomer_2 = Isotopomer(name='C13B', sites=[C13B], abundance=56)
    >>> isotopomer_3 = Isotopomer(name='H1', sites=[H1], abundance=100)
    >>> isotopomer_4 = Isotopomer(name='O17', sites=[O17], abundance=1)

Here, we have created four isotopomers, each with a single site.


Method object
-------------
Likewise, we can create a :class:`~mrsimulator.methods.BlochDecaySpectrum`
object following,

.. doctest::

    >>> from mrsimulator.methods import BlochDecaySpectrum
    >>> method_1 = BlochDecaySpectrum(
    ...     channels=["13C"],
    ...     dimensions = [dict(count=2046, spectral_width=25000)]
    ... )

The above method, ``method_1``, is defined to record :math:`^{13}\text{C}`
resonances over 25 kHz spectral width using 2046 points. The
unspecified attributes, such as `rotor_frequency`, `rotor_angle`,
`magnetic_flux_density`, assume their default value.


Simulator object
----------------

The use of the simulator object is the same as described in the previous
section.

.. doctest::

    >>> sim = Simulator()
    >>> # add isotopomers
    >>> sim.isotopomers += [isotopomer_1, isotopomer_2, isotopomer_3, isotopomer_4]
    >>> # add methods
    >>> sim.methods += [method_1]

A quick run
-----------

Let's do a quick run of the simulator object, and observe the spectrum. But before,
here is the plotting method we'll use to plot the spectrum for all further examples.

.. doctest::

    >>> import matplotlib.pyplot as plt
    >>> def plot(csdm):
    ...     # Convert the coordinates unit from Hz to ppm
    ...     csdm.dimensions[0].to('ppm', 'nmr_frequency_ratio')
    ...     x, y = csdm.to_list()
    ...     plt.figure(figsize=(4.5, 2.5))
    ...     plt.plot(x, y, linewidth=1)
    ...     plt.xlim([x.value.max(), x.value.min()])
    ...     plt.xlabel(f"frequency ratio / {str(x.unit)}")
    ...     plt.grid(color='gray', linestyle='--', linewidth=1.0, alpha=0.25)
    ...     plt.tight_layout()
    ...     plt.show()

And now, a quick run.

.. doctest::

    >>> sim.run()
    >>> plot(sim.methods[0].simulation) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(*sim.methods[0].simulation.to_list(), 'example_1')

.. figure:: _images/example_1.*
    :figclass: figure-polaroid

    An example of the solid-state :math:`^{13}\text{C}` isotropic lineshape
    simulation.

Tweak the sites and re-simulate
*******************************

Let's add shielding tensors to sites ``C13A`` and ``C13B``.

.. doctest::

    >>> C13A.shielding_symmetric = dict(zeta=80, eta=0.5)
    >>> C13B.shielding_symmetric = dict(zeta=-100, eta=0.25)
    >>> sim.run()
    >>> plot(sim.methods[0].simulation) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(*sim.methods[0].simulation.to_list(), 'example_2')

.. figure:: _images/example_2.*
    :figclass: figure-polaroid

    An example of the static-solid state :math:`^{13}\text{C}` CSA lineshape
    simulation.

.. note::
    Because the objects in python are passed by reference, we were able to
    modify the ``C13A`` and ``C13B`` Site objects without having to reassemble
    the isotopomer or method objects.

Turn up the rotor frequency and re-simulate
*******************************************

Let's turn up the rotor frequency from 0 Hz to 1 kHz.

.. doctest::

    >>> # Update the method object at index 0.
    >>> sim.methods[0] = BlochDecaySpectrum(
    ...     channels=["13C"],
    ...     rotor_frequency=1000, # in Hz
    ...     dimensions=[dict(count=2046, spectral_width=25000)]
    ... )
    >>> sim.run()
    >>> plot(sim.methods[0].simulation) # doctest:+SKIP

Here, we update the method object at index 0.

.. .. testsetup::
..     >>> plot_save(*sim.methods[0].simulation.to_list(), 'example_3')

.. figure:: _images/example_3.*
    :figclass: figure-polaroid

    An example of the solid-state :math:`^{13}\text{C}` MAS sideband simulation.

Change the rotor angle and re-simulate
**************************************

Let's also set the rotor angle from magic angle to 90 degrees.

.. doctest::

    >>> # Update the method object at index 0.
    >>> sim.methods[0] = BlochDecaySpectrum(
    ...     channels=["13C"],
    ...     rotor_frequency=1000, # in Hz.
    ...     rotor_angle=90*3.1415926/180, # 90 degree in radians.
    ...     dimensions=[dict(count=2046, spectral_width=25000)]
    ... )
    >>> sim.run()
    >>> plot(sim.methods[0].simulation) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(*sim.methods[0].simulation.to_list(), 'example_4')

.. figure:: _images/example_4.*
    :figclass: figure-polaroid

    An example of the solid-state :math:`^{13}\text{C}` VAS sideband simulation.

Switch to 1H and simulate
*************************

To switch the channel, update the value of the `channels` attribute of the
method. Here, we change the channel from `13C` to `1H`.

.. doctest::

    >>> # Update the method object at index 0.
    >>> sim.methods[0] = BlochDecaySpectrum(
    ...     channels=["1H"],
    ...     rotor_frequency=1000, # in Hz.
    ...     rotor_angle=90*3.1415926/180, # 90 degree in radians.
    ...     dimensions=[dict(count=2046, spectral_width=25000)]
    ... )
    >>> sim.run()
    >>> plot(sim.methods[0].simulation) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(*sim.methods[0].simulation.to_list(), 'example_5')

.. figure:: _images/example_5.*
    :figclass: figure-polaroid

    An example of solid-state :math:`^{1}\text{H}` VAS sideband simulation.

Switch to 17O and simulate
**************************

Likewise, update the `channels` to `17O`.

.. doctest::

    >>> sim.methods[0] = BlochDecaySpectrum(
    ...     channels=["17O"],
    ...     rotor_frequency= 15000, # in Hz.
    ...     rotor_angle = 0.9553166, # magic angle is rad.
    ...     dimensions = [dict(count=2046, spectral_width=25000)]
    ... )
    >>> sim.run()
    >>> plot(sim.methods[0].simulation) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(*sim.methods[0].simulation.to_list(), 'example_6')

.. figure:: _images/example_6.*
    :figclass: figure-polaroid

    An example of the solid-state :math:`^{17}\text{O}` MAS central-transition
    simulation.

Notice, the spinning sidebands from the satellite transition in the case of 17O
line-shapes. The BlochDecaySpectrum method selects all :math:`\Delta m = -1`
transitions from the isotopomer. For example,

.. doctest::

    >>> print(sim.methods[0].get_transition_pathways(isotopomer_4))
    [[|-2.5⟩⟨-1.5|]
     [|-1.5⟩⟨-0.5|]
     [|-0.5⟩⟨0.5|]
     [|0.5⟩⟨1.5|]
     [|1.5⟩⟨2.5|]]

To only select central transition, use the :class:`~mrsimulator.methods.BlochDecayCentralTransitionSpectrum`
object.

.. doctest::

    >>> from mrsimulator.methods import BlochDecayCentralTransitionSpectrum
    >>> method_2 = BlochDecayCentralTransitionSpectrum(
    ...     channels=["17O"],
    ...     rotor_frequency= 15000, # in Hz.
    ...     rotor_angle = 0.9553166, # magic angle is rad.
    ...     dimensions = [dict(count=2046, spectral_width=25000)]
    ... )
    >>> print(method_2.get_transition_pathways(isotopomer_4))
    [[|-0.5⟩⟨0.5|]]

Now, you can simulate central transition spectrum.

    >>> sim.methods += [method_2]
    >>> sim.run()
    >>> plot(sim.methods[1].simulation) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(*sim.methods[1].simulation.to_list(), 'example_7')

.. figure:: _images/example_7.*
    :figclass: figure-polaroid
