

.. _using_objects:

.. .. image:: https://mybinder.org/badge_logo.svg
..  :target: https://mybinder.org/v2/gh/DeepanshS/mrsimulator/master?filepath=jupyternotebooks%2F

=================================================
Getting started with `Mrsimulator`: Using objects
=================================================

In the previous section on getting started, we show an example where
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

Let's create a few more sites.

.. doctest::

    >>> C13B = Site(isotope='13C', isotropic_chemical_shift=-10)
    >>> H1 = Site(isotope='1H', shielding_symmetric=dict(zeta=5.1, eta=0.1))
    >>> O17 = Site(isotope='17O', isotropic_chemical_shift=41.7, quadrupolar=dict(Cq=5.15e6, eta=0.21))

The site, ``C13B``, is the second :math:`^{13}\text{C}` site with an isotropic
chemical shift of -10 ppm.

In creating the site, ``H1``, we use the dictionary object to
describe a traceless symmetric second-rank irreducible nuclear shielding
tensor, using the attributes `zeta` and `eta`, respectively.
The parameter `zeta` and `eta` are defined as per
Haeberlen convention and describes the anisotropy and asymmetry parameter of
the tensor, respectively.
The default unit of the attributes from the `shielding_symmetric`
may be found with the ``property_units`` attribute, such as

.. doctest::

    >>> H1.shielding_symmetric.property_units
    {'zeta': 'ppm', 'alpha': 'rad', 'beta': 'rad', 'gamma': 'rad'}

For site, ``O17``, we once again make use of the dictionary object, only this time
to describe a traceless symmetric second-rank irreducible electric quadrupole
tensor, using the attributes `Cq` and `eta`, respectively. The parameter `Cq`
is the quadrupole coupling constant, and `eta` is the asymmetry parameters of
the quadrupole tensor, respectively.
The default unit of these attributes is once again found with the ``property_units``
attribute,

.. doctest::

    >>> O17.quadrupolar.property_units
    {'Cq': 'Hz', 'alpha': 'rad', 'beta': 'rad', 'gamma': 'rad'}


Isotopomer object
-----------------

An isotopomer object contains sites and couplings along with the abundance
of the respective isotopomer. In this version, we focus on isotopomers with a
single site, and therefore the couplings are irrelevant.

Let's use the sites we have already created to set up four isotopomers.

.. doctest::

    >>> isotopomer_1 = Isotopomer(name='C13A', sites=[C13A], abundance=20)
    >>> isotopomer_2 = Isotopomer(name='C13B', sites=[C13B], abundance=56)
    >>> isotopomer_3 = Isotopomer(name='H1', sites=[H1], abundance=100)
    >>> isotopomer_4 = Isotopomer(name='O17', sites=[O17], abundance=1)


Method object
-------------
Likewise, we can create a :class:`~mrsimulator.methods.BlochDecaySpectrum`
object following,

.. doctest::

    >>> from mrsimulator.methods import BlochDecaySpectrum
    >>> method_1 = BlochDecaySpectrum(
    ...     channels=["13C"],
    ...     spectral_dimensions = [dict(count=2046, spectral_width=25000)] # spectral_width is in Hz.
    ... )

The above method, ``method_1``, is declared to record :math:`^{13}\text{C}`
resonances over 25 kHz spectral width using 2046 points. The
unspecified attributes, such as `rotor_frequency`, `rotor_angle`,
`magnetic_flux_density`, assume their default value.


Simulator object
----------------

The use of the simulator object is the same as described in the previous
section.

.. doctest::

    >>> sim = Simulator()
    >>> sim.isotopomers += [isotopomer_1, isotopomer_2, isotopomer_3, isotopomer_4] # add isotopomers
    >>> sim.methods += [method_1] # add method

A quick run
-----------

Let's do a quick run of the simulator object, and observe the spectrum. But before,
here is the plotting script we'll use to plot the spectrum for all further examples.

.. doctest::

    >>> import matplotlib.pyplot as plt
    >>> def plot(csdm):
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

Notice, we added four isotopomers to the ``sim`` object, two :math:`^{13}\text{C}`, one
:math:`^1\text{H}`, and one :math:`^{17}O` site along with a BlochDecaySpectrum method
with a :math:`^{13}\text{C}` channel. When you run the simulation, only the resonance
from the given channel will be recorded, as seen from the above plot, where just the two
:math:`^{13}\text{C}` isotropic chemical shifts are observed.


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

Let's turn up the rotor frequency from 0 Hz to 1 kHz. Note, that we do not add another
method to the ``sim`` object, but update the method at index 0 with a new method.

.. doctest::

    >>> # Update the method object at index 0.
    >>> sim.methods[0] = BlochDecaySpectrum(
    ...     channels=["13C"],
    ...     rotor_frequency=1000, # in Hz
    ...     spectral_dimensions=[dict(count=2046, spectral_width=25000)] # spectral_width is in Hz.
    ... )
    >>> sim.run()
    >>> plot(sim.methods[0].simulation) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(*sim.methods[0].simulation.to_list(), 'example_3')

.. figure:: _images/example_3.*
    :figclass: figure-polaroid

    An example of the solid-state :math:`^{13}\text{C}` MAS sideband simulation.

Change the rotor angle and re-simulate
**************************************

Let's also set the rotor angle from magic angle to 90 degrees. Again, we update the
method at index 0.

.. doctest::

    >>> # Update the method object at index 0.
    >>> sim.methods[0] = BlochDecaySpectrum(
    ...     channels=["13C"],
    ...     rotor_frequency=1000, # in Hz.
    ...     rotor_angle=90*3.1415926/180, # 90 degree in radians.
    ...     spectral_dimensions=[dict(count=2046, spectral_width=25000)] # spectral_width is in Hz.
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
    ...     spectral_dimensions=[dict(count=2046, spectral_width=25000)]
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

Likewise, update the value of the `channels` attribute to `17O`.

.. doctest::

    >>> sim.methods[0] = BlochDecaySpectrum(
    ...     channels=["17O"],
    ...     rotor_frequency= 15000, # in Hz.
    ...     rotor_angle = 0.9553166, # magic angle is rad.
    ...     spectral_dimensions = [dict(count=2046, spectral_width=25000)]
    ... )
    >>> sim.run()
    >>> plot(sim.methods[0].simulation) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(*sim.methods[0].simulation.to_list(), 'example_6')

.. figure:: _images/example_6.*
    :figclass: figure-polaroid

    An example of the solid-state :math:`^{17}\text{O}` MAS central-transition
    simulation.

If you are familiar with the quadrupolar line-shapes, you may immediately associate
this simulation to a second-order quadrupolar line-shape of the central transition.
You may also notice some unexpected resonances around 50 ppm and -220 ppm. These
unexpected resonances are the spinning sidebands of the satellite transitions.
Note, the BlochDecaySpectrum method computes resonances from all transitions with
:math:`p = \Delta m = -1`.

Let's see what transition pathways are used in our simulation. Use the
:meth:`~mrsimulator.Method.get_transition_pathways` function of the Method instance to
see the list of transition pathways, for example,

.. doctest::

    >>> print(sim.methods[0].get_transition_pathways(isotopomer_1)) # 13C
    [[|-0.5⟩⟨0.5|]]

    >>> print(sim.methods[0].get_transition_pathways(isotopomer_2)) # 13C
    [[|-0.5⟩⟨0.5|]]

    >>> print(sim.methods[0].get_transition_pathways(isotopomer_3)) # 1H
    [[|-0.5⟩⟨0.5|]]

    >>> print(sim.methods[0].get_transition_pathways(isotopomer_4)) # 17O
    [[|-2.5⟩⟨-1.5|]
     [|-1.5⟩⟨-0.5|]
     [|-0.5⟩⟨0.5|]
     [|0.5⟩⟨1.5|]
     [|1.5⟩⟨2.5|]]

Notice, there are five transition pathways for the :math:`^{17}\text{O}` site, one
associated with the central-transition, two with the inner-satellites, and two with
the outer-satellites. For central transition selective simulation, use the
:class:`~mrsimulator.methods.BlochDecayCentralTransitionSpectrum` method.

.. doctest::

    >>> from mrsimulator.methods import BlochDecayCentralTransitionSpectrum
    >>> method_2 = BlochDecayCentralTransitionSpectrum(
    ...     channels=["17O"],
    ...     rotor_frequency= 15000, # in Hz.
    ...     rotor_angle = 0.9553166, # magic angle is rad.
    ...     spectral_dimensions = [dict(count=2046, spectral_width=25000)]
    ... )

    >>> # the transition pathways
    >>> print(method_2.get_transition_pathways(isotopomer_4)) # 17O
    [[|-0.5⟩⟨0.5|]]

Now, you may simulate the central transition selective spectrum.

    >>> sim.methods += [method_2]
    >>> sim.run()
    >>> plot(sim.methods[1].simulation) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(*sim.methods[1].simulation.to_list(), 'example_7')

.. figure:: _images/example_7.*
    :figclass: figure-polaroid
