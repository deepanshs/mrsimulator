

.. _config_simulator:

============================
Configuring Simulator object
============================

Up till now, we have been using the simulator object with the default setting.
In `Mrsimulator`, we choose the default settings such that it applies to a wide
range of simulations including, static, magic angle spinning (MAS), and
variable angle spinning (VAS) lineshapes. In certain situations, however, the
default settings are not sufficient to accurately represent the lineshapes. In
such cases, the user is advised to modify these settings as required. In the
following section, we briefly describe the configuration settings.

The :ref:`simulator_api` class is configured using the
:attr:`~mrsimulator.simulator.Simulator.config` attribute. The default value
of the config attributes is as follows,

.. doctest::

    >>> from mrsimulator import Simulator, Isotopomer, Dimension, Site
    >>> from mrsimulator.methods import one_d_spectrum
    >>> the_simulator = Simulator()

    >>> the_simulator.config
    ConfigSimulator(number_of_sidebands=64, integration_volume=octant, integration_density=70, decompose=False)

Here, the configurable attributes are ``number_of_sidebands``,
``integration_volume``, ``integration_density``, and ``decompose``.


Number of sidebands
-------------------
The value of this attribute is the number of sidebands
requested in evaluating the lineshapes. The default value is 64 and is
sufficient for most cases, as seen from our previous examples. In certain
circumstances, especially when the anisotropy is large or the rotor spin
frequency is low, 64 sidebands might not be sufficient. The user is thus
advised to increase this value as required. Consider the following example.

.. doctest::

    >>> simulator_1 = Simulator()

    >>> # create a site with a large anisotropy, 100 ppm.
    >>> Si29site = Site(isotope='29Si', shielding_symmetric={'zeta': 100, 'eta': 0.2})

    >>> # create a dimension with a low rotor frequency, 200 Hz
    >>> dimension = Dimension(isotope='29Si', spectral_width=25000, rotor_frequency=200)

    >>> simulator_1.isotopomers = [Isotopomer(sites=[Si29site])]
    >>> simulator_1.dimensions = [dimension]

    >>> # simulate and plot
    >>> x, y = simulator_1.run(method=one_d_spectrum)
    >>> plot(x, y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_sidebands_1')

.. figure:: _images/example_sidebands_1.*
    :figclass: figure-polaroid

    Spinning sidebands evaluated with a relatively low number of sidebands,
    resulting in an inaccurate sideband simulation.

If you are familiar with NMR lineshapes, you may notice that the above sideband
simulation is inaccurate, as evident from the abrupt termination of the
sideband amplitudes at the edges. As mentioned earlier, this
inaccuracy arises from evaluating a small number of sidebands relative to
the given anisotropy. Let's increase the number of sidebands to `90` and
observe.

.. doctest::

    >>> # set the number of sidebands to 90.
    >>> simulator_1.config.number_of_sidebands = 90
    >>> x, y = simulator_1.run(method=one_d_spectrum)
    >>> plot(x, y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_sidebands_2')

.. figure:: _images/example_sidebands_2.*
    :figclass: figure-polaroid

    Spinning sideband spectrum evaluated with a large number of sidebands.

Integration volume
------------------

Integration volume refers to the volume of the sphere over which the lineshape
is integrated. The default value is `octant`, i.e., the lineshape is integrated
over the positive octant of the sphere.
`Mrsimulator` enables the user to exploit the orientational symmetry of the
problem, and thus optimize the simulation by performing a partial integration
---`octant` or `hemisphere`. To learn more about the orientational symmetries,
please refer to Eden et. al. [#f4]_

In previous examples, we exploited the :math:`\text{D}_{2h}` symmetry
of the problem and therefore were able to simulate the spectrum by integrating
the line-shape over an octant. Observe what happens when this symmetry breaks.

Consider the :math:`^{29}\text{Si}` site, ``Si29site``, from the previous
example. This site has a symmetric shielding tensor with `zeta` and `eta` as
100 ppm and 0.2, respectively, giving a :math:`\text{D}_{2h}` symmetry to the
problem. We can break this symmetry by assigning Euler angles to this symmetric
shielding tensor, as follows,

.. doctest::

    >>> # add Euler angles to the shielding tensor.
    >>> Si29site.shielding_symmetric.alpha = 1.563 # in rad
    >>> Si29site.shielding_symmetric.beta = 1.2131 # in rad
    >>> Si29site.shielding_symmetric.gamma = 2.132 # in rad

    >>> # Let's observe the static spectrum which is more intuitive.
    >>> dimension.rotor_frequency = 0 # in Hz

    >>> # simulate and plot
    >>> x, y = simulator_1.run(method=one_d_spectrum)
    >>> plot(x, y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_integration_volume_1')

.. figure:: _images/example_integration_volume_1.*
    :figclass: figure-polaroid

    An example of an incomplete lineshape integration, lineshape simulation
    resulting from the frequency contributions evaluated over the positive
    octant.

Clearly, the above spectrum is incorrect. To fix this, set the integration
volume to `hemisphere` and re-simulate.

.. doctest::

    >>> # set integration volume to `hemisphere`.
    >>> simulator_1.config.integration_volume = 'hemisphere'

    >>> # simulate and plot
    >>> x, y = simulator_1.run(method=one_d_spectrum)
    >>> plot(x, y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_integration_volume_2')

.. figure:: _images/example_integration_volume_2.*
    :figclass: figure-polaroid

    The lineshape resulting from the frequency contributions evaluted over the
    top hemisphere.

Integration density
-------------------

Integration density controls the number of orientational points sampled over
the given volume. The NMR resonance frequency is then evaluated at these
orientations. The number of orientation is related to the value of this
attribute, `n`, following

``number_of_orientational_points = number_of_octants * (n + 1)(n + 2)/2``

where `number_of_octants` is the number of octants from the integration volume
attribute.

The default value, ``70``, produces 2556 orientations at which the NMR
frequency contribution is evaluated. The user may increase or decrease this
value as required by the problem.


Decompose
---------

Decompose is a boolean, if true, produces a series of spectra, each
arising from an individual isotopomer. For example,

.. doctest::

    >>> # Create two sites
    >>> site_A = Site(isotope='1H', shielding_symmetric={'zeta': 5, 'eta': 0.1})
    >>> site_B = Site(isotope='1H', shielding_symmetric={'zeta': -2, 'eta': 0.83})

    >>> # Create dimension object
    >>> dimension = Dimension(isotope='1H', spectral_width=10000)

    >>> # Create simulator object
    >>> sim = Simulator()
    >>> sim.isotopomers = [Isotopomer(sites=[s]) for s in [site_A, site_B]]
    >>> sim.dimensions = [dimension]

    >>> # simulate and run.
    >>> x, y = sim.run(method=one_d_spectrum)
    >>> plot(x, y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_decompose_1')

.. figure:: _images/example_decompose_1.*
    :figclass: figure-polaroid

    By default, the spectrum is an integration of the spectra from individual
    isotopomers.

Now, that we have a spectrum from two isotopomers, try setting the value of the
decompose attribute to ``True`` and observe.

.. doctest::

    >>> # set decompose to true.
    >>> sim.config.decompose = True

    >>> # simulate.
    >>> x, y = sim.run(method=one_d_spectrum)

Here, ``y`` is an ordered list of numpy arrays corresponding to the ordered
list of isotopomers. In this example, ``y`` is a list of two numpy arrays.

.. doctest::

    >>> # plot the two spectrum
    >>> plt.plot(x, y[0]) # arising from site_A # doctest:+SKIP
    >>> plt.plot(x, y[1]) # arising from site_B # doctest:+SKIP

.. .. testsetup::
..     >>> import numpy as np
..     >>> plot_save(x, np.asarray(y).T, 'example_decompose_2')

.. figure:: _images/example_decompose_2.*
    :figclass: figure-polaroid

    Spectrum from individual isotopomers when the value of the `decompose`
    config is True.

.. [#f4] Edén, M. and Levitt, M. H. Computation of orientational averages in
         solid-state nmr by gaussian spherical quadrature. J. Mag. Res.,
         **132**, *2*, 220–239, 1998. `doi:10.1006/jmre.1998.1427 <https://doi.org/10.1006/jmre.1998.1427>`_.
