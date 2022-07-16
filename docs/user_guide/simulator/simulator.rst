.. _simulator_documentation:

=========
Simulator
=========

Overview
--------

The :ref:`simulator_api` is the top-level object in **mrsimulator**. Each :ref:`simulator_api`
object holds a list of :ref:`spin_sys_api` and :ref:`method_api` objects.
A simulator object also holds a :ref:`config_api` object, which can be modified to change aspects
of the simulation as discussed below.

Consider the following simulator set up.

.. plot::
    :context: reset

    from mrsimulator import Site, Simulator, SpinSystem
    from mrsimulator.spin_system.tensors import SymmetricTensor
    from mrsimulator.method import SpectralDimension
    from mrsimulator.method.lib import BlochDecaySpectrum

    # Setup the spin system and method objects
    Si29_site = Site(
        isotope="29Si",
        shielding_symmetric=SymmetricTensor(
            zeta=100,  # in ppm
            eta=0.2,
            alpha=1.563,  # in rads
            beta=1.2131,  # in rads
            gamma=2.132,  # in rads
        )
    )
    system = SpinSystem(sites=[Si29_site])

    method = BlochDecaySpectrum(
        channels=["29Si"],
        rotor_frequency=0,  # in Hz
        spectral_dimensions=[SpectralDimension(count=1024, spectral_width=25000)]
    )

    # Create the Simulator object
    sim = Simulator(spin_systems=[system], methods=[method])

Here, ``sim`` is a :ref:`simulator_api` object which holds one spin system and one method.
See :ref:`spin_system_documentation` and :ref:`method_documentation` documentation for more
information on the respective classes.

The following code is used to create the plots in this section.

.. plot::
    :context: close-figs

    import matplotlib.pyplot as plt
    import matplotlib as mpl

    mpl.rcParams["figure.figsize"] = (6, 3.5)
    mpl.rcParams["font.size"] = 11

    # function to render figures.
    def plot(csdm_object):
        ax = plt.subplot(projection="csdm")
        ax.plot(csdm_object.real, linewidth=1.5)
        ax.invert_xaxis()
        plt.tight_layout()
        plt.show()

.. _config_simulator:

Simulator configs
-----------------

In mrsimulator, the default configuration settings applies to a wide
range of simulations, including static, magic angle spinning (MAS), and variable angle
spinning (VAS) spectra. In certain situations, however, the default settings are insufficient
to represent the spectrum accurately.


----

Integration Volume
''''''''''''''''''

The attribute :py:attr:`~mrsimulator.simulator.ConfigSimulator.integration_volume` is an
enumeration with two string literals,
``octant`` and ``hemisphere``. The integration volume refers to the volume of a unit sphere over
which the NMR frequencies are integrated. The default value is ``"octant"``, i.e., the spectrum
comprises of integrated frequencies from the positive octant of a unit sphere.
**mrsimulator** can exploit the problem's orientational symmetry, thus optimizing the simulation
by performing a partial integration.

To learn more about the orientational symmetries, refer to Eden et al. [#f4]_

Consider the :math:`^{29}\text{Si}` site, ``Si29_site``, from the above set up. This
site has a symmetric shielding tensor with ``zeta`` and ``eta`` as 100 ppm and 0.2,
respectively. With only ``zeta`` and ``eta``, we can exploit the symmetry of the problem
and evaluate the frequency integral over the octant, which is equivalent to the
integration over the sphere. By adding the Euler angles to this tensor, we break the
symmetry, and the integration over the octant is no longer accurate.

.. skip: next

.. plot::
    :context: close-figs
    :caption: Inaccurate simulation resulting from integrating over an octant when the
        spin system contains non-zero Euler angles.

    sim.run()
    plot(sim.methods[0].simulation)

To fix this inaccurate spectrum, set the integration volume to ``"hemisphere"`` and re-simulate.

.. skip: next

.. plot::
    :context: close-figs
    :caption: Accurate CSA spectrum resulting from the frequency contributions evaluated over
        the top hemisphere.

    sim.config.integration_volume = "hemisphere"
    sim.run()
    plot(sim.methods[0].simulation)

Integration Density
'''''''''''''''''''

The attribute :py:attr:`~mrsimulator.simulator.ConfigSimulator.integration_density`
controls the number of orientations sampled over the given
volume. The resulting spectrum is an integration of the NMR resonance frequency
evaluated at these orientations. The total
number of orientations, :math:`\Theta_\text{count}`, is given as

.. math::

    \Theta_\text{count} = M (n + 1)(n + 2)/2

where :math:`M` is the number of octants and :math:`n` is value of this attribute. The
number of octants is deciphered from the value of the ``integration_volume`` attribute.
The default value of this attribute, 70, produces 2556 orientations at which the NMR
frequency contributions are evaluated.

.. skip: next

.. plot::
    :context: close-figs
    :caption: Low quality simulation from reduced integration density.

    print(sim.config.integration_density)  # default is 70
    sim.config.integration_density = 10
    sim.run()
    plot(sim.methods[0].simulation)

.. skip: next

.. plot::
    :context: close-figs
    :caption: High quality simulation from increased integration density.

    sim.config.integration_density = 100
    sim.run()
    plot(sim.methods[0].simulation)

Decreasing the integration density may decrease simulation time for computationally intensive
simulations but will also reduce the quality of the spectrum. Similarly, increasing integration
density will improve spectrum quality but at the cost of increased computation time.

Number of Sidebands
'''''''''''''''''''

The number of sidebands simulated is determined by the attribute
:py:attr:`~mrsimulator.simulator.ConfigSimulator.number_of_sidebands`. The default value is 64 and
is sufficient for most cases.

In certain circumstances, especially when the anisotropy is large or the rotor spin frequency is
low, 64 sidebands might not be sufficient.

.. skip: next

.. plot::
    :context: close-figs
    :caption: Inaccurate simulation resulting from computing low number of sidebands.

    sim.methods[0] = BlochDecaySpectrum(
        channels=["29Si"],
        rotor_frequency=200,
        spectral_dimensions=[SpectralDimension(count=1024, spectral_width=25000)],
    )
    sim.run()
    plot(sim.methods[0].simulation)

Looking at the spinning sideband patterns, we see an abrupt termination of the sideband
amplitudes at the edges. This inaccurate simulation arises from evaluating a small number of
sidebands relative to the given anisotropy. Increasing the number of sidebands to 90 should
resolve the issue.

.. skip: next

.. plot::
    :context: close-figs
    :caption: Accurate simulation after increasing number of sidebands computed.

    sim.config.number_of_sidebands = 90
    sim.run()
    plot(sim.methods[0].simulation)

Conversely, 64 sidebands might be redundant, in which case reduce the number of sidebands.
Reducing the number of sidebands will significantly improve performance, which might save
computation time when used in iterative algorithms, such as least-squares minimization.

Decompose Spectrum
''''''''''''''''''

The attribute :py:attr:`~mrsimulator.simulator.ConfigSimulator.decompose_spectrum`
is an enumeration with two string literals, ``None`` and ``spin_system``. The default value is ``None``.

If the value is ``None`` (default), the resulting simulation is a single spectrum
where the frequency contributions from all the spin systems are co-added. Consider the
following example.

.. skip: next

.. plot::
    :context: close-figs
    :caption: The frequency contributions from each individual spin systems are
        combined into one spectrum.

    # Create two distinct sites
    site_A = Site(
        isotope="1H",
        shielding_symmetric=SymmetricTensor(zeta=5, eta=0.1),
    )
    site_B = Site(
        isotope="1H",
        shielding_symmetric=SymmetricTensor(zeta=-2, eta=0.83),
    )

    # Create two single site spin systems
    sys_A = SpinSystem(sites=[site_A], name="System A")
    sys_B = SpinSystem(sites=[site_B], name="System B")

    # Create a method representing a simple 1-pulse acquire experiment
    method = BlochDecaySpectrum(
        channels=["1H"], spectral_dimensions=[SpectralDimension(count=1024, spectral_width=10000)]
    )

    # Create simulator object, simulate, and plot
    sim = Simulator(spin_systems=[sys_A, sys_B], methods=[method])
    sim.run()
    plot(sim.methods[0].simulation)

When :py:attr:`~mrsimulator.simulator.ConfigSimulator.decompose_spectrum` is set to
``spin_system``, the resulting simulation
is a series of spectra each arising from a single spin system. The number of spectra is the
same as the number of spin systems within the simulator object. Consider the same
system as above, but change the decomposition to ``spin_system``.

.. skip: next

.. plot::
    :context: close-figs
    :caption: Each spin system's frequency contributions are held in separate spectra.

    # sim already has the two spin systems and method; no need to reconstruct
    sim.config.decompose_spectrum = "spin_system"
    sim.run()
    plot(sim.methods[0].simulation)

Isotropic interpolation
'''''''''''''''''''''''

The attribute :py:attr:`~mrsimulator.simulator.ConfigSimulator.isotropic_interpolation`
is an enumeration with two string literals, ``linear`` and ``gaussian``. The default value is ``linear``.

The value specifies the interpolation scheme used in binning isotropic contributions.

Attribute Summaries
-------------------

.. cssclass:: table-bordered table-striped centered
.. _table_simulator:
.. list-table:: The attributes of a Simulator object
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - spin_systems
    - ``list``
    - An *optional* list of :ref:`spin_sys_api` objects whose spectra are simulated when running
      a simulation.

  * - methods
    - ``list``
    - An *optional* list of :ref:`method_api` objects describing the NMR methods used to simulate
      the spectra of the contained spin systems.

  * - config
    - ``dict`` or :py:class:`~mrsimulator.simulator.config.ConfigSimulator`
    - An *optional* ConfigSimulator object, or its dictionary representation, specifying how to
      run a simulation.

.. cssclass:: table-bordered table-striped centered
.. _table_sim_config:
.. list-table:: The attributes of a Simulator object
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - number_of_sidebands
    - ``int``
    - An *optional* integer greater than zero specifying the number of sidebands to simulate. The
      default is ``64`` sidebands.

  * - integration_volume
    - ``str``
    - An *optional* string specifying what portion of the unit sphere NMR frequencies are
      integrated over. The allowed strings are ``"octant"`` and ``"hemisphere"``. The default
      is ``"octant"``.

  * - integration_density
    - ``int``
    - An *optional* integer greater than zero specifying the number of orientations to sample over
      the given volume according to the equation :math:`\Theta_\text{count} = M (n + 1)(n + 2)/2`
      where :math:`M` is the number of octants. The default value is ``70``.

  * - decompose_spectrum
    - ``str``
    - An *optional* string specifying how to decompose the spectrum. The allowed strings are
      ``"none"`` and ``"spin_system"``. The value of ``"none"`` will add the spectral from
      individual spin systems into one dependent variable whereas ``"spin_system"`` will leave
      contributions from each spin system as its own dependent variable. The default is ``"none"``.

  * - isotropic_interpolation
    - ``str``
    - An *optional* string specifying an interpolation scheme when binning isotropic contributions.
      The allowed strings are ``"linear"`` and ``"gaussian"``. The default is ``"linear"``.


----

.. [#f4] Edén, M. and Levitt, M. H. Computation of orientational averages in
    solid-state nmr by gaussian spherical quadrature. J. Mag. Res.,
    **132**, *2*, 220-239, 1998. `doi:10.1006/jmre.1998.1427 <https://doi.org/10.1006/jmre.1998.1427>`_.
