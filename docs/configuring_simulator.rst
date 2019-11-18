

.. _config_simulator:

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

============================
Configuring Simulator object
============================

Up till now, we have been using the simulator object with the default setting.
We choose the default settings such that it applies to a wide variety of
simulations including, static, magic angle spinning (MAS), and variable angle
spinning (VAS) lineshapes. In certain situations, however, the default settings
are not sufficient to accurately represent the lineshape. In such cases, users
are advised to modify these settings as required. In the following section,
we will walk through these configuration settings.

The :ref:`simulator_api` class is configured using the
:attr:`~mrsimulator.simulator.Simulator.config` attribute. The default value
of the config attributes is as follows,

.. doctest::

    >>> from mrsimulator import Simulator, Isotopomer, Dimension, Site
    >>> from mrsimulator.methods import one_d_spectrum
    >>> the_simulator = Simulator()

    >>> the_simulator.config
    ConfigSimulator(number_of_sidebands=64, integration_volume=octant, integration_density=70, decompose=False)

Number of sidebands
-------------------
The value of the attribute ``number_of_sidebands`` is the number of sidebands
requested in evaluated the lineshapes. The default value 64 is sufficient for
most cases, as seen from our previous examples. In certain cases, especially
when the anisotropy is large or the rotor spin frequency is low, 64 sidebands
might not be sufficient, and the user is advised to increase this number as
required. Consider the following example,

.. doctest::

    >>> simulator_1 = Simulator()

    >>> # create a site with a large anisotropy, 100 ppm.
    >>> site = Site(isotope='29Si', shielding_symmetric={'zeta': 100, 'eta': 0.2})

    >>> # create a dimension with a low rotor frequency, 200 Hz
    >>> dimension = Dimension(isotope='29Si', spectral_width=25000, rotor_frequency=200)

    >>> simulator_1.isotopomers = [Isotopomer(sites=[site])]
    >>> simulator_1.dimensions = [dimension]

    >>> # simulate and plot
    >>> x, y = simulator_1.run(method=one_d_spectrum)
    >>> plt.plot(x, y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_sidebands_1')

.. figure:: _images/example_sidebands_1.*
    :figclass: figure-polaroid

If you are familiar with NMR lineshapes, you may notice that the above sideband
simulation is incorrect. This inaccuracy is predominant from the abrupt
termination of the sideband amplitudes at the edges. As mentioned earlier, this
inaccuracy arises from evaluating a relatively small number of sidebands for
the given anisotropy. Let's increase the number of sidebands to ``90`` and
observe.

.. doctest::

    >>> simulator_1.config.number_of_sidebands = 90
    >>> x, y = simulator_1.run(method=one_d_spectrum)
    >>> plt.plot(x, y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_sidebands_2')

.. figure:: _images/example_sidebands_2.*
    :figclass: figure-polaroid

.. If you are unsure on how many sideband to use in simulation, you may run a
.. convergence test
.. You will know if lineshape has converged if a further increase in the number of
.. sidebands has minimal effect on the lineshape. In this case, ``90`` sidebands
.. are sufficient to describe the spectrum as we see no significant changes to the
.. lineshape when the number of sidebands is increased further.

.. .. doctest::

..     >>> simulator_1.config.number_of_sidebands = 128
..     >>> x, y = simulator_1.run(method=one_d_spectrum)
..     >>> plt.plot(x, y) # doctest:+SKIP

.. .. testsetup::

..     >>> plot_save(x, y, 'example_sidebands_3')

.. .. figure:: _images/example_sidebands_3.*
..     :figclass: figure-polaroid
