

.. .. _examples:

.. Examples
.. --------

.. The following examples utilizes the `one_d_spectrum` method to simulate
.. NMR spectrum under various condition

.. - :ref:`static_example`
.. - :ref:`mas_example`
.. - :ref:`vas_example`

.. To begin, set up an instance of the :ref:`simulator_api` class and import the
.. ``one_d_spectrum`` method.

.. .. doctest::

..     >>> from mrsimulator import Simulator, Dimension
..     >>> from mrsimulator.methods import one_d_spectrum
..     >>> sim = Simulator()

.. Let's add the isotopomers to the ``sim`` object from a JSON serialized
.. isotopomers file, see :ref:`load_isotopomers` for details,

.. .. doctest::

..     >>> filename = 'https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers_test.json'

..     >>> sim.load_isotopomers(filename)
..     Downloading '/DeepanshS/mrsimulator-test/master/isotopomers_test.json'
..     from 'raw.githubusercontent.com' to file 'isotopomers_test.json'.
..     [████████████████████████████████████]

.. .. testsetup::

..     >>> import os
..     >>> os.remove('isotopomers_test.json')

.. Once the isotopomers are assigned, you may use the
.. :meth:`~mrsimulator.Simulator.get_isotopes` method to list the
.. unique isotope symbols present in the list of isotopomers.

.. .. doctest::

..     >>> print(sim.get_isotopes(spin=0.5)) # doctest: +SKIP
..     ['1H', '13C', '29Si']

.. In this example, the list of isotopomers contain three unique isotopes,
.. :math:`^{13}\mathrm{C}`, :math:`^{29}\mathrm{Si}`, and :math:`^{1}\mathrm{H}`.


.. .. _static_example:

.. Static spectrum
.. '''''''''''''''

.. To simulate a static spectrum, set up with the following
.. :ref:`dimension_api` object,

.. .. doctest::

..     >>> dim = dict(isotope="13C", magnetic_flux_density="9.4 T", rotor_frequency="0 kHz",
..     ...         rotor_angle="54.735 deg", number_of_points=8192, spectral_width="5 kHz",
..     ...         reference_offset="0 Hz")
..     >>> sim.dimensions = [Dimension.parse_dict_with_units(dim)]

.. The above spectroscopic dimension will simulate a :math:`^{13}\mathrm{C}`
.. static spectrum at 9.4 T magnetic field over 5 kHz frequency-bandwidth using
.. 8192 points.

.. Now, generate the line-shape using the :meth:`~mrsimulator.Simulator.run`
.. method as

.. .. doctest::

..     >>> freq, amp = sim.run(one_d_spectrum)

.. The simulator object goes through every isotopomer in the list and
.. simulates the line-shape corresponding to the :math:`^{13}\mathrm{C}` isotopes.
.. In this example, there are two isotopomers with :math:`^{13}\mathrm{C}` sites.

.. You may visualize the spectrum using any plotting library of choice. We use
.. matplotlib in our examples.

..     >>> import matplotlib.pyplot as plt
..     >>> def plot(x, y):
..     ...     plt.figure(figsize=(4, 3))
..     ...     plt.plot(x, y)
..     ...     plt.xlabel(f"frequency ratio / {str(x.unit)}")
..     ...     plt.xlim([x.value.max(), x.value.min()]) # for reverse axis
..     ...     plt.grid()
..     ...     plt.tight_layout()
..     ...     plt.show()

..     >>> plot(freq, amp) # doctest: +SKIP

.. .. .. testsetup::

.. ..    >>> plot_save(freq, amp, '13C_static.pdf')

.. .. figure:: _images/13C_static.*
..     :figclass: figure-polaroid

.. .. _mas_example:

.. Magic angle spinning (MAS) spectrum
.. '''''''''''''''''''''''''''''''''''

.. To simulate a magic angle spinning (MAS) spectrum, set the rotor frequency
.. of the spectroscopic dimension object to the desired value, and set the
.. rotor_angle to :math:`54.735^\circ`. In the following example,
.. the spectroscopic dimension object is similar to the one from the
.. previous examples, except for the value of the ``rotor_frequency`` which
.. is set to 100 Hz.

.. .. doctest::

..     >>> dim = dict(isotope="13C", magnetic_flux_density="9.4 T", rotor_frequency="100 Hz",
..     ...         rotor_angle="54.735 deg", number_of_points=8192, spectral_width="5 kHz",
..     ...         reference_offset="0 Hz")
..     >>> sim.dimensions = [Dimension.parse_dict_with_units(dim)]

.. Now compute the line-shape as before.

.. .. doctest::

..     >>> freq, amp = sim.run(one_d_spectrum)
..     >>> plot(freq, amp) # doctest: +SKIP

.. .. .. testsetup::

.. ..    >>> plot_save(freq, amp, '13C_mas_100Hz.pdf')

.. .. figure:: _images/13C_mas_100Hz.*
..     :figclass: figure-polaroid


.. .. _vas_example:

.. Variable angle spinning (VAS) spectrum
.. ''''''''''''''''''''''''''''''''''''''

.. To simulate a variable angle spinning spectrum, set the rotor angle
.. and the rotor frequency to the desired value. The rotor angle may be
.. set to any value from :math:`0^\circ` to :math:`90^\circ`.
.. In the following example, the spectroscopic dimension object is similar
.. to spectroscopic dimension from the previous example, except the
.. ``rotor_angle`` and ``rotor_frequency`` is set to :math:`90^\circ` and
.. and 500 Hz, respectively.

.. .. doctest::

..     >>> dim = dict(isotope="13C", magnetic_flux_density="9.4 T", rotor_frequency="500 Hz",
..     ...         rotor_angle="90 deg", number_of_points=8192, spectral_width="5 kHz",
..     ...         reference_offset="0 Hz")
..     >>> sim.dimensions = [Dimension.parse_dict_with_units(dim)]

.. The simulated lineshape.

.. .. doctest::

..     >>> freq, amp = sim.run(one_d_spectrum)
..     >>> plot(freq, amp) # doctest: +SKIP

.. .. .. testsetup::

.. ..    >>> plot_save(freq, amp, '13C_vas_100Hz_90.pdf')

.. .. figure:: _images/13C_vas_100Hz_90.*
..     :figclass: figure-polaroid


.. Switching to a different isotope
.. ''''''''''''''''''''''''''''''''

.. Up till now, we were simulating a one-dimensional :math:`^{13}\mathrm{C}`
.. spectrum under conditions. Notice, however, there are three unique isotopes,
.. :math:`^{13}\mathrm{C}`, :math:`^{29}\mathrm{Si}`, and :math:`^{1}\mathrm{H}`,
.. in the list of isotopomers.
.. To simulate, for example, a :math:`^{29}\mathrm{Si}` dimensions, create a new
.. spectroscopic dimension with "29Si" as the value of the ``isotope`` key.

.. .. doctest::

..     >>> dim = dict(isotope="29Si", magnetic_flux_density="9.4 T", rotor_frequency="1 kHz",
..     ...         rotor_angle="54.735 deg", number_of_points=8192, spectral_width="30 kHz",
..     ...         reference_offset="-5 kHz")
..     >>> sim.dimensions = [Dimension.parse_dict_with_units(dim)]

.. Run the simulation.

.. .. doctest::

..     >>> freq, amp = sim.run(one_d_spectrum)
..     >>> plot(freq, amp) # doctest: +SKIP

.. .. .. testsetup::

.. ..    >>> plot_save(freq, amp, '29Si_mas_1kHz.pdf')

.. .. figure:: _images/29Si_mas_1kHz.*
..     :figclass: figure-polaroid


.. In this another examples, we simulate a :math:`^1\mathrm{H}` dimensions.

.. .. doctest::

..     >>> dim = dict(isotope="1H", magnetic_flux_density="9.4 T", rotor_frequency="2 kHz",
..     ...         rotor_angle="54.735 deg", number_of_points=8192, spectral_width="50 kHz",
..     ...         reference_offset="0 Hz")
..     >>> sim.dimensions = [Dimension.parse_dict_with_units(dim)]

.. The line-shape simulation

.. .. doctest::

..     >>> freq, amp = sim.run(one_d_spectrum)
..     >>> plot(freq, amp) # doctest: +SKIP

.. .. .. testsetup::
.. ..    >>> plot_save(freq, amp, '1H_mas_2kHz.pdf')

.. .. figure:: _images/1H_mas_2kHz.*
..     :figclass: figure-polaroid
