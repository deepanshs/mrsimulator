

.. _examples:

========
Examples
========

Setup an instance of the :ref:`simulator_api` class and import the
``one_d_spectrum`` method.

.. doctest::

    >>> from mrsimulator import Simulator
    >>> from mrsimulator.methods import one_d_spectrum
    >>> sim = Simulator()

Now add the isotopomers to the ``sim`` object.
See :ref:`load_isotopomers` for details. Here, we import the list of
isotopomers from a JSON serialized isotopomers file.

.. doctest::

    >>> filename = 'https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers_ppm.json'
    >>> sim.load_isotopomers(filename)
    Downloading '/DeepanshS/mrsimulator-test/master/isotopomers_ppm.json' from 'raw.githubusercontent.com' to file 'isotopomers.json'.
    [█████████████████████████████████████████████████████████████████████████]

.. testcleanup::

    import os
    os.remove('isotopomers_ppm.json')

Use the :attr:`~mrsimulator.Simulator.isotope_list` attribute of the ``sim``
instance to get the list of unique isotopes from the list of isotopomers.

.. doctest::

    >>> print(sim.isotope_list)
    ['1H', '13C', '29Si']

In this example, the list of isotopomers contains three unique isotopes,
:math:`^{13}\mathrm{C}`, :math:`^{29}\mathrm{Si}`, and :math:`^{1}\mathrm{H}`.

---------------
Static spectrum
---------------

To generate a static spectrum, set up with the following :ref:`spectrum`
object,

.. doctest::

    >>> sim.spectrum = {
    ...     "direct_dimension": {
    ...         "nucleus": "13C",
    ...         "magnetic_flux_density": "9.4 T",
    ...         "rotor_frequency": "0 kHz",
    ...         "rotor_angle": "54.735 deg",
    ...         "number_of_points": 8192,
    ...         "spectral_width": "5 kHz",
    ...         "reference_offset": "0 Hz",
    ...     }
    ... }

The above spectrum object is set to simulate a :math:`^{13}\mathrm{C}` static
spectrum at 9.4 T magnetic field over 5 kHz frequency-bandwidth using 8192
points.

Now, generate the lineshape with the :meth:`~mrsimulator.Simulator.run`
method as

.. doctest::

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T.
    Setting rotation angle to 0.9553059660790962 rad.
    Setting rotation frequency to 0.0 Hz.
    Detecting 13C(I=0.5, precession frequency = 100.65896 MHz) isotope.
    Recording 13C spectrum with 8192 points over a 5000.0 Hz bandwidth and a reference offset of 0.0 Hz.
    <BLANKLINE>
    13C site 0 from isotopomer 0 @ 100.0% abundance
    -----------------------------------------------
    Isotropic chemical shift = 1.0 ppm
    Shielding anisotropy = -3.89 ppm
    Shielding asymmetry = 0.25
    <BLANKLINE>
    13C site 0 from isotopomer 1 @ 100.0% abundance
    -----------------------------------------------
    Isotropic chemical shift = 1.0 ppm
    Shielding anisotropy = 8.2 ppm
    Shielding asymmetry = 0.0

The simulator goes through every isotopomer in the list and simulates the
lineshape corresponding to the :math:`^{13}\mathrm{C}` isotopes. In this
example, there are two :math:`^{13}\mathrm{C}` isotopes in the list of the
isotopomers.

You may visualize the spectrum using any plotting library of choise. We use
matplotlib in our examples.

.. doctest::

    >>> import matplotlib.pyplot as plt
    >>> def plot(x, y):
    ...     plt.plot(x,y)
    ...     plt.xlabel(f'frequency / {x.unit}')
    ...     plt.show()

    >>> plot(freq, amp)

.. image:: /_static/13C_static.png


-----------------------------------
Magic angle spinning (MAS) spectrum
-----------------------------------

To generate a magic angle spinning (MAS) spectrum, set the rotor frequency
of the :ref:`spectrum` object to the desired value. In the following example,
the spectrum object is similar to the one from the previous examples, except
for the ``rotor_frequency`` which is set to 100 Hz.

.. doctest::

    >>> sim.spectrum = {
    ...     "direct_dimension": {
    ...         "nucleus": "13C",
    ...         "magnetic_flux_density": "9.4 T",
    ...         "rotor_frequency": "100 Hz",
    ...         "rotor_angle": "54.735 deg",
    ...         "number_of_points": 8192,
    ...         "spectral_width": "5 kHz",
    ...         "reference_offset": "0 Hz",
    ...     }
    ... }

Now compute the lineshape as before.

.. doctest::

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T.
    Setting rotation angle to 0.9553059660790962 rad.
    Setting rotation frequency to 100.0 Hz.
    Detecting 13C(I=0.5, precession frequency = 100.65896 MHz) isotope.
    Recording 13C spectrum with 8192 points over a 5000.0 Hz bandwidth and a reference offset of 0.0 Hz.
    <BLANKLINE>
    13C site 0 from isotopomer 0 @ 100.0% abundance
    -----------------------------------------------
    Isotropic chemical shift = 1.0 ppm
    Shielding anisotropy = -3.89 ppm
    Shielding asymmetry = 0.25
    <BLANKLINE>
    13C site 0 from isotopomer 1 @ 100.0% abundance
    -----------------------------------------------
    Isotropic chemical shift = 1.0 ppm
    Shielding anisotropy = 8.2 ppm
    Shielding asymmetry = 0.0

.. doctest::

    >>> plot(freq, amp)

.. image:: /_static/13C_mas_1kHz.png


-----------------------------
Switch to a different isotope
-----------------------------

Generate a new :ref:`spectrum` object with a different isotope. The isotope
is specified with the `nucleus` key, as shown below. In the following
example, a :math:`^1\mathrm{H}` spectrum is simulated at 9.4 T field, spinning
at the magic angle at 2 kHz frequency, and sampled over 100 kHz frequency
bandwidth with 8192 points.

.. doctest::

    >>> sim.spectrum = {
    ...     "direct_dimension": {
    ...         "nucleus": "1H",
    ...         "magnetic_flux_density": "9.4 T",
    ...         "rotor_frequency": "2 kHz",
    ...         "rotor_angle": "54.735 deg",
    ...         "number_of_points": 8192,
    ...         "spectral_width": "50 kHz",
    ...         "reference_offset": "0 Hz",
    ...     }
    ... }

Now compute the lineshape.

.. doctest::

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T.
    Setting rotation angle to 0.9553059660790962 rad.
    Setting rotation frequency to 2000.0 Hz.
    Detecting 1H(I=0.5, precession frequency = 400.228301848 MHz) isotope.
    Recording 1H spectrum with 8192 points over a 50000.0 Hz bandwidth and a reference offset of 0.0 Hz.
    <BLANKLINE>
    1H site 0 from isotopomer 2 @ 100.0% abundance
    ----------------------------------------------
    Isotropic chemical shift = 3.0 ppm
    Shielding anisotropy = 23.2 ppm
    Shielding asymmetry = 0.0
    <BLANKLINE>
    1H site 0 from isotopomer 6 @ 100.0% abundance
    ----------------------------------------------
    Isotropic chemical shift = 5.6 ppm
    Shielding anisotropy = 13.2 ppm
    Shielding asymmetry = 0.0

.. doctest::

    >>> plot(freq, amp)

.. image:: /_static/1H_mas_2kHz.png


In this example, we simulate the MAS lineshape of :math:`^{29}\mathrm{Si}`
isotope.

.. doctest::

    >>> sim.spectrum = {
    ...     "direct_dimension": {
    ...         "nucleus": "29Si",
    ...         "magnetic_flux_density": "9.4 T",
    ...         "rotor_frequency": "1 kHz",
    ...         "rotor_angle": "54.735 deg",
    ...         "number_of_points": 8192,
    ...         "spectral_width": "30 kHz",
    ...         "reference_offset": "5 kHz",
    ...     }
    ... }

The simulated lineshape.

.. doctest::

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T.
    Setting rotation angle to 0.9553059660790962 rad.
    Setting rotation frequency to 1000.0 Hz.
    Detecting 29Si(I=0.5, precession frequency = -79.571 MHz) isotope.
    Recording 29Si spectrum with 8192 points over a 30000.0 Hz bandwidth and a reference offset of 5000.0 Hz.
    <BLANKLINE>
    29Si site 0 from isotopomer 3 @ 100.0% abundance
    ------------------------------------------------
    Isotropic chemical shift = -100.0 ppm
    Shielding anisotropy = 1.36 ppm
    Shielding asymmetry = 0.0
    <BLANKLINE>
    29Si site 0 from isotopomer 4 @ 100.0% abundance
    ------------------------------------------------
    Isotropic chemical shift = -100.0 ppm
    Shielding anisotropy = 70.36 ppm
    Shielding asymmetry = 0.0
    <BLANKLINE>
    29Si site 0 from isotopomer 5 @ 100.0% abundance
    ------------------------------------------------
    Isotropic chemical shift = -90.0 ppm
    Shielding anisotropy = 80.36 ppm
    Shielding asymmetry = 0.5

.. doctest::

    >>> plot(freq, amp)

.. image:: /_static/29Si_mas_1kHz.png


-----------------------
Variable angle spinning
-----------------------

The rotor angle may be set to any value from :math:`0^\circ` to
:math:`90^\circ`. In the following example, the :ref:`spectrum`
object is the same are from the previous example, except the
``rotor_angle`` is set to :math:`90^\circ`.

.. doctest::

    >>> sim.spectrum = {
    ...     "direct_dimension": {
    ...         "nucleus": "1H",
    ...         "magnetic_flux_density": "9.4 T",
    ...         "rotor_frequency": "2 kHz",
    ...         "rotor_angle": "90 deg",
    ...         "number_of_points": 8192,
    ...         "spectral_width": "50 kHz",
    ...         "reference_offset": "0 Hz",
    ...     }
    ... }

The simulated lineshape.

.. doctest::

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T.
    Setting rotation angle to 1.5707963267948966 rad.
    Setting rotation frequency to 2000.0 Hz.
    Detecting 1H(I=0.5, precession frequency = 400.228301848 MHz) isotope.
    Recording 1H spectrum with 8192 points over a 50000.0 Hz bandwidth and a reference offset of 0.0 Hz.
    <BLANKLINE>
    1H site 0 from isotopomer 2 @ 100.0% abundance
    ----------------------------------------------
    Isotropic chemical shift = 3.0 ppm
    Shielding anisotropy = 23.2 ppm
    Shielding asymmetry = 0.0
    <BLANKLINE>
    1H site 0 from isotopomer 6 @ 100.0% abundance
    ----------------------------------------------
    Isotropic chemical shift = 5.6 ppm
    Shielding anisotropy = 13.2 ppm
    Shielding asymmetry = 0.0

.. doctest::

    >>> plot(freq, amp)

.. image:: /_static/1H_mas_2khz_90deg.png
