

.. _examples:

========
Examples
========

Setup a :ref:`simulator_api` instance and load the isotopomers.
See :ref:`load_isotopomers` for details.

.. doctest::
    :skipif: None is None

    >>> from mrsimulator import Simulator
    >>> from mrsimulator.methods import one_d_spectrum

    >>> filename = 'https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers.json'
    >>> sim = Simulator()
    >>> sim.load_isotopomers(filename)
    Downloading '/DeepanshS/mrsimulator-test/master/isotopomers.json' from 'raw.githubusercontent.com' to file 'isotopomers.json'.
    [█████████████████████████████████████████████████████████████████████████]

Use the :attr:`~mrsimulator.Simulator.isotope_list` attribute of the ``sim``
instance to get the list of unique isotopes from the list of isotopomers.

.. doctest::
    :skipif: None is None

    >>> print(sim.isotope_list)
    ['13C', '29Si', '1H']

In this example, the list of isotopomers contains three unique isotopes,
:math:`^{13}\mathrm{C}`, :math:`^{29}\mathrm{Si}`, and :math:`^{1}\mathrm{H}`.

---------------
Static spectrum
---------------

To generate a static spectrum, set up with the following :ref:`spectrum`
object,

.. doctest::
    :skipif: None is None

    >>> sim.spectrum = {
    ...     "direct_dimension": {
    ...         "nucleus": "13C",
    ...         "magnetic_flux_density": "9.4 T",
    ...         "rotor_frequency": "0 kHz",
    ...         "rotor_angle": "54.735 deg",
    ...         "number_of_points": 8192,
    ...         "spectral_width": "30 kHz",
    ...         "reference_offset": "0 Hz",
    ...     }
    ... }

The above spectrum object is set to simulate a :math:`^{13}\mathrm{C}` static
spectrum at 9.4 T magnetic field over 30 kHz frequency-bandwidth using 8192
points.

Now, generate the lineshape with the :meth:`~mrsimulator.Simulator.run`
method as

.. doctest::
    :skipif: None is None

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T
    Setting rotation angle to 0.9553059660790962 rad
    Setting rotation frequency to 0.0 Hz
    Detecting 13C(I=0.5, precession frequency = 100.65896 MHz) isotope
    Recording 13C spectrum with 8192 points over a 30000.0 Hz bandwidth and a reference offset of 0.0 Hz.
    <BLANKLINE>
    13C site 0 in isotopomer 0 @ 12.0% abundance
    --------------------------------------------
    isotropic chemical shift = 1.0 Hz
    chemical shift anisotropy = -3890.0 Hz
    chemical shift asymmetry = 0.25
    <BLANKLINE>
    13C site 0 in isotopomer 1 @ 100.0% abundance
    ---------------------------------------------
    isotropic chemical shift = 1000.0 Hz
    chemical shift anisotropy = 8200.0 Hz
    chemical shift asymmetry = 0.0
    <BLANKLINE>
    Execution time 0.03153 s

The simulator goes through every isotopomer in the list and simulates the
lineshape corresponding to the :math:`^{13}\mathrm{C}` isotopes. In this
example, there are two :math:`^{13}\mathrm{C}` isotopes in the list of the
isotopomers.

.. image:: /_static/13C_static.pdf


------------
MAS spectrum
------------

To generate a MAS spectrum, set the rotor frequency of the :ref:`spectrum`
object to the desired value. For example,

.. doctest::
    :skipif: None is None

    >>> sim.spectrum = {
    ...     "direct_dimension": {
    ...         "nucleus": "13C",
    ...         "magnetic_flux_density": "9.4 T",
    ...         "rotor_frequency": "1 kHz",
    ...         "rotor_angle": "54.735 deg",
    ...         "number_of_points": 8192,
    ...         "spectral_width": "30 kHz",
    ...         "reference_offset": "0 Hz",
    ...     }
    ... }

Compute the lineshape.

.. doctest::
    :skipif: None is None

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T
    Setting rotation angle to 0.9553059660790962 rad
    Setting rotation frequency to 1000.0 Hz
    Detecting 13C(I=0.5, precession frequency = 100.65896 MHz) isotope
    Recording 13C spectrum with 8192 points over a 30000.0 Hz bandwidth and a reference offset of 0.0 Hz.
    <BLANKLINE>
    13C site 0 in isotopomer 0 @ 12.0% abundance
    --------------------------------------------
    isotropic chemical shift = 1.0 Hz
    chemical shift anisotropy = -3890.0 Hz
    chemical shift asymmetry = 0.25
    <BLANKLINE>
    13C site 0 in isotopomer 1 @ 100.0% abundance
    ---------------------------------------------
    isotropic chemical shift = 1000.0 Hz
    chemical shift anisotropy = 8200.0 Hz
    chemical shift asymmetry = 0.0
    <BLANKLINE>
    Execution time 0.027652 s

.. image:: /_static/13C_mas_1kHz.pdf


-------------------------
Switch to another isotope
-------------------------

Generate a new :ref:`spectrum` object with a different isotope. The isotope
may be specified with the `nucleus` key, as shown below. In the following
example, a :math:`^1\mathrm{H}` spectrum is simulated at 9.4 T field, spinning
at the magic angle at 2 kHz frequency, and sampled over 100 kHz frequency
bandwidth with 8192 points.

.. doctest::
    :skipif: None is None

    >>> sim.spectrum = {
    ...     "direct_dimension": {
    ...         "nucleus": "1H",
    ...         "magnetic_flux_density": "9.4 T",
    ...         "rotor_frequency": "2 kHz",
    ...         "rotor_angle": "54.735 deg",
    ...         "number_of_points": 8192,
    ...         "spectral_width": "100 kHz",
    ...         "reference_offset": "0 Hz",
    ...     }
    ... }

Now compute the lineshape.

.. doctest::
    :skipif: None is None

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T
    Setting rotation angle to 0.9553059660790962 rad
    Setting rotation frequency to 2000.0 Hz
    Detecting 1H(I=0.5, precession frequency = 400.228301848 MHz) isotope
    Recording 1H spectrum with 8192 points over a 100000.0 Hz bandwidth and a reference offset of 0.0 Hz.
    <BLANKLINE>
    1H site 0 in isotopomer 2 @ 100.0% abundance
    --------------------------------------------
    isotropic chemical shift = 3000.0 Hz
    chemical shift anisotropy = 23200.0 Hz
    chemical shift asymmetry = 0.0
    <BLANKLINE>
    1H site 0 in isotopomer 6 @ 100.0% abundance
    --------------------------------------------
    isotropic chemical shift = 5600.0 Hz
    chemical shift anisotropy = 13200.0 Hz
    chemical shift asymmetry = 0.0
    <BLANKLINE>
    Execution time 0.03747 s

.. image:: /_static/1H_mas_2kHz.pdf


-------------------
Setting rotor angle
-------------------

The rotor angle may be set to any value from :math:`0^\circ` to
:math:`90^\circ`. In the following example, the :ref:`spectrum`
object is the same are from the previous example, except the
``rotor_angle`` is set to :math:`90^\circ`.

.. doctest::
    :skipif: None is None

    >>> sim.spectrum = {
    ...     "direct_dimension": {
    ...         "nucleus": "1H",
    ...         "magnetic_flux_density": "9.4 T",
    ...         "rotor_frequency": "2 kHz",
    ...         "rotor_angle": "90 deg",
    ...         "number_of_points": 8192,
    ...         "spectral_width": "100 kHz",
    ...         "reference_offset": "0 Hz",
    ...     }
    ... }

The simulated lineshape.

.. doctest::
    :skipif: None is None

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T
    Setting rotation angle to 1.5707963267948966 rad
    Setting rotation frequency to 2000.0 Hz
    Detecting 1H(I=0.5, precession frequency = 400.228301848 MHz) isotope
    Recording 1H spectrum with 8192 points over a 100000.0 Hz bandwidth and a reference offset of 0.0 Hz.
    <BLANKLINE>
    1H site 0 in isotopomer 2 @ 100.0% abundance
    --------------------------------------------
    isotropic chemical shift = 3000.0 Hz
    chemical shift anisotropy = 23200.0 Hz
    chemical shift asymmetry = 0.0
    <BLANKLINE>
    1H site 0 in isotopomer 6 @ 100.0% abundance
    --------------------------------------------
    isotropic chemical shift = 5600.0 Hz
    chemical shift anisotropy = 13200.0 Hz
    chemical shift asymmetry = 0.0
    <BLANKLINE>
    Execution time 0.050539 s

.. image:: /_static/1H_mas_2khz_90deg.pdf


-------------------------
Switch to another nucleus
-------------------------

In this example, we simulate the MAS lineshape of :math:`^{29}\mathrm{Si}`
isotope.

.. doctest::
    :skipif: None is None

    >>> sim.spectrum = {
    ...     "direct_dimension": {
    ...         "nucleus": "29Si",
    ...         "magnetic_flux_density": "9.4 T",
    ...         "rotor_frequency": "2 kHz",
    ...         "rotor_angle": "90 deg",
    ...         "number_of_points": 8192,
    ...         "spectral_width": "150 kHz",
    ...         "reference_offset": "20 Hz",
    ...     }
    ... }

The simulated lineshape.

.. doctest::
    :skipif: None is None

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T
    Setting rotation angle to 0.9553059660790962 rad
    Setting rotation frequency to 2000.0 Hz
    Detecting 29Si(I=0.5, precession frequency = -79.571 MHz) isotope
    Recording 29Si spectrum with 8192 points over a 150000.0 Hz bandwidth and a reference offset of 20000.0 Hz.
    <BLANKLINE>
    29Si site 0 in isotopomer 3 @ 100.0% abundance
    ----------------------------------------------
    isotropic chemical shift = 1640.0 Hz
    chemical shift anisotropy = 7360.0 Hz
    chemical shift asymmetry = 0.0
    <BLANKLINE>
    29Si site 0 in isotopomer 4 @ 100.0% abundance
    ----------------------------------------------
    isotropic chemical shift = 43000.0 Hz
    chemical shift anisotropy = 8360.0 Hz
    chemical shift asymmetry = 0.5
    <BLANKLINE>
    29Si site 0 in isotopomer 5 @ 100.0% abundance
    ----------------------------------------------
    isotropic chemical shift = 10000.0 Hz
    chemical shift anisotropy = 6360.0 Hz
    chemical shift asymmetry = 0.0
    <BLANKLINE>
    Execution time 0.045870999999999995

.. image:: /_static/29Si_mas_2kHz.pdf
