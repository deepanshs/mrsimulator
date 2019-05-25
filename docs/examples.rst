

.. _examples:

========
Examples
========

Setup a :ref:`simulator` instance and load the isotopmers.

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
instance to get the list of unique isotopes from the lsit of isotopmers.

.. doctest::
    :skipif: None is None

    >>> print(sim.isotope_list)
    ['13C', '29Si', '1H']

-------------------
Generating spectrum
-------------------

Static spectrum
---------------

Generate a static spectrum with the following :ref:`spectrum` object,

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

Generate the spectrum using the :meth:`~mrsimulator.Simulator.run` method as

.. doctest::
    :skipif: None is None

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    <BLANKLINE>
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

The simulator goes through every isotopomer in the list and simulates spectrum
corresponding to the :math:`^{13}\mathrm{C}` isotopes.

.. image:: /_static/13C_static.pdf


MAS spectrum
------------

To generate a MAS spectrum, simply set the rotor frequency to the desired
value. For example,

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

and compute.

.. doctest::
    :skipif: None is None

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    <BLANKLINE>
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


Switch to another isotope
-------------------------

Generate the spectrum corresponding to another isotope by specifying
the desired isotope in the :ref:`spectrum` object as follows,

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

and compute.

.. doctest::
    :skipif: None is None

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    <BLANKLINE>
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


Seting rotor angle
------------------

Set rotor angle to any value between 0 degree and 90 degree.

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

and compute.

.. doctest::
    :skipif: None is None

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    <BLANKLINE>
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


Yet another nuclei
------------------

Swich to another nuclei. In this case :math:`^{29}\mathrm{Si}`.

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

and compute.

.. doctest::
    :skipif: None is None

    >>> freq, amp = sim.run(one_d_spectrum, verbose=1)
    <BLANKLINE>
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
