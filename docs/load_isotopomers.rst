

.. _load_isotopomers:

==========================
Setting up the isotopomers
==========================

The isotopomers may either be specified as a list of :ref:`isotopomer`
objects or directly imported from a JSON serialized isotopomers file.

**Using python list:**

.. doctest::

    >>> isotopomers = [
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "13C",
    ...                 "isotropic_chemical_shift": "1 Hz",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "-3.89 kHz",
    ...                     "asymmetry": 0.25
    ...                 }
    ...             }
    ...         ],
    ...         "abundance": "12.1 %"
    ...     },
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "1H",
    ...                 "isotropic_chemical_shift": "1 kHz",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "8.2 kHz",
    ...                     "asymmetry": 0.0
    ...                 }
    ...             }
    ...         ]
    ...     },
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "1H",
    ...                 "isotropic_chemical_shift": "1 kHz",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "8.2 kHz",
    ...                     "asymmetry": 0.0
    ...                 }
    ...             }
    ...         ]
    ...     },
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "1H",
    ...                 "isotropic_chemical_shift": "3 kHz",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "23.2 kHz",
    ...                     "asymmetry": 0.0
    ...                 }
    ...             }
    ...         ]
    ...     },
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "29Si",
    ...                 "isotropic_chemical_shift": "1.64 kHz",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "7.36 kHz",
    ...                     "asymmetry": 0.0
    ...                 }
    ...             }
    ...         ],
    ...         "abundance": "4.67 %"
    ...     },
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "29Si",
    ...                 "isotropic_chemical_shift": "43 kHz",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "8.36 kHz",
    ...                     "asymmetry": 0.5
    ...                 }
    ...             }
    ...         ],
    ...         "abundance": "4.67 %"
    ...     },
    ... ]

Create an instance of the :ref:`simulator_api` class,

.. doctest::

    >>> from mrsimulator import Simulator

and then assign the isotopomers using either

.. doctest::

    >>> sim1 = Simulator(isotopomers)

or

.. doctest::

    >>> sim1 = Simulator()
    >>> sim1.isotopomers = isotopomers



**Import from JSON serialized isotopomers file**


.. doctest::
    :skipif: None is None

    >>> from pprint import pprint
    >>> filename = 'https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers.json'
    >>> st2 = Simulator()
    >>> st2.load_isotopomers(filename)
    Downloading '/DeepanshS/mrsimulator-test/master/isotopomers.json' from 'raw.githubusercontent.com' to file 'isotopomers.json'.
    [█████████████████████████████████████████████████████████████████████████]

.. doctest::
    :skipif: None is None

    >>> pprint(st2.isotopomers)
    [{'abundance': '12%',
      'sites': [{'isotope_symbol': '13C',
                 'isotropic_chemical_shift': '1 Hz',
                 'shielding_symmetric': {'anisotropy': '-3.89 kHz',
                                         'asymmetry': 0.25}}]},
     {'sites': [{'isotope_symbol': '13C',
                 'isotropic_chemical_shift': '1 kHz',
                 'shielding_symmetric': {'anisotropy': '8.2 kHz',
                                         'asymmetry': 0.0}}]},
     {'sites': [{'isotope_symbol': '1H',
                 'isotropic_chemical_shift': '3 kHz',
                 'shielding_symmetric': {'anisotropy': '23.2 kHz',
                                         'asymmetry': 0.0}}]},
     {'sites': [{'isotope_symbol': '29Si',
                 'isotropic_chemical_shift': '1.64 kHz',
                 'shielding_symmetric': {'anisotropy': '7.36 kHz',
                                         'asymmetry': 0.0}}]},
     {'sites': [{'isotope_symbol': '29Si',
                 'isotropic_chemical_shift': '43 kHz',
                 'shielding_symmetric': {'anisotropy': '8.36 kHz',
                                         'asymmetry': 0.5}}]},
     {'sites': [{'isotope_symbol': '29Si',
                 'isotropic_chemical_shift': '10 kHz',
                 'shielding_symmetric': {'anisotropy': '6.36 kHz',
                                         'asymmetry': 0.0}}]},
     {'sites': [{'isotope_symbol': '1H',
                 'isotropic_chemical_shift': '5.6 kHz',
                 'shielding_symmetric': {'anisotropy': '13.2 kHz',
                                         'asymmetry': 0.0}}]}]
