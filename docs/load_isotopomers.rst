

.. _load_isotopomers:

==========================
Setting up the isotopomers
==========================

The isotopomers may either be specified as a list of :ref:`isotopomer`
objects or directly imported from a JSON serialized isotopomers file.

**Using python list of** :ref:`isotopomer` **objects:**

.. doctest::

    >>> isotopomers = [
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "13C",
    ...                 "isotropic_chemical_shift": "1 ppm",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "-3.89 ppm",
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
    ...                 "isotropic_chemical_shift": "1 ppm",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "8.2 ppm",
    ...                     "asymmetry": 0.0
    ...                 }
    ...             }
    ...         ]
    ...     },
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "1H",
    ...                 "isotropic_chemical_shift": "1 ppm",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "8.2 ppm",
    ...                     "asymmetry": 0.0
    ...                 }
    ...             }
    ...         ]
    ...     },
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "1H",
    ...                 "isotropic_chemical_shift": "3 ppm",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "23.2 ppm",
    ...                     "asymmetry": 0.0
    ...                 }
    ...             }
    ...         ]
    ...     },
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "29Si",
    ...                 "isotropic_chemical_shift": "-90 ppm",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "1 mHz/Hz",
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
    ...                 "isotropic_chemical_shift": "-100 ppm",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "80.36 µHz/Hz",
    ...                     "asymmetry": 0.0
    ...                 }
    ...             }
    ...         ],
    ...         "abundance": "4.67 %"
    ...     },
    ... ]

To load this list of isotopomers, first, create an instance of
the :ref:`simulator_api` class,

.. doctest::

    >>> from mrsimulator import Simulator

and then assign the list using either

.. doctest::

    >>> sim1 = Simulator(isotopomers)

or

.. doctest::

    >>> sim1 = Simulator()
    >>> sim1.isotopomers = isotopomers



**Import the list of isotopomers from JSON serialized file**

The list of isotopomers may directly be assigned to an instance of a
:ref:`simulator_api` class from a JSON serialized isotopomers file.
In the following example, we load an example
`JSON <https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers_ppm.json>`_
serialized isotopomers file. For this, we make use of the
meth:`~mrsimulator.Simulator.load_isotopomers` method as follows,


.. doctest::

    >>> filename = 'https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers_ppm.json'
    >>> st2 = Simulator()
    >>> st2.load_isotopomers(filename)
    Downloading '/DeepanshS/mrsimulator-test/master/isotopomers_ppm.json' from 'raw.githubusercontent.com' to file 'isotopomers_ppm_0.json'.
    [████████████████████████████████████████████████████████████████████]

.. testcleanup::

    import os
    os.remove('isotopomers_ppm.json')

The list of isotopomers from this file are

.. doctest::

    >>> from pprint import pprint
    >>> pprint(st2.isotopomers)
    [{'abundance': '100 %',
      'sites': [{'isotope_symbol': '13C',
                 'isotropic_chemical_shift': '1 ppm',
                 'shielding_symmetric': {'anisotropy': '-3.89 ppm',
                                         'asymmetry': 0.25}}]},
     {'sites': [{'isotope_symbol': '13C',
                 'isotropic_chemical_shift': '1 ppm',
                 'shielding_symmetric': {'anisotropy': '8.2 ppm',
                                         'asymmetry': 0.0}}]},
     {'sites': [{'isotope_symbol': '1H',
                 'isotropic_chemical_shift': '3 ppm',
                 'shielding_symmetric': {'anisotropy': '23.2 ppm',
                                         'asymmetry': 0.0}}]},
     {'sites': [{'isotope_symbol': '29Si',
                 'isotropic_chemical_shift': '-100 ppm',
                 'shielding_symmetric': {'anisotropy': '1.36 ppm',
                                         'asymmetry': 0.0}}]},
     {'sites': [{'isotope_symbol': '29Si',
                 'isotropic_chemical_shift': '-100 ppm',
                 'shielding_symmetric': {'anisotropy': '70.36 ppm',
                                         'asymmetry': 0.0}}]},
     {'sites': [{'isotope_symbol': '29Si',
                 'isotropic_chemical_shift': '-90 ppm',
                 'shielding_symmetric': {'anisotropy': '80.36 ppm',
                                         'asymmetry': 0.5}}]},
     {'sites': [{'isotope_symbol': '1H',
                 'isotropic_chemical_shift': '5.6 ppm',
                 'shielding_symmetric': {'anisotropy': '13.2 ppm',
                                         'asymmetry': 0.0}}]}]
