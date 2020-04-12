

.. _load_isotopomers:


Exporting/Importing isotopomers to/from JSON file
=================================================

Exporting isotopomers to a JSON file
------------------------------------

The isotopomers and dimensions from the :ref:`simulator_api` class object may
be serialized to a JSON compliant python dictionary object using the
:meth:`~mrsimulator.simulator.Simulator.to_dict_with_units` method. Consider
the simulator object from the coesite example.

.. doctest::

    >>> py_dict = sim_coesite.to_dict_with_units()
    >>> pprint(py_dict)
    {'isotopomers': [{'abundance': '0.83%',
                      'description': '',
                      'name': '',
                      'sites': [{'isotope': '17O',
                                 'isotropic_chemical_shift': '29.0 ppm',
                                 'quadrupolar': {'Cq': '6050000.0 Hz',
                                                 'eta': 0.0}}]},
                     {'abundance': '1.05%',
                      'description': '',
                      'name': '',
                      'sites': [{'isotope': '17O',
                                 'isotropic_chemical_shift': '41.0 ppm',
                                 'quadrupolar': {'Cq': '5430000.0 Hz',
                                                 'eta': 0.166}}]},
                     {'abundance': '2.16%',
                      'description': '',
                      'name': '',
                      'sites': [{'isotope': '17O',
                                 'isotropic_chemical_shift': '57.0 ppm',
                                 'quadrupolar': {'Cq': '5450000.0 Hz',
                                                 'eta': 0.168}}]},
                     {'abundance': '2.05%',
                      'description': '',
                      'name': '',
                      'sites': [{'isotope': '17O',
                                 'isotropic_chemical_shift': '53.0 ppm',
                                 'quadrupolar': {'Cq': '5520000.0 Hz',
                                                 'eta': 0.169}}]},
                     {'abundance': '1.9%',
                      'description': '',
                      'name': '',
                      'sites': [{'isotope': '17O',
                                 'isotropic_chemical_shift': '58.0 ppm',
                                 'quadrupolar': {'Cq': '5160000.0 Hz',
                                                 'eta': 0.292}}]}]}

The :meth:`~mrsimulator.simulator.Simulator.to_dict_with_units` method returns
a python dictionary, which can then be serialized to a file using the JSON
module, as follows

.. doctest::

    >>> import json
    >>> filename = 'my_serialized_file.json'
    >>> with open(filename, 'w') as f:
    ...     json.dump(py_dict, f)

.. testsetup::
    >>> import os
    >>> os.remove('my_serialized_file.json')

Importing isotopomers from a JSON file
--------------------------------------

Similarly, a list of isotopomers can be directly imported from a JSON
serialized isotopomers file. Consider this `JSON isotopomers <https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers_test.json>`_ file.
To import the isotopomers from this file, use the
:meth:`~mrsimulator.simulator.Simulator.load_isotopomers`
method of the :ref:`simulator_api` class, as follows

.. doctest::

    >>> from mrsimulator import Simulator
    >>> sim = Simulator()

    >>> filename = 'https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers_test.json'

    >>> sim.load_isotopomers(filename)  # doctest:+SKIP
    Downloading '/DeepanshS/mrsimulator-test/master/isotopomers_test.json'
    from 'raw.githubusercontent.com' to file 'isotopomers_test.json'.
    [████████████████████████████████████]

    >>> # The seven isotopomers from the file are added to the isotopomers
    >>> # attribute of the simulator class.
    >>> len(sim.isotopomers) # doctest:+SKIP
    7

.. .. testsetup::
..     >>> import os
..     >>> os.remove('isotopomers_test.json')
