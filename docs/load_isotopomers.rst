

.. _load_isotopomers:


Importing isotopomers from JSON file
------------------------------------

A list of isotopomers may also be imported from a `JSON <https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers_ppm.json>`_
serialized isotopomers file. Consider `JSON isotopomers <https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers_test.json>`_ file.

To import the isotopomers from this file, use the
:meth:`~mrsimulator.Simulator.load_isotopomers`
method of the :ref:`simulator_api` class, as follows,

.. doctest::

    >>> from mrsimulator import Simulator
    >>> sim = Simulator()

    >>> filename = 'https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers_test.json'

    >>> sim.load_isotopomers(filename)
    Downloading '/DeepanshS/mrsimulator-test/master/isotopomers_test.json'
    from 'raw.githubusercontent.com' to file 'isotopomers_test.json'.
    [████████████████████████████████████]

.. testsetup::

    >>> import os
    >>> os.remove('isotopomers_test.json')
