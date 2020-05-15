

.. _load_isotopomers:


Exporting/Importing isotopomers
===============================

Exporting isotopomers to a JSON file
------------------------------------

The isotopomers from the :ref:`simulator_api` class object may be serialized to a JSON
file using the :meth:`~mrsimulator.Simulator.export_isotopomers` method. Consider
the simulator object from the coesite example.

.. doctest::

    >>> py_dict = sim_coesite.export_isotopomers('coesite_isotopomer.json')

.. testsetup::
    >>> import os
    >>> os.remove('coesite_isotopomer.json')

Importing isotopomers from a JSON file
--------------------------------------

Similarly, a list of isotopomers can be directly imported from a JSON
serialized isotopomers file. Consider this `JSON isotopomers <https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers_test.json>`_ file.
To import the isotopomers from this file, use the
:meth:`~mrsimulator.Simulator.load_isotopomers`
method of the :ref:`simulator_api` class, as follows

.. testsetup::
    >>> from mrsimulator import Simulator
    >>> sim = Simulator()
    >>> filename = 'https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers_test.json'
    >>> sim.load_isotopomers(filename)
    Downloading '/DeepanshS/mrsimulator-test/master/isotopomers_test.json'
    from 'raw.githubusercontent.com' to file 'isotopomers_test.json'.
    [████████████████████████████████████]

.. doctest::

    >>> from mrsimulator import Simulator # doctest:+SKIP
    >>> sim = Simulator() # doctest:+SKIP

    >>> filename = 'https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers_test.json'

    >>> sim.load_isotopomers(filename) # doctest:+SKIP
    Downloading '/DeepanshS/mrsimulator-test/master/isotopomers_test.json'
    from 'raw.githubusercontent.com' to file 'isotopomers_test.json'.
    [████████████████████████████████████]

    >>> # The seven isotopomers from the file are added to the isotopomers
    >>> # attribute of the simulator class.
    >>> len(sim.isotopomers) # doctest:+SKIP
    7

.. testsetup::
    >>> import os
    >>> os.remove('isotopomers_test.json')
