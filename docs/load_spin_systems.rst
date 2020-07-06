

.. _load_spin_systems:


Exporting/Importing spin systems
================================

**Exporting spin systems to a JSON file**

The spin systems from the :ref:`simulator_api` class object may be serialized to a JSON
compliant file-format using the :meth:`~mrsimulator.Simulator.export_spin_systems`
method with the following syntax,

.. doctest::

    >>> sim_coesite.export_spin_systems('coesite_sample.json')

.. testsetup::
    >>> import os
    >>> os.remove('coesite_sample.json')

**Importing spin systems from a JSON file**

Similarly, a list of spin systems can be directly imported from a JSON serialized
file. To import the spin systems, use the
:meth:`~mrsimulator.Simulator.load_spin_systems` method of the :ref:`simulator_api`
class with the syntax

.. code-block:: python

    sim.load_spin_systems(filename)

Here is an example.

.. testsetup::
    >>> from mrsimulator import Simulator
    >>> sim = Simulator()
    >>> filename = 'https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/spin_systems_v0.3.json'
    >>> sim.load_spin_systems(filename)
    Downloading '/DeepanshS/mrsimulator-test/master/spin_systems_v0.3.json'
    from 'raw.githubusercontent.com' to file 'spin_systems_v0.3.json'.
    [███████████████████████]

.. doctest::

    >>> from mrsimulator import Simulator # doctest:+SKIP
    >>> sim = Simulator() # doctest:+SKIP

    >>> filename = 'https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/spin_systems_v0.3.json'

    >>> sim.load_spin_systems(filename) # doctest:+SKIP
    Downloading '/DeepanshS/mrsimulator-test/master/spin_systems_v0.3.json'
    from 'raw.githubusercontent.com' to file 'spin_systems_v0.3.json'.
    [███████████████████████]

    >>> # The seven spin systems from the file are added to the sim object.
    >>> len(sim.spin_systems)
    7

.. testsetup::
    >>> import os
    >>> os.remove('spin_systems_v0.3.json')
