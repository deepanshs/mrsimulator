

.. _load_spin_systems:

``mrsimulator`` I/O
===================

Simulator object
----------------

**Export simulator object to a JSON file**

To serialize a :ref:`simulator_api` object to a JSON-compliant file, use the
:py:meth:`~mrsimulator.Simulator.save` method of the object.

.. doctest::

    >>> sim_coesite.save('sample.mrsim')

where ``sim_coesite`` is a :ref:`simulator_api` object.
By default, the attribute values are serialized as physical quantities, represented
as a string with a value and a unit.

.. You may also serialize the file without the
.. units, in which case, follow

.. .. doctest::

..     >>> sim_coesite.save('sample_no_units.json', with_units=False)


**Load simulator object from a JSON file**

To load a JSON-compliant :ref:`simulator_api` serialized file, use the
:py:meth:`~mrsimulator.Simulator.load` method of the class. By default, the load method
parses the file for units.

.. doctest::

    >>> from mrsimulator import Simulator
    >>> sim_load = Simulator.load('sample.mrsim')
    >>> sim_coesite == sim_load
    True

.. If the file is serialized without the units, you may load the file as follows

.. .. doctest::

..     >>> sim_load_no_units = Simulator.load('sample_no_units.json', parse_units=False)
..     >>> sim_coesite == sim_load_no_units
..     True

.. testsetup::
    >>> import os
    >>> os.remove('sample.mrsim')


----


Spin systems objects from Simulator class
-----------------------------------------

**Export spin systems to a JSON file**

You may also serialize the spin system objects from the :ref:`simulator_api` object to
a JSON-compliant file using the :py:meth:`~mrsimulator.Simulator.export_spin_systems`
method as

.. doctest::

    >>> sim_coesite.export_spin_systems('coesite_spin_systems.mrsys')


**Import spin systems from a JSON file**

Similarly, a list of spin systems can be directly imported from a JSON serialized
file. To import the spin systems, use the
:py:meth:`~mrsimulator.Simulator.load_spin_systems` method of the :ref:`simulator_api`
class as

.. doctest::

    >>> sim.load_spin_systems('coesite_spin_systems.mrsys')

.. testsetup::
    >>> import os
    >>> os.remove('coesite_spin_systems.mrsys')

**Importing spin systems from URL**

.. doctest::

    >>> from mrsimulator import Simulator
    >>> sim = Simulator()
    >>> filename = 'https://raw.githubusercontent.com/deepanshs/mrsimulator-examples/master/spin_systems.json'
    >>> sim.load_spin_systems(filename)
    >>> # The seven spin systems from the file are added to the sim object.
    >>> len(sim.spin_systems)
    7

.. testsetup::
    >>> import os
    >>> os.remove('spin_systems.json')


Serialize simulation object from Method class as CSDM compliant file
--------------------------------------------------------------------

**Export simulation to a JSON file**

You may serialize the simulation object from the method object to a CSDM compliant JSON file using the
save function as follows,

.. doctest::

    >>> sim_coesite.method[0].simulation.save('coesite_simulation.csdf') # doctest:+SKIP


Serialize Simulator, SignalProcessor object to file
---------------------------------------------------

**Export Simulator, SignalProcessor objects to a JSON file**

You may serialize the Simulator, a list of SignalProcessor objects to a *.mrsim* file
as follows. The order of SignalProcessor objects is the order of the methods in the Simulator
object.

.. doctest::

    >>> from mrsimulator import save
    >>> save('coesite.mrsim', sim_coesite, processors) # doctest:+SKIP

**Load Simulator, SignalProcessor objects from a JSON file**

.. doctest::

    >>> from mrsimulator import load
    >>> sim_coesite, processors, _ = load('coesite.mrsim') # doctest:+SKIP
