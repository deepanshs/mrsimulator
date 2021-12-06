.. _simulator_documentation:

====================
The Simulator Object
====================

The :ref:`simulator_api` object is the core of the ``mrsimulator`` library. Each **Simulator**
object holds a list of :ref:`spin_sys_api` objects and a list of :ref:`method_api` objects.
For a complete list of attributes, see ((table ref))

Setting up a Simulator Object
-----------------------------

Setting up a simulator object and running a simulation is simple. Below we add some arbitrary
spin systems and methods to a simulator object.

.. code-block:: python

    from mrsimulator import Site, Simulator, spin_system
    from mrsimulator.methods import BlochDecaySpectrum

    # Setup the spin system and method objects
    system1 = SpinSystem(sites=[Site(isotope="1H")])  # Proton spin system
    system2 = SpinSystem(sites=[Site(isotope="17O")])  # Oxygen spin system
    system3 = SpinSystem(sites=[Site(isotope="29Si")])  # Silicon spin system
    method1 = BlochDecaySpectrum(channels=["1H"])
    method2 = BlochDecaySpectrum(channels=["29Si"])

    # Create the Simulator object
    sim = Simulator()
    sim.spin_systems = [system1, system2, system3]  # Add list of spin systems
    sim.methods = [method1, method2]  # add list of methods

``sim`` is a **Simulator** object which holds three spin systems and two methods. See
:ref:`spin_system_documentation` and :ref:`method_documentation` documentation for more
information on the respective classes.

Running a Simulation
--------------------

To simulate the NMR spectrum of the given spin systems using each method, call the simulator
class method :meth:`~mrsimulator.Simulator.run`.

.. code-block:: python

    sim.run()

The simulated spectrum is stored as a CSDM object in each method object under the
``simulation`` attribute. For more information on the Core Scientific Data Model (CSDM),
see the `csdmpy documentation <https://csdmpy.readthedocs.io/en/stable/>`_.
Below we put the simulated spectra of the method at index 0 into the variable ``data_0``

.. code-block:: python

    data_0 = sim.methods[0].simulation
    # data_n = sim.methods[n].simulation (for multiple methods)

Specifying which methods to simulate
''''''''''''''''''''''''''''''''''''

By default, :meth:`~mrsimulator.Simulator.run` simulates the spectrum of the given spin systems
over all methods. You may specify which methods to simulate using the ``method_index`` argument.
``method_index`` accepts a list of integers specifying the index of methods to simulate. The code
below simulates the first and third methods in ``sim``

.. code-block:: python

    sim.run(method_index=[0, 2])

Packing the data as Numpy array
'''''''''''''''''''''''''''''''

By default the simulated spectrum is packed into a CSDM object. The spectrum can also be packed
as a numpy array by using the ``pack_as_csdm`` argument.

.. code-block:: python

    sim.run(pack_as_csdm=False)

Although this packing the simulated spectrum as a numpy array is possible,
**we strongly recommend against it since this breaks serialization**.

Sub-class Documentation
-----------------------

The following sections cover the usage of **Simulator** subclasses as well as
configuring properties of the simulation and serialization.

.. toctree::
    :maxdepth: 2
    :caption: Sections

    spin_system
    method
    configuring_simulator
