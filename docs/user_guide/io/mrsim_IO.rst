.. _IO_documentation:

===============
mrsimulator I/O
===============

We offer a range of serialization options based on a JSON structure demonstrated below.

Dictionary Representation of Objects
------------------------------------

All **mrsimulator** objects can be serialized into a JSON format. Calling the
``json()`` method on an object will return a Python dictionary representing the object
in JSON format.
Below we call the :meth:`~mrsimulator.Site.json` method of the :ref:`site_api` class.

.. code-block:: python

    from mrsimulator import Site, SpinSystem
    from mrsimulator.spin_system.tensors import SymmetricTensor

    Si29_site = Site(
        isotope="29Si",
        isotropic_chemical_shift=-89.0,
        shielding_symmetric=SymmetricTensor(
            zeta=59.8,
            eta=0.62,
        ),
    )

    py_dict = Si29_site.json()
    print(py_dict)
    # {
    #     'isotope': '29Si',
    #     'isotropic_chemical_shift': '-89.0 ppm',
    #     'shielding_symmetric': {'zeta': '59.8 ppm', 'eta': 0.62}
    # }

All values are serialized with units when applicable, but you may call ``json(units=False)``
if you wish to serialize values without units.

Similarly, all **mrsimulator** objects can be loaded from a dictionary representation. Here we
construct the same site as a dictionary and call :meth:`~mrsimulator.Site.parse_dict_with_units`
to create a :ref:`site_api` object from a dictionary.

.. code-block:: python

    site_dict = {
        "isotope": "29Si",
        "isotropic_chemical_shift": "-89.0 ppm",
        "shielding_symmetric": {
            "zeta": "59.8 ppm",
            "eta": 0.62,
        },
    }

    Si29_site_from_dict = Site().parse_dict_with_units(site_dict)
    print(Si29_site_from_dict == Si29_site)
    # True

We see that both these sites are equivalent. Values in dictionaries can be given as a
number and a unit in a string. However, passing values with units increases overhead and
will throw errors if the units cannot be converted into the expected units for a field.
For this reason, we recommend instantiating objects directly from classes.

.. _load_spin_systems:

Saving and Loading Spin Systems from a File
-------------------------------------------

A list of spin systems in a :ref:`simulator_api` object can be serialized to a file. Here we create
a simulator with three distinct :math:`^{29}\text{Si}` spin systems and serialize these spin
systems to a file by calling :meth:`~mrsimulator.Simulator.export_spin_systems`.

.. code-block:: python

    from mrsimulator import Site, SpinSystem, Simulator
    from mrsimulator.spin_system.tensors import SymmetricTensor

    # Create the spin systems
    Si29_1 = SpinSystem(
        sites=[
            Site(
                isotope="29Si",
                isotropic_chemical_shift=-89.0,
                shielding_symmetric=SymmetricTensor(zeta=59.8, eta=0.62),
            )
        ]
    )
    Si29_2 = SpinSystem(
        sites=[
            Site(
                isotope="29Si",
                isotropic_chemical_shift=-89.5,
                shielding_symmetric=SymmetricTensor(zeta=52.1, eta=0.68),
            )
        ]
    )
    Si29_3 = SpinSystem(
        sites=[
            Site(
                isotope="29Si",
                isotropic_chemical_shift=-87.8,
                shielding_symmetric=SymmetricTensor(zeta=69.4, eta=0.60),
            )
        ]
    )

    # Create the Simulator object
    sim = Simulator(spin_systems=[Si29_1, Si29_2, Si29_3])

    # Save spin systems to file
    sim.export_spin_systems("example.mrsys")

Now the file ``example.mrsys`` holds a JSON representation of the spin system objects. The
extension of the file is irrelevant; however, we strongly encourage using ``.mrsys`` to
adhere to the convention.

Just as spin systems can be saved to a file, spin systems can be loaded from a file. Loading spin
systems is useful when working with a large number of spin systems over multiple Python scripts. Here
we load the spin system file, ``example.mrsys``, into a new simulator using the method
:meth:`~mrsimulator.Simulator.load_spin_systems`.

.. code-block:: python

    new_sim = Simulator()
    new_sim.load_spin_systems("example.mrsys")
    print(len(new_sim.spin_systems))
    # 3

Saving and Loading Methods from a File
--------------------------------------

A list of methods in a :ref:`simulator_api` object can be serialized to a file. Here we create a
custom DAS method and serialize it to a file using the method
:meth:`~mrsimulator.Simulator.export_methods`.

.. code-block:: python

    from mrsimulator import Simulator
    from mrsimulator.method import Method
    from mrsimulator.method import SpectralDimension, SpectralEvent

    # Create DAS method
    das = Method(
        name="DAS of 17O",
        channels=["17O"],
        magnetic_flux_density=11.744,
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=10000,
                reference_offset=-1220.9,
                origin_offset=67793215,
                label="Isotropic dimension",
                events=[
                    SpectralEvent(
                        fraction=0.5,
                        rotor_angle=37.38 * 3.14159 / 180,
                        transition_query=[{"ch1": {"P": [-1], "D": [0]}}],
                    ),
                    SpectralEvent(
                        fraction=0.5,
                        rotor_angle=79.19 * 3.14159 / 180,
                        transition_query=[{"ch1": {"P": [-1], "D": [0]}}],
                    ),
                ],
            ),
            # The last spectral dimension block is the direct-dimension
            SpectralDimension(
                count=256,
                spectral_width=11001,
                reference_offset=-1228,
                origin_offset=67793215,
                label="MAS dimension",
                events=[
                    SpectralEvent(
                        rotor_angle=54.735 * 3.14159 / 180,
                        transition_query=[{"ch1": {"P": [-1], "D": [0]}}],
                    )
                ],
            ),
        ],
    )

    # Create simulator with das method
    sim = Simulator(methods=[das])

    # Save methods to file
    sim.export_methods("example.mrmtd")

Now the file ``example.mrmtd`` holds a JSON representation of the method object. If multiple
methods are present, e.g., at different spinning speeds, they will also be serialized. The file's extension
is not essential; however, we strongly encourage using ``.mrmtd`` to adhere to the convention.

Just like spin systems, methods can also be loaded from a file. Here we load the DAS
method into a new simulator object by calling the method
:meth:`~mrsimulator.Simulator.load_methods`.

.. code-block:: python

    new_sim = Simulator()
    new_sim.load_methods("example.mrmtd")
    print(new_sim.methods[0].name)
    # DAS of 17O

Loading complex methods from a file, like the DAS example above, can reduce complex code.
Methods representing actual experiments can be saved to a file to later be loaded into a script
as needed.

Serializing a Simulator Object
------------------------------

The entire :ref:`simulator_api` object may be serialized to a JSON-compliant file using the
:meth:`~mrsimulator.Simulator.save` Python method.
By default, the attribute values are serialized as physical quantities represented as a
string with a value and a unit.

.. code-block:: python

    sim = Simulator()
    # ... Setup Simulator object
    sim.save("sample.mrsim")

Now the file ``sample.mrsim`` holds the JSON representation of ``sim``, a :ref:`simulator_api` object.
To load a simulator from a file, call the class method :meth:`~mrsimulator.Simulator.load`.
By default, the load method parses the file for units.

.. code-block:: python

    new_sim = Simulator.load("sample.mrsim")

Serialize simulation from a Method to a CSDM Compliant File
-----------------------------------------------------------

The simulated spectrum may be exported to a CSDM-compliant JSON file using the following code:

.. skip: next
.. code-block:: python

    sim_coesite.methods[0].simulation.save("coesite_simulation.csdf")


For more information on the CSDM format, see the
`csdmpy documentation <https://csdmpy.readthedocs.io/en/stable/>`__.

Serialize Simulator and SignalProcessor object
----------------------------------------------

The :ref:`simulator_api` object and a list of :ref:`signal_processor_api` objects
can both be serialized within the same file by calling the :meth:`~mrsimulator.save`
method.

.. code-block:: python

    from mrsimulator import save
    from mrsimulator import Simulator
    from mrsimulator import signal_processor as sp

    sim = Simulator()
    processor1 = sp.SignalProcessor()
    processor2 = sp.SignalProcessor()

    save(
        filename="example.mrsim",
        simulator=sim,
        signal_processors=[processor1, processor2],
    )

All attribute values are serialized with units by default, but you may serialize without
units by passing ``with_units=False`` to the method.
Additionally, a metadata dictionary can be passed using the ``application`` keyword.

To load a simulator and signal processors from a file, call the :meth:`~mrsimulator.load`
method. This method will return an ordered list of a :ref:`simulator_api` object, a list of
:ref:`signal_processor_api` objects, and a metadata dictionary

.. code-block:: python

    from mrsimulator import load

    sim, processors, application = load("example.mrsim")

.. note::

    The serialization structure has been updated in **mrsimulator** v0.7. **Mrsimulator** should automatically
    take care of this update when loading files from v0.6 and below. However, you can use the
    :py:meth:`~mrsimulator.update_old_file_struct` method
    to convert older files to the new format.


.. plot::
    :include-source: False

    import os
    from os.path import isfile

    if isfile("example.mrmtd"): os.remove("example.mrmtd")
    if isfile("example.mrsim"): os.remove("example.mrsim")
    if isfile("example.mrsys"): os.remove("example.mrsys")
    if isfile("sample.mrsim"): os.remove("sample.mrsim")
