

.. _dictionary_objects:

+++++++++++++++++++++++++++++++++++++++++
Creating instances from python dictionary
+++++++++++++++++++++++++++++++++++++++++

.. _orientation:

Orientation
^^^^^^^^^^^

We define `Orientation` as an object with attributes representing the
Euler angles---:math:`\alpha, \beta, \gamma`.


.. cssclass:: table-bordered table-hover
.. list-table::
  :widths: 15 20 65
  :header-rows: 1

  * - Attributes
    - Type
    - Description
  * - ``alpha``
    - ScalarQuantity
    - Euler angle, :math:`\alpha`. For example, ``2.1 rad``. The default value is
      ``0 rad``.
  * - ``beta``
    - ScalarQuantity
    - Euler angle, :math:`\beta`. For example, ``23.5 deg``. The default value is
      ``0 rad``.
  * - ``gamma``
    - ScalarQuantity
    - Euler angle, :math:`\gamma`. For example, ``0.5 rad``. The default value is
      ``0 rad``.

**Example**

Using python `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`_
object, Orientation may be represented as,

.. code-block:: py

    orientation_example = dict(alpha="0.5 rad", beta="0.23 rad", gamma="2.54 rad")


.. _symmetric_tensor:

SymmetricTensor
^^^^^^^^^^^^^^^

We define `SymmetricTensor` as an object with attributes representing the
traceless second rank symmetric spatial irreducible tensor.

.. cssclass:: table-bordered table-hover

.. list-table::
  :widths: 15 20 65
  :header-rows: 1

  * - Attributes
    - Type
    - Description
  * - ``zeta``
    - ScalarQuantity
    - The strength of the anisotropy calculated using Haeberlen
      convention. The value is a physical quantity given as either frequency,
      ``4.2 kHz``, or as a dimensionless frequency ratio, ``10 ppm`` or
      ``10 µHz/Hz``. The default value is ``0 ppm``.
  * - ``eta``
    - Float
    - The asymmetry parameter calculated using Haeberlen convention, for
      example, ``0.75``. The default value is ``0``.
  * - ``orientation``
    - :ref:`orientation`
    - The Euler angles (:math:`\alpha, \beta, \gamma`) that rotates the second rank
      irreducible symmetric tensor from the principal axis system (PAS) to the
      frame of reference using wigner rotations.

**Example**

Using python `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`_
object, SymmetricTensor may be represented as,

.. code-block:: py

    symmetric_tensor_example = dict(
        anisotropy="10.3 ppm", asymmetry=0.5, orientation=orientation_example
    )

where ``orientation_example`` is the dict object with Euler angles from the
previous example.





.. _site:

Site
^^^^

We define `Site` as an object with attributes representing an isolated nuclear
site.

.. cssclass:: table-bordered table-hover
.. list-table::
  :widths: 25 20 55
  :header-rows: 1

  * - Attributes
    - Type
    - Description
  * - ``isotope``
    - String
    - The isotope symbol of the NMR active nucleus, for example, ``13C``.
      This is a required key.
  * - ``isotropic_chemical_shift``
    - ScalarQuantity
    - The isotropic chemical shift of the site. The value is a physical quantity given as either frequency,
      ``4.2 kHz``, or as a dimensionless frequency ratio, ``10 ppm`` or
      ``10 µHz/Hz``. The default value is ``0 ppm``.
  * - ``shielding_symmetric``
    - :ref:`symmetric_tensor`
    - See the description for the :ref:`symmetric_tensor` object.

**Example**

Using python `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`_,
object, Site may be represented as,

.. code-block:: py

  site_example = dict(
      isotope="13C",
      isotropic_chemical_shift="15 ppm",
      shielding_symmetric=symmetric_tensor_example
      )
  )

where ``symmetric_tensor_example`` is the dict object with symmetric tensor
attributes from the previous example.


.. _isotopomer:

Isotopomer
^^^^^^^^^^

An `Isotopomer` object is a python
`dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`__
object which represents an isotopomer.
In `mrsimulator`, each `isotopomer` is treated as a :math:`n`-coupled spin
system where :math:`n` is the number of sites in the isotopomer.
It is recommended that if the sites are uncoupled, it be specified as
individual isotopomers with a single site object, rather than a single
isotopomer with multiple sites.

The key-value pairs of the `Isotopomer` object follow,

.. cssclass:: table-bordered table-hover
.. list-table::
  :widths: 15 15 70
  :header-rows: 1

  * - Attributes
    - Type
    - Description
  * - ``sites``
    - A `list <https://docs.python.org/3/library/stdtypes.html#list>`__
    - A list of :ref:`site` objects.
  * - ``abundance``
    - A `string <https://docs.python.org/3/library/functions.html#float>`__
    - The abundance of the isotopomer. The abundance is given as
      percent, for example, '25.4 %'. This key-value is useful when
      simulating multiple isotopomers. The default value is '100 %'.

..  * - ``coulpings``
..    - Not yet implemented.


*An example of Isotopomer object.*

.. code-block:: py

  isotopomer_object = dict(sites=[site_example], abundance="15.3 %")

where `site_example` is the dict object from the previous example.


.. _direct_dimension:

DirectDimension
^^^^^^^^^^^^^^^

A `DirectDimension` object is a python
`dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`__
object with the following key-value pairs.

.. cssclass:: table-bordered table-hover
.. list-table::
  :widths: 25 25 50
  :header-rows: 1

  * - Attributes
    - Type
    - Description
  * - ``isotope``
    - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
    - The isotope symbol of the nuclei. The recorded spectrum a histogram of
      frequencies corresponding to this nuclear isotope. An example may
      be '29Si'.
  * - ``magnetic_flux_density``
    - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
      containing a physical quantity
    - The strength of the external static magnetic field of the spectrometer,
      for example, '14.1 T'. The default value is '9.4 T'.
  * - ``rotor_frequency``
    - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
      containing a physical quantity
    - The sample spinning frequency, for example, '10 kHz'. The default value
      is '0 Hz'.
  * - ``rotor_angle``
    - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
      containing a physical quantity
    - The angle between the sample rotation axis and the external magnetic
      field, for example, ‘90 deg’. The default value is ‘54.735 deg’.
  * - ``number_of_points``
    - An `integer <https://docs.python.org/3.3/library/functions.html#int>`__
    - The number of points used in sampling the spectrum, for example, 8192.
      The default value is 1024.
  * - ``spectral_width``
    - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
      containing a physical quantity
    - The frequency spectral width over which the spectrum is evaluated,
      for example, '500 kHz'. The default value is '100 kHz'.
  * - ``reference_offset``
    - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
      containing a physical quantity
    - The reference offset of the spectrum, for example, '1.4587 kHz'.
      The default value is '0 Hz'.

.. Note::
    All physical quantities are specified as strings containing a numerical
    value and a unit.
